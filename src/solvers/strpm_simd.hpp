/*
 * Copyright 2017-2018 Tom van Dijk, Johannes Kepler University Linz
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef STRPM_SIMD_HPP
#define STRPM_SIMD_HPP

#include "oink/solver.hpp"
#include <experimental/simd>
#include <cstring>
#include <cstddef>

namespace stdx = std::experimental;
using simd_uint8      = stdx::fixed_size_simd<uint8_t, 8>;
using simd_uint8_mask = stdx::fixed_size_simd_mask<uint8_t, 8>;
using simd_int16      = stdx::fixed_size_simd<int16_t, 8>;
using simd_int16_mask = stdx::fixed_size_simd_mask<int16_t, 8>;

// Top is represented by levels[0] == -1 (matching the scalar strpm solver).

// Bit-parallel popcount for all 8 uint8 lanes simultaneously.
// Uses the standard Hamming-weight (sideways addition) algorithm:
//   Step 1: Pair adjacent bits: count = b1+b0 for each 2-bit group
//           x - ((x >> 1) & 0x55) works because for 2-bit value ab:
//           ab - 0b = ab (i.e. a+b when both ≤1), handles carry correctly.
//   Step 2: Sum adjacent pairs into 4-bit nibble counts via masking.
//   Step 3: Sum adjacent nibbles; result fits in low nibble, mask off high.
inline simd_uint8 simd_popcount8(simd_uint8 x) noexcept {
    x = x - ((x >> 1) & simd_uint8(0x55));          // 2-bit sums
    x = (x & simd_uint8(0x33)) + ((x >> 2) & simd_uint8(0x33)); // 4-bit sums
    return (x + (x >> 4)) & simd_uint8(0x0F);        // 8-bit sum (0..8)
}

// Precomputed lane index vectors [0,1,2,...,7] for uint8 and int16.
static const simd_uint8 LANE_INDICES{[](uint8_t i){ return i; }};
static const simd_int16 LANE_INDICES_I16{[](int16_t i){ return i; }};

// Per-node progress-measure record.
// Packing all four fields into one struct collapses the four separate SoA
// vectors into a single contiguous allocation: one cache-line fetch (≤2 in
// the worst case) per to_tmp()/from_tmp() call instead of four.
//
// Layout (40 bytes, naturally 2-byte aligned):
//   offset  0 :  bits[8]   — uint8_t, one bit-string value per SIMD lane
//   offset  8 :  masks[8]  — uint8_t, which bit positions are occupied
//   offset 16 :  levels[8] — int16_t, tree depth per lane (−1 == Top sentinel)
//   offset 32 :  nlanes    — uint8_t, number of active lanes (0 … k-1)
//   offset 33 :  _pad[7]   — padding → sizeof == 40
struct NodePM {
    uint8_t  bits[8];    // SIMD-ready: 8 contiguous uint8 → copy_from/copy_to
    uint8_t  masks[8];   // SIMD-ready: 8 contiguous uint8
    int16_t  levels[8];  // 8 contiguous int16_t; natural 2-byte alignment at offset 16
    uint8_t  nlanes;
    uint8_t  _pad[7];    // bring sizeof to 40
};
static_assert(sizeof(NodePM) == 40, "NodePM must be 40 bytes");
static_assert(offsetof(NodePM, bits)   ==  0);
static_assert(offsetof(NodePM, masks)  ==  8);
static_assert(offsetof(NodePM, levels) == 16);
static_assert(offsetof(NodePM, nlanes) == 32);

namespace pg {

class STRPM_SIMDSolver : public Solver
{
public:
    STRPM_SIMDSolver(Oink& oink, Game& game);
    virtual ~STRPM_SIMDSolver();

    virtual void run();

protected:
    /**
     * Parameters: U^k_{t, h}
     *      - k: Strahler-number
     *      - t: number of bits
     *      - h: height
     */
    int k, t, h;

    // AoS layout: one NodePM per node, all four fields contiguous.
    // Each to_tmp()/from_tmp() call touches a single 40-byte region (≤2 cache
    // lines) instead of four separate vectors at four different addresses.
    // Levels stored as int16_t (values in [0, h-2], sentinel -1 for Top).
    // h is asserted < 32767 at startup so int16_t always suffices.
    std::vector<NodePM>   pm;

    simd_uint8 tmp_bits;
    simd_uint8 tmp_masks;
    alignas(16) int16_t tmp_levels[8]; // alignas(16) for simd copy_from/copy_to
    uint8_t tmp_nlanes;

    simd_uint8 best_bits;
    simd_uint8 best_masks;
    alignas(16) int16_t best_levels[8];
    uint8_t best_nlanes;

    uintqueue Q;
    bitset dirty;

    bool always_reset = false;

    uint64_t *lift_counters;

    // Reusable buffer for collecting valid successors in lift(),
    // pre-reserved to nodecount() to avoid per-call heap allocation.
    std::vector<int> succs;

    // Copy pm[idx] into tmp
    inline void to_tmp(int idx) {
        tmp_bits.copy_from(&pm[idx].bits[0], stdx::element_aligned);
        tmp_masks.copy_from(&pm[idx].masks[0], stdx::element_aligned);
        std::memcpy(tmp_levels, &pm[idx].levels[0], 8 * sizeof(int16_t));
        tmp_nlanes = pm[idx].nlanes;
    }
    // Copy tmp into pm[idx]
    inline void from_tmp(int idx) {
        tmp_bits.copy_to(&pm[idx].bits[0], stdx::element_aligned);
        tmp_masks.copy_to(&pm[idx].masks[0], stdx::element_aligned);
        std::memcpy(&pm[idx].levels[0], tmp_levels, 8 * sizeof(int16_t));
        pm[idx].nlanes = tmp_nlanes;
    }
    // Copy pm[idx] into best
    inline void to_best(int idx) {
        best_bits.copy_from(&pm[idx].bits[0], stdx::element_aligned);
        best_masks.copy_from(&pm[idx].masks[0], stdx::element_aligned);
        std::memcpy(best_levels, &pm[idx].levels[0], 8 * sizeof(int16_t));
        best_nlanes = pm[idx].nlanes;
    }
    // Copy best into pm[idx]
    inline void from_best(int idx) {
        best_bits.copy_to(&pm[idx].bits[0], stdx::element_aligned);
        best_masks.copy_to(&pm[idx].masks[0], stdx::element_aligned);
        std::memcpy(&pm[idx].levels[0], best_levels, 8 * sizeof(int16_t));
        pm[idx].nlanes = best_nlanes;
    }
    // Copy tmp into best
    inline void tmp_to_best() {
        best_bits = tmp_bits;
        best_masks = tmp_masks;
        std::memcpy(best_levels, tmp_levels, 8 * sizeof(int16_t));
        best_nlanes = tmp_nlanes;
    }

    // Zero out inactive lanes' bits and masks (levels are scalar, not touched)
    inline void fill_inactive_tmp() {
        simd_uint8_mask inactive = LANE_INDICES >= simd_uint8(tmp_nlanes);
        stdx::where(inactive, tmp_bits) = simd_uint8(0);
        stdx::where(inactive, tmp_masks) = simd_uint8(0);
    }

    // Render pm[idx] to given ostream
    void stream_pm(std::ostream &out, int idx);
    // Render SIMD to given ostream
    void stream_simd(std::ostream &out, simd_uint8& bits, simd_uint8& masks, int16_t* levels, uint8_t nlanes);

    // Compare tmp to best
    int compare(int pindex);

    void prog_tmp(int pindex, int h);

    // Lift node, triggered by change to target
    bool lift(int node, int target, int &str, int pl);

    inline void todo_push(int node) {
        if (dirty[node]) return;
        Q.push(node);
        dirty[node] = true;
#ifndef NDEBUG
        if (trace >= 2) logger << "push(" << node << ")" << std::endl;
#endif
    }

    inline int todo_pop() {
        int node = Q.pop();
        dirty[node] = false;
#ifndef NDEBUG
        if (trace >= 2) logger << "pop() => " << node << std::endl;
#endif
        return node;
    }

    int lift_count = 0;
    int lift_attempt = 0;

    void run(int t_val, int k_val, int depth, int player);
};

}

#endif
