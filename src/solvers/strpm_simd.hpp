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

// Architecture detection for native SIMD prefix sum
#if defined(__SSE2__)
#include <emmintrin.h>
#define STRPM_SIMD_SSE2 1
#elif defined(__ARM_NEON) || defined(__aarch64__)
#include <arm_neon.h>
#define STRPM_SIMD_NEON 1
#endif

namespace stdx = std::experimental;

// Unified uint16 SIMD type: bits, masks, and levels all share the same lane
// width.  This eliminates the costly int16→uint8 narrowing conversion that was
// previously required to blend level predicates with bitstring operations.
// Top sentinel: 0xFFFF (was -1 as int16_t).
using simd_u16x8      = stdx::fixed_size_simd<uint16_t, 8>;
using simd_u16x8_mask = stdx::fixed_size_simd_mask<uint16_t, 8>;

// Bit-parallel popcount for all 8 uint16 lanes simultaneously.
// Bitstring values never exceed 8 bits, so upper bytes are always zero;
// the extra reduction step (vs the old uint8 version) is one SIMD op.
inline simd_u16x8 simd_popcount_u16x8(simd_u16x8 x) noexcept {
    x = x - ((x >> 1) & simd_u16x8(0x5555));                          // 2-bit sums
    x = (x & simd_u16x8(0x3333)) + ((x >> 2) & simd_u16x8(0x3333));  // 4-bit sums
    x = (x + (x >> 4)) & simd_u16x8(0x0F0F);                          // 8-bit sums
    return (x + (x >> 8)) & simd_u16x8(0x00FF);                       // 16-bit sum (0..16)
}

// Inclusive prefix sum over 8 uint16 lanes: result[i] = sum(x[0..i]).
// Three-step parallel scan (log2(8) dependent adds) vs 7 sequential adds.
// SSE2 and NEON paths use native intrinsics; fallback uses scalar element access.
inline simd_u16x8 simd_prefix_sum_inclusive_u16x8(simd_u16x8 x) noexcept {
#if defined(STRPM_SIMD_SSE2)
    // Marshal simd_u16x8 → __m128i via aligned buffer (compiler elides the memops).
    alignas(16) uint16_t buf[8];
    x.copy_to(buf, stdx::vector_aligned);
    __m128i v = _mm_load_si128(reinterpret_cast<const __m128i*>(buf));
    // Step 1: add neighbor — shift left by 1 lane (2 bytes), add
    v = _mm_add_epi16(v, _mm_slli_si128(v, 2));
    // Step 2: add pair — shift left by 2 lanes (4 bytes), add
    v = _mm_add_epi16(v, _mm_slli_si128(v, 4));
    // Step 3: add quad — shift left by 4 lanes (8 bytes), add
    v = _mm_add_epi16(v, _mm_slli_si128(v, 8));
    // Marshal __m128i → simd_u16x8
    _mm_store_si128(reinterpret_cast<__m128i*>(buf), v);
    simd_u16x8 result;
    result.copy_from(buf, stdx::vector_aligned);
    return result;
#elif defined(STRPM_SIMD_NEON)
    // Marshal simd_u16x8 → uint16x8_t via aligned buffer (compiler elides the memops).
    alignas(16) uint16_t buf[8];
    x.copy_to(buf, stdx::vector_aligned);
    uint16x8_t v = vld1q_u16(buf);
    uint16x8_t zero = vdupq_n_u16(0);
    // Step 1: shift left by 1 lane, add
    v = vaddq_u16(v, vextq_u16(zero, v, 7));
    // Step 2: shift left by 2 lanes, add
    v = vaddq_u16(v, vextq_u16(zero, v, 6));
    // Step 3: shift left by 4 lanes, add
    v = vaddq_u16(v, vextq_u16(zero, v, 4));
    // Marshal uint16x8_t → simd_u16x8
    vst1q_u16(buf, v);
    simd_u16x8 result;
    result.copy_from(buf, stdx::vector_aligned);
    return result;
#else
    // Fallback: sequential prefix sum via std::experimental::simd element access
    simd_u16x8 result;
    result[0] = x[0];
    for (int i = 1; i < 8; i++)
        result[i] = result[i-1] + x[i];
    return result;
#endif
}

// Precomputed lane index vector [0,1,2,...,7] for uint16.
static const simd_u16x8 LANE_INDICES{[](uint16_t i){ return i; }};

static constexpr uint16_t TOP_SENTINEL = 0xFFFF;

// Per-node progress-measure record.
// Packing all four fields into one struct collapses the four separate SoA
// vectors into a single contiguous allocation: one cache-line fetch
// per to_tmp()/from_tmp() call instead of four.
//
// Layout (56 bytes, fits in one 64-byte cache line):
//   offset  0 :  bits[8]   — uint16_t, one bit-string value per SIMD lane
//   offset 16 :  masks[8]  — uint16_t, which bit positions are occupied
//   offset 32 :  levels[8] — uint16_t, tree depth per lane (0xFFFF == Top)
//   offset 48 :  nlanes    — uint8_t, number of active lanes (0 … k-1)
//   offset 49 :  _pad[7]   — padding → sizeof == 56
struct NodePM {
    uint16_t bits[8];    // SIMD-ready: 8 contiguous uint16 → copy_from/copy_to
    uint16_t masks[8];   // SIMD-ready: 8 contiguous uint16
    uint16_t levels[8];  // 8 contiguous uint16_t; 0xFFFF == Top sentinel
    uint8_t  nlanes;
    uint8_t  _pad[7];    // bring sizeof to 56
};
static_assert(sizeof(NodePM) == 56, "NodePM must be 56 bytes");
static_assert(offsetof(NodePM, bits)   ==  0);
static_assert(offsetof(NodePM, masks)  == 16);
static_assert(offsetof(NodePM, levels) == 32);
static_assert(offsetof(NodePM, nlanes) == 48);

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
    // Each to_tmp()/from_tmp() call touches a single 56-byte region (one cache
    // line) instead of four separate vectors at four different addresses.
    // Levels stored as uint16_t (values in [0, h-2], sentinel 0xFFFF for Top).
    // h is asserted < 32767 at startup so uint16_t always suffices.
    std::vector<NodePM>   pm;

    simd_u16x8 tmp_bits;
    simd_u16x8 tmp_masks;
    simd_u16x8 tmp_levels;
    uint8_t tmp_nlanes;

    simd_u16x8 best_bits;
    simd_u16x8 best_masks;
    simd_u16x8 best_levels;
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
        tmp_levels.copy_from(&pm[idx].levels[0], stdx::element_aligned);
        tmp_nlanes = pm[idx].nlanes;
    }
    // Copy tmp into pm[idx]
    inline void from_tmp(int idx) {
        tmp_bits.copy_to(&pm[idx].bits[0], stdx::element_aligned);
        tmp_masks.copy_to(&pm[idx].masks[0], stdx::element_aligned);
        tmp_levels.copy_to(&pm[idx].levels[0], stdx::element_aligned);
        pm[idx].nlanes = tmp_nlanes;
    }
    // Copy pm[idx] into best
    inline void to_best(int idx) {
        best_bits.copy_from(&pm[idx].bits[0], stdx::element_aligned);
        best_masks.copy_from(&pm[idx].masks[0], stdx::element_aligned);
        best_levels.copy_from(&pm[idx].levels[0], stdx::element_aligned);
        best_nlanes = pm[idx].nlanes;
    }
    // Copy best into pm[idx]
    inline void from_best(int idx) {
        best_bits.copy_to(&pm[idx].bits[0], stdx::element_aligned);
        best_masks.copy_to(&pm[idx].masks[0], stdx::element_aligned);
        best_levels.copy_to(&pm[idx].levels[0], stdx::element_aligned);
        pm[idx].nlanes = best_nlanes;
    }
    // Copy tmp into best
    inline void tmp_to_best() {
        best_bits = tmp_bits;
        best_masks = tmp_masks;
        best_levels = tmp_levels;
        best_nlanes = tmp_nlanes;
    }

    // Zero out inactive lanes' bits and masks
    inline void fill_inactive_tmp() {
        simd_u16x8_mask inactive = LANE_INDICES >= simd_u16x8(tmp_nlanes);
        stdx::where(inactive, tmp_bits) = simd_u16x8(0);
        stdx::where(inactive, tmp_masks) = simd_u16x8(0);
    }

    // Render pm[idx] to given ostream
    void stream_pm(std::ostream &out, int idx);
    // Render SIMD to given ostream
    void stream_simd(std::ostream &out, simd_u16x8& bits, simd_u16x8& masks, simd_u16x8& levels, uint8_t nlanes);

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
