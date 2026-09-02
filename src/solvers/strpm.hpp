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

#ifndef STRPM_HPP
#define STRPM_HPP

#include "oink/solver.hpp"

#include <vector>
#include <cstdint>
#include <boost/container/small_vector.hpp>

namespace pg {

// Progress-measure storage. The level array is non-decreasing and the whole
// measure length is bounded by ~k-1+t, which stays within STRPM_INLINE for the
// SIMD-representable range and for most practical games. Small-buffer
// optimization keeps per-node measures and the tmp/best/test working copies
// allocation-free on the common (small) case; larger measures spill to the heap
// transparently. This removes the per-node heap allocations and whole-vector
// copies that dominated the scalar solver's cost.
inline constexpr std::size_t STRPM_INLINE = 16;
using BitVec = boost::container::small_vector<uint8_t, STRPM_INLINE>;
using LevVec = boost::container::small_vector<int, STRPM_INLINE>;

class STRPMSolver : public Solver
{
public:
    STRPMSolver(Oink& oink, Game& game);
    virtual ~STRPMSolver();

    virtual void run();

protected:
    /**
     * Parameters: U^k_{t, h}
     *      - k: Strahler-number
     *      - t: number of bits
     *      - h: height
     */
    int k, t, h;
    std::vector<BitVec> pm_b;
    std::vector<LevVec> pm_d;

    BitVec tmp_b;
    LevVec tmp_d;

    BitVec best_b;
    LevVec best_d;

    BitVec test_b;
    LevVec test_d;

    uintqueue Q;
    bitset dirty;
    bitset unstable;

    bool always_reset = false;

    int *cap; // caps!
    uint64_t *lift_counters;

    // Copy pm[idx] into tmp
    void to_tmp(int idx);
    // Copy tmp into pm[idx]
    void from_tmp(int idx);
    // Copy pm[idx] into best
    void to_best(int idx);
    // Copy best into pm[idx]
    void from_best(int idx);
    // Copy tmp into best
    void tmp_to_best();
    // Copy tmp into test
    void tmp_to_test();

    // Render pm[idx] to given ostream
    void stream_pm(std::ostream &out, int idx);
    // Render all occuring pms to ostream
    void print_full_pm(std::ostream &out);
    // Render tmp to given ostream
    void stream_tmp(std::ostream &out, int h);
    // Render best to given ostream
    void stream_best(std::ostream &out, int h);

    // Compare tmp to best
    int compare(int pindex, BitVec& other_b, LevVec& other_d);
    // Compare tmp to test
    int compare_test(int pindex);

    // Bump tmp given priority p
    void trunc_tmp(int pindex);
    int skipUntilNextLevel (LevVec& curr_d, int i);
    void prog_tmp(int pindex, int h);
    void prog_cap_tmp(int pindex);

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

class NaiveSTRPMSolver : public STRPMSolver
{
public:
    NaiveSTRPMSolver(Oink& oink, Game& game) : STRPMSolver(oink, game) { always_reset = true; }
    virtual ~NaiveSTRPMSolver() { }
};

}

#endif 
