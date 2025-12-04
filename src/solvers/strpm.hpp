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

namespace pg {

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
    std::vector<std::vector<bool>> pm_b;
    std::vector<std::vector<int>> pm_d;

    std::vector<bool> tmp_b;
    std::vector<int> tmp_d;

    std::vector<bool> best_b;
    std::vector<int> best_d;

    std::vector<bool> test_b;
    std::vector<int> test_d;

    uintqueue Q;
    bitset dirty;
    bitset unstable;

    bool bounded = false;

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
    // Render tmp to given ostream
    void stream_tmp(std::ostream &out, int h);
    // Render best to given ostream
    void stream_best(std::ostream &out, int h);

    // Compare tmp to best
    int compare(int pindex, std::vector<bool>& other_b, std::vector<int>& other_d);
    // Compare tmp to test
    int compare_test(int pindex);

    // Bump tmp given priority p
    void trunc_tmp(int pindex);
    int skipUntilNextLevel (std::vector<int>& curr_d, int i);
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

    void run(int nbits, int depth, int player);
};

class BoundedSTRPMSolver : public STRPMSolver
{
public:
    BoundedSTRPMSolver(Oink& oink, Game& game) : STRPMSolver(oink, game) { bounded = true; }
    virtual ~BoundedSTRPMSolver() { }
};

}

#endif 
