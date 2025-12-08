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

#include <iomanip>
#include <unordered_set>

#include "strpm.hpp"

#define ODDFIRST 1

namespace pg {

STRPMSolver::STRPMSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

STRPMSolver::~STRPMSolver()
{
}

void
STRPMSolver::to_tmp(int idx)
{
    tmp_b = pm_b[idx];
    tmp_d = pm_d[idx];
}

void
STRPMSolver::from_tmp(int idx)
{
    pm_b[idx] = tmp_b;
    pm_d[idx] = tmp_d;
}

void
STRPMSolver::to_best(int idx)
{
    best_b = pm_b[idx];
    best_d = pm_d[idx];
}

void
STRPMSolver::from_best(int idx)
{
    pm_b[idx] = best_b;
    pm_d[idx] = best_d;
}

void
STRPMSolver::tmp_to_best()
{
    best_b = tmp_b;
    best_d = tmp_d;
}

void
STRPMSolver::tmp_to_test()
{
    test_b = tmp_b;
    test_d = tmp_d;
}

/**
 * Set tmp := min { m | m ==_p tmp }
 */
void
STRPMSolver::trunc_tmp(int pindex)
{
    if (tmp_d[0] == -1) return; // already Top
    // compute the lowest pindex >= p
    // [pindex],.,...,.. => [pindex],000
    // if pindex is the bottom, then this simply "buries" the remainder
    for (int i=tmp_d.size()-1; i>=0 and tmp_d[i] > pindex; i--) {
        tmp_b[i] = 0;
        tmp_d[i] = pindex+1;
    }
}

/**
 * Helper: skip over bits until the level changes
 */
int 
STRPMSolver::skipUntilNextLevel (std::vector<int>& curr_d, int i) 
{
    while ((i >= 0 && curr_d[i] == curr_d[i+1]) || i == curr_d.size() - 1) 
    {
        tmp_b[i] = 0;
        i --;
    }
    return i;
}

/**
 * Set tmp := min { m | m >_p tmp }
 */
void
STRPMSolver::prog_tmp(int pindex, int h)
{
    // Simple case 1: Top >_p Top
    if (tmp_d[0] == -1) return; // already Top

    bool skipLevel = false;
    int i = tmp_d.size() - 1;
#ifndef NDEBUG
    if (trace >= 2) 
    {
    logger << "Start i in " << i << std::endl;
    logger << "Skipping bits below p\n";
    }
#endif
    // skip bits "below p"
    while (i >= 0 && tmp_d[i] > pindex) 
    {
        tmp_b[i] = 0;
        i --;
    }
#ifndef NDEBUG
    if (trace >= 2) 
    {
    logger << "After skipping i is " << i << std::endl;

    logger << "Calculating NES ";
    }
#endif
    // Calculate number of Non-Empty Strings (NES): count unique values in tmp_d up until there
    size_t nes = std::unordered_set<int>( tmp_d.begin(), tmp_d.begin() + i + 1 ).size();
#ifndef NDEBUG
    if (trace >= 2) logger << nes << std::endl;
#endif
    // Subtract the string that we currently look at (A, B only refers to strings "above")
    if (nes - 1 == k - 1) {
        // A: No next sibling on this layer
#ifndef NDEBUG
        if (trace >= 2) logger << "Skipping a level\n";
#endif
        i = skipUntilNextLevel(tmp_d, i);
    }

    while (i >= 0) {
        // check if there was a level change
        if (i == tmp_d.size()-1 or tmp_d[i] != tmp_d[i+1]) 
        {
#ifndef NDEBUG
            if (trace >= 2) logger << "Handling level change\n";
#endif
            // Calculate the Non-Leading Bits (NLB): take the complete length and subtract the number of NES (every NES has one leading bit)
            int nlb = (i + 1) - nes;
            nes --;
#ifndef NDEBUG
            if (trace >= 2) logger << "NLB " << nlb << " i " << i << " NES " << nes << std::endl;
#endif
            if (nlb < t) 
            {
#ifndef NDEBUG
                if (trace >= 2) logger <<  "Smaller than t\n";
#endif
                int new_index = i == tmp_d.size() - 1 ? tmp_d[i] : tmp_d[i+1] - 1;
                i ++;
                if ((i + t - nlb + 1)  > tmp_b.size())
                {
                    tmp_b.insert(tmp_b.end(), t - nlb, 0);
                    tmp_b[i] = 1;
                    tmp_d.insert(tmp_d.end(), t - nlb, new_index);
                    i += t - nlb;
                }
                else 
                {
                    tmp_b[i] = 1;
                    tmp_d[i] = new_index;
                    if (i != 0 and tmp_d[i-1] == tmp_d[i]) nlb++;
                    else nes++;
                    int j = 1;
                    while (nlb + j <= t) 
                    {
                        tmp_b [i + j] = 0;
                        tmp_d [i + j] = new_index;
                        j ++;
                    }
                    // remember the last changed position
                    i += j;
                }
                break;
            }
            // Case B: Check if the current string is only leading bit (so all NLB are in strings 0 to r-1)
            else if (i == 0 or tmp_d[i-1] != tmp_d[i])
            {
                tmp_b[i] = 0;
                i --;
                continue;
            }
        }
        // For all following cases we know NLB == t        
        // Have to always check 0s - even when e.g. the first case applies
        if (tmp_b[i] == 0) 
        {
            if (i == 0 || tmp_d[i - 1] != tmp_d[i])
            {
#ifndef NDEBUG
                if (trace >= 2) logger << "Found a 0 in the beginning\n";
#endif
                // The 0 is either the first bit in total, or it is the first bit of that level
                int strings_after_current = std::unordered_set<int> (tmp_d.begin() + i, tmp_d.end()).size();
                if (strings_after_current == (h-1) - tmp_d[i])
                {
                    // All bitstrings after the current level are non-empty, we simply move on
                    // C: No sibling on this layer
                    i --;
                    continue;
                }
                else
                {
                    // F: We have to remove one NES: The string is 01^j, so the new string would be empty
                    nes--;
                    skipLevel = true;
                }
                // We can use this level: Break out of the loop and start changing at the current i
                break;
            }
            else {
                assert (tmp_d[i - 1] == tmp_d[i]);
#ifndef NDEBUG
                if (trace >= 2) logger << "A zero in the middle!\n";
#endif
                break;
            }
        }
        else 
        {
#ifndef NDEBUG
            if (trace >= 2) logger << "start setting things to 0\n";
#endif
            // We can already start setting everything there to 0 and change where it belongs later
            tmp_b[i] = 0;
            i --;
        }
    }

    // A base case, we have reached the root: we are top
    // Special case of D
    if (i == -1) 
    {
        if (tmp_d[0] == 0)
        {
#ifndef NDEBUG
            if (trace >= 2) logger << "We are at top\n";
#endif
            tmp_d[0] = -1;
            return;
        }
        else
        {
            // Special case: There are still empty strings "left" that can be filled instead
            i = 0;
            tmp_b[i] = 1;
            // TODO: maybe find more elegant way to avoid setting the index too high as opposed to starting out one lower...
            tmp_d[i] = tmp_d[i] - 2;
            skipLevel = true;
            nes--;
        }
    }

    // Change where the bits belong
    #ifndef NDEBUG
    if (trace >= 2)
    {
    logger << "Adjusting bit level\n";
    logger << "starting at " << i << std::endl; 
    }
    #endif
    
    int no_of_needed_nes = (k-1) - (nes+1);
    if (no_of_needed_nes == 0) {
        // special case: Don't add anything, remove everything after the current position
        tmp_d.resize(i);
        tmp_b.resize(i);
    }
    else {
        int set_index = tmp_d[(skipLevel ? i : i-1)] + 1;
        // Fill up with just the next one as long as we still have bits
    #ifndef NDEBUG
        if (trace >= 2) logger << "Needed nes: " << no_of_needed_nes << std::endl;
    #endif
        while (tmp_b.size() - i >= no_of_needed_nes) 
        {
            assert (i >= 0 && i < tmp_d.size());
            tmp_d[i] = set_index;
            i ++;
        }
#ifndef NDEBUG
        if (trace >= 2) logger << "Filling singles\n";
#endif
        // Now assign the rest of the bits one level a piece
        while (i < tmp_b.size())
        {
            set_index ++;
            tmp_d[i] = set_index;
            i ++;
        }
    }
}

/**
 * Write pm to ostream.
 */
void
STRPMSolver::stream_pm(std::ostream &out, int idx)
{
    if (pm_d[idx][0] == -1) {
        out << " \033[1;33mTop\033[m";
    } else {
        out << " { ";
        int j=0;
        for (int i=0; i<h; i++) {
            if (i>0) out << ",";
            int c=0;
            while (j<pm_d[idx].size() and pm_d[idx][j] == i) {
                c++;
                out << pm_b[idx][j];
                j++;
            }
            if (c == 0) out << "ε";
        }
        out << " }";
    }
}

/**
 * Write tmp to ostream.
 */
void
STRPMSolver::stream_tmp(std::ostream &out, int h)
{
    if (tmp_d[0] == -1) {
        out << " \033[1;33mTop\033[m";
    } else {
        out << " { ";
        int j=0;
        for (int i=0; i<h; i++) {
            if (i>0) out << ",";
            int c=0;
            while (j<tmp_b.size() and tmp_d[j] == i) {
                c++;
                out << tmp_b[j];
                j++;
            }
            if (c == 0) out << "ε";
        }
        out << " }";

        /* Not sure we want to keep this
        out << " {";

        // compute value
        int i=0;
        for (int d=0; d<h; d++) {
            int val = 0;

            for (; i<tmp_b.size(); i++) {
                if (tmp_d[i] != d) {
                    // e found
                    val |= ((1 << (tmp_b.size()-i)) - 1);
                    break;
                }

                if (tmp_b[i]) val |= (1 << (tmp_b.size()-i));
            }

            logger << " " << val;
        }

        out << " }";*/
    }
}

/**
 * Write best to ostream.
 */
void
STRPMSolver::stream_best(std::ostream &out, int h)
{
    if (best_d[0] == -1) {
        out << " \033[1;33mTop\033[m";
    } else {
        out << " { ";
        int j=0;
        for (int i=0; i<h; i++) {
            if (i>0) out << ",";
            int c=0;
            while (j<best_b.size() and best_d[j] == i) {
                c++;
                out << best_b[j];
                j++;
            }
            if (c == 0) out << "ε";
        }
        out << " }";
    }
}

/**
 * Compare tmp and other
 * res := -1 :: tmp < other
 * res := 0  :: tmp = other
 * res := 1  :: tmp > other
 */
int
STRPMSolver::compare(int pindex, std::vector<bool>& other_b, std::vector<int>& other_d)
{
    // cases involving Top
    if (tmp_d[0] == -1 and other_d[0] == -1) return 0;
    if (tmp_d[0] == -1) return 1;
    if (other_d[0] == -1) return -1;

    for (int i=0; i<std::max(tmp_d.size(), other_d.size()); i++) {
        if (i >= other_b.size())
        {
            return tmp_b[i] == 0 ? -1: 1;
        }
        else if (i >= tmp_b.size())
        {
            return other_b[i] == 0 ? 1 : -1;
        }
        else if (tmp_d[i] > pindex and other_d[i] > pindex) {
            // equal until pindex, return 0
            return 0;
        } else if (tmp_d[i] < other_d[i]) {
            // equal until other has [eps]
            return tmp_b[i] == 0 ? -1 : 1;
        } else if (tmp_d[i] > other_d[i]) {
            // equal until tmp has [eps]
            return other_b[i] == 0 ? 1: -1;
        } else if (tmp_b[i] < other_b[i]) {
            // equal until tmp<other
            return -1;
        } else if (tmp_b[i] > other_b[i]) {
            // equal until tmp>other
            return 1;
        }
    }
    return 0;
}

bool
STRPMSolver::lift(int v, int target, int &str, int pl)
{
    // check if already Top
    if (pm_d[v][0] == -1) return false; // already Top

    const int pr = priority(v);
    const int pindex = pl == 0 ? h-(pr+1)/2-1 : h-pr/2-1;

#ifndef NDEBUG
    if (trace >= 2) {
        logger << "\033[37;1mupdating vertex " << label_vertex(v) << " (" << pr << " => " << pindex << ")" << (owner(v)?" (odd)":" (even)") << "\033[m with current measure";
        stream_pm(logger, v);
        logger << std::endl;
    }
#endif

    // if even owns and target is set, just check if specific target is better
    if (owner(v) == pl and target != -1) {
        to_tmp(target);
#ifndef NDEBUG
            if (trace >= 2) {
                logger << "to target " << label_vertex(target) << ":";
                stream_tmp(logger, h);
                logger << " =>";
            }
#endif
        if (pl == (pr&1)) prog_tmp(pindex, h);
        else trunc_tmp(pindex);
#ifndef NDEBUG
            if (trace >= 2) {
                stream_tmp(logger, h);
                logger << std::endl;
            }
#endif
        to_best(v);
        if (compare(pindex, best_b, best_d) > 0) {
            from_tmp(v);
#ifndef NDEBUG
            if (trace >= 1) {
                logger << "\033[32;1mnew measure\033[m of \033[36;1m" << label_vertex(v) << "\033[m:";
                stream_tmp(logger, h);
                logger << " (to " << label_vertex(target) << ")\n";
            }
#endif
            return true;
        } else {
            return false;
        }
    }

    // compute best measure
    bool first = true;
    for (auto curedge = outs(v); *curedge != -1; curedge++) {
        int to = *curedge;
        if (disabled[to]) continue;
        to_tmp(to);
#ifndef NDEBUG
        if (trace >= 2) {
            logger << "to successor " << label_vertex(to) << " from";
            stream_tmp(logger, h);
            logger << " =>";
        }
#endif
        if (pl == (pr&1)) prog_tmp(pindex, h);
        else trunc_tmp(pindex);
#ifndef NDEBUG
        if (trace >= 2) {
            stream_tmp(logger, h);
            logger << std::endl;
        }
#endif
        if (first) {
            tmp_to_best();
            str = to;
        } else if (owner(v) == pl) {
            // we want the max!
            if (compare(pindex, best_b, best_d) > 0) {
                tmp_to_best();
                str = to;
            }
        } else {
            // we want the min!
            if (compare(pindex, best_b, best_d) < 0) {
                tmp_to_best();
                str = to;
            }
        }
        first = false;
    }

    // set best to pm if higher
    to_tmp(v);
    if (compare(pindex, best_b, best_d) < 0) {
#ifndef NDEBUG
        if (trace >= 1) {
            logger << "\033[1;32mnew measure\033[m of \033[36;1m" << label_vertex(v) << "\033[m:";
            stream_best(logger, h);
            logger << " (to " << label_vertex(str) << ")\n";
        }
#endif
        from_best(v);
        return true;
    } else {
        return false;
    }
}

static int
ceil_log2(unsigned long long x)
{
    static const unsigned long long t[6] = {
        0xFFFFFFFF00000000ull,
        0x00000000FFFF0000ull,
        0x000000000000FF00ull,
        0x00000000000000F0ull,
        0x000000000000000Cull,
        0x0000000000000002ull
    };

    int y = (((x & (x - 1)) == 0) ? 0 : 1);
    int j = 32;
    int i;

    for (i = 0; i < 6; i++) {
        int k = (((x & t[i]) == 0) ? 0 : j);
        y += k;
        x >>= k;
        j >>= 1;
    }

    return y;
}

static int
floor_log2 (unsigned long long x)
{
    static const unsigned long long t[6] = {
        0xFFFFFFFF00000000ull,
        0x00000000FFFF0000ull,
        0x000000000000FF00ull,
        0x00000000000000F0ull,
        0x000000000000000Cull,
        0x0000000000000002ull
    };

    int y = 0;             // no +1 for non-powers of two[2][1]
    int j = 32;

    for (int i = 0; i < 6; i++) {
        int k = (((x & t[i]) == 0) ? 0 : j);
        y += k;
        x >>= k;
        j >>= 1;
    }

    return y;
}

void
STRPMSolver::run(int n_bits, int depth, int player)
{
    // Marcin's word: think of h as the number of priorities of the
    // opponent... PLUS ONE!
    t = n_bits;
    h = depth + 1;  // FIXME: This is Guillermo's hack, the +1
    k = std::min(n_bits + 2, h);  // Maybe possible: std::min(t + 2, h);

    logger << "Strahler-tree parameters for player " << player << ": k = " << k << ", t = " << t << ", h = " << h << std::endl;

    // initialize progress measures - Every node is set to the smallest leaf in the tree
    pm_b = std::vector<std::vector<bool>> (nodecount(), std::vector<bool>(k-1+t, 0));
    std::vector<int> initial_d (k-1+t, 0);
    for (size_t i = t + 2; i < initial_d.size(); i++)
    {
        initial_d[i] = initial_d[i-1] + 1;
    }
    pm_d = std::vector<std::vector<int>> (nodecount(), initial_d);

#ifndef NDEBUG
    if (trace >= 1)
    {
        logger << "Initial PM: " << std::endl;
        stream_pm(logger, 0);
        logger << std::endl;
    }
#endif

    for (int n=nodecount()-1; n>=0; n--) {
        if (disabled[n]) continue;
        lift_attempt++;
        int s;
        if (lift(n, -1, s, player)) {
            lift_count++;
            // lift_counters[n]++;
            for (auto curedge = ins(n); *curedge != -1; curedge++) {
                int from = *curedge;
                if (disabled[from]) continue;
                lift_attempt++;
                int s;
                if (lift(from, n, s, player)) {
                    lift_count++;
                    // lift_counters[from]++;
                    todo_push(from);
                }
            }
        }
    }

    while (!Q.empty()) {
        int n = todo_pop();
        for (auto curedge = ins(n); *curedge != -1; curedge++) {
            int from = *curedge;
            if (disabled[from]) continue;
            lift_attempt++;
            int s;
            if (lift(from, n, s, player)) {
                lift_count++;
                // lift_counters[from]++;
                todo_push(from);
            }
        }
    }

    /**
     * Derive strategies.
     */

    for (int v=0; v<nodecount(); v++) {
        if (disabled[v]) continue;
        if (pm_d[v][0] != -1) {
            if (owner(v) != player) {
                // TODO: don't rely on the strategy array in the Game class
                if (lift(v, -1, game.getStrategy()[v], player)) logger << "error: " << v << " is not progressive!" << std::endl;
            }
        }
    }

    if (trace) {
        for (int v=0; v<nodecount(); v++) {
            if (disabled[v]) continue;

            logger << "\033[1m" << label_vertex(v) << (owner(v)?" (odd)":" (even)") << "\033[m:";
            stream_pm(logger, v);

            if (pm_d[v][0] != -1) {
                if (owner(v) != player) {
                    logger << " => " << label_vertex(game.getStrategy(v));
                }
            }

            logger << std::endl;
        }
    }

    /**
     * Mark solved.
     */

    for (int v=0; v<nodecount(); v++) {
        if (disabled[v]) continue;
        if (pm_d[v][0] != -1) Solver::solve(v, 1-player, game.getStrategy(v));
    }

    Solver::flush();
}

void
STRPMSolver::run()
{
    int max_prio = priority(nodecount()-1);

    // compute ml (max l) and the h for even/odd
    int ml = floor_log2(nodecount());
    int h0 = (max_prio/2)+1;
    int h1 = (max_prio+1)/2;

    // create datastructures
    Q.resize(nodecount());
    dirty.resize(nodecount());
    unstable.resize(nodecount());

    // if running bounded STRPM, start with 1-bounded adaptive counters
    int i = bounded ? 1 : ml;

    for (; i<=ml; i++) {
        int _l = lift_count, _a = lift_attempt;
        uint64_t _c = game.count_unsolved(), c;

        if (ODDFIRST) {
            // run odd counters
            run(i, h1, 1);
            c = game.count_unsolved();
            logger << "after odd  with i=" << i << ", " << std::setw(9) << lift_count-_l << " lifts, " << std::setw(9) << lift_attempt-_a << " lift attempts, " << c << " unsolved left." << std::endl;

            // if now solved, no need to run odd counters
            if (c == 0) break;

            // run even counters
            run(i, h0, 0);
            c = game.count_unsolved();
            logger << "after even with i=" << i << ", " << std::setw(9) << lift_count-_l << " lifts, " << std::setw(9) << lift_attempt-_a << " lift attempts, " << c << " unsolved left." << std::endl;
        } else {
            // run even counters
            run(i, h0, 0);
            c = game.count_unsolved();
            logger << "after even with i=" << i << ", " << std::setw(9) << lift_count-_l << " lifts, " << std::setw(9) << lift_attempt-_a << " lift attempts, " << c << " unsolved left." << std::endl;

            // if now solved, no need to run odd counters
            if (c == 0) break;

            // run odd counters
            run(i, h1, 1);
            c = game.count_unsolved();
            logger << "after odd  with i=" << i << ", " << std::setw(9) << lift_count-_l << " lifts, " << std::setw(9) << lift_attempt-_a << " lift attempts, " << c << " unsolved left." << std::endl;
        }

        if (i != ml) {
            // if i == ml then we are guaranteed to be done
            // otherwise check if done
            if (c == 0) break;
            if (_c != c) i--; // do not increase i if we solved vertices with current i
        } else {
            break; // do not count higher pls
        }
    }

    logger << "solved with " << lift_count << " lifts, " << lift_attempt << " lift attempts, max l " << i << "." << std::endl;
}

}
