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
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <boost/functional/hash.hpp>
#include <stack>
#include <utility>

#include "strpm.hpp"

#define ODDFIRST 1
#define ALWAYS_RESET 0

namespace pg {

STRPMSolver::STRPMSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

STRPMSolver::~STRPMSolver()
{
}

struct RatioCompare {
    bool operator()(const std::pair<int,int>& lhs,
                    const std::pair<int,int>& rhs) const 
    {
        auto lhs_ratio = std::max(lhs.first, lhs.second)/std::min(lhs.first, lhs.second);
        auto rhs_ratio = std::max(rhs.first, rhs.second)/std::min(rhs.first, rhs.second);  
        
        return lhs_ratio > rhs_ratio;
    }
};

constexpr inline size_t binom(size_t n, size_t k) noexcept
{
    return
      (        k> n  )? 0 :          // out of range
      (k==0 || k==n  )? 1 :          // edge
      (k==1 || k==n-1)? n :          // first
      (     k+k < n  )?              // recursive:
      (binom(n-1,k-1) * n)/k :       //  path to k=1   is faster
      (binom(n-1,k) * n)/(n-k);      //  path to k=n-1 is faster
}

struct ApproxSizeCompare {
    int h;

    double approxSize(int k, int t) const 
    {
        if (k == 1) return 1;

        double approximation = (k-1) * (std::log(static_cast<double>(h-1)/(k-1)) + 1) - 
                               (std::log(2* M_PI * (k-1))/2) - ((k-1)*(k-1))/static_cast<double>(2*(h-1));
    
        return (k + t)*std::log(2.0) + std::log(binom(t + k - 2, k + 2)) + approximation;
    }

    bool operator()(const std::pair<int,int>& lhs,
                    const std::pair<int,int>& rhs) const 
    {
        auto lhs_size = approxSize(lhs.first, lhs.second);
        auto rhs_size = approxSize(rhs.first, rhs.second);

        return lhs_size > rhs_size;
    }
};

// Count distinct values in the non-decreasing range [first, last).
// tmp_d (the level array) is maintained non-decreasing, so a distinct-value
// count equals 1 + (number of adjacent changes) — no hash set / allocation
// needed. This replaces the per-call std::unordered_set constructions that
// previously dominated prog_tmp's allocation traffic.
static inline int
count_distinct_nondecreasing(LevVec::const_iterator first,
                             LevVec::const_iterator last)
{
    assert(std::is_sorted(first, last));
    if (first == last) return 0;
    int distinct = 1;
    for (auto it = first + 1; it != last; ++it)
        if (*it != *(it - 1)) ++distinct;
    return distinct;
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
    int  i=tmp_d.size()-1;
    while (i>=0 and tmp_d[i] > pindex) i--;
    tmp_b.resize(i + 1);
    tmp_d.resize(i + 1);
    assert ((tmp_d.size() - std::unordered_set<int>(tmp_d.begin(), tmp_d.end()).size()) <= t);
}

/**
 * Helper: skip over bits until the level changes
 */
int 
STRPMSolver::skipUntilNextLevel (LevVec& curr_d, int i)
{
    while ((i == curr_d.size() - 1) || (i >= 0 && curr_d[i] == curr_d[i+1])) 
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
        logger << "Calculating successor of: ";
        for (auto &&bit : tmp_b)
        {
            logger << bit;
        }
        logger << " ";
        for (auto &&loc : tmp_d)
        {
            logger << loc;
        }
        logger << " with k = " << k << ", t = " << t << std::endl;
    logger << "Start i in " << i << std::endl;
    logger << "Skipping bits below p " << pindex << std::endl;
    }
#endif
    assert (*std::max_element(tmp_d.begin(), tmp_d.end()) < h-1);
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
    int nes = count_distinct_nondecreasing(tmp_d.begin(), tmp_d.begin() + i + 1);
    assert (nes >= 0);
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
                int new_index = std::min(std::cmp_equal(i+1, tmp_d.size()) ? tmp_d[i] : tmp_d[i+1] - 1, pindex);
                assert(new_index < h-1);
                i ++;

                // Either we add the one at the end of the existing string, which means we add one NLB, or we add a NES!
                bool isNewString = false;
                if (i != 0 and (i == tmp_d.size() or tmp_d[i-1] == new_index)) nlb++;
                else {
                    nes++;
                    isNewString = true;
#ifndef NDEBUG
                    if (trace >= 2) logger <<  "Appending to an empty string...\n";
#endif
                }
                
                int total_nes = std::min((nes + h - 1 - new_index), k - 1);
                if (std::cmp_greater(total_nes + t, tmp_b.size()))
                {
                    size_t newBits = total_nes + t - tmp_d.size();
#ifndef NDEBUG
                    if (trace >= 2) logger <<  "Resizing to fit bits (" << newBits << " additional bits)\n";
#endif
                    
                    tmp_b.insert(tmp_b.end(), newBits, 0);
                    tmp_b[i] = 1;
                    for (size_t j = i; j < tmp_d.size(); j++) tmp_d[j] = new_index;
                    tmp_d.insert(tmp_d.end(), newBits, new_index);
                    i += t - nlb + 1;
                }
                else 
                {
#ifndef NDEBUG
                  if (trace >= 2) logger <<  "Adding bits\n";
#endif
                    tmp_b[i] = 1;
                    tmp_d[i] = new_index;
                    
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
                int strings_after_current = count_distinct_nondecreasing(tmp_d.begin() + i, tmp_d.end());
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
            tmp_d[i] = std::min(tmp_d[i], pindex+1) - 2;
            if (std::cmp_greater(t + k - 1, tmp_b.size()))
            {
                size_t newBits = t + k - 1 - tmp_d.size();
#ifndef NDEBUG
                if (trace >= 2) logger << "Vector too small, appending " << newBits << "bits\n";
#endif
                tmp_b.insert(tmp_b.end(), newBits, 0);
                for (size_t j = i; j < tmp_d.size(); j++) tmp_d[j] = tmp_d[i];
                tmp_d.insert(tmp_d.end(), newBits, tmp_d[i]);
            }
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
    
    // We can only fill as many NES as we have height... Too many empty strings "in the middle" = less tham k-1 strings total
    int available_nes = h - 1 - (tmp_d[(skipLevel ? i : i-1)] + 1);
    int no_of_needed_nes = std::min((k-1) - (nes+1), available_nes);
    
    if (no_of_needed_nes == 0) {
        // special case: Don't add anything, remove everything after the current position
        tmp_d.resize(i);
        tmp_b.resize(i);
    }
    else {
        int set_index = tmp_d[(skipLevel ? i : i-1)] + 1;
        // Fill up with just the next one as long as we still have bits
    #ifndef NDEBUG
        if (trace >= 2) logger << "Needed nes: " << no_of_needed_nes << ", available levels: " << available_nes << std::endl;
    #endif
        while (std::cmp_greater_equal(tmp_b.size(), i + no_of_needed_nes))
        {
            assert (i >= 0 and std::cmp_less(i, tmp_d.size()));
            assert(set_index < h-1);
#ifndef NDEBUG
            if (trace >= 2) logger << "Set bit " << i << " to " << set_index << std::endl;
#endif
            tmp_d[i] = set_index;
            i ++;
        }
#ifndef NDEBUG
        if (trace >= 2) logger << "Filling singles\n";
#endif
        // Now assign the rest of the bits one level a piece
        while (std::cmp_less(i, tmp_b.size()))
        {
            set_index ++;
#ifndef NDEBUG
            if (trace >= 2) logger << "Set bit " << i << " to " << set_index << std::endl;
#endif
            assert (set_index < h-1);
            tmp_d[i] = set_index;
            i ++;
        }
    }
    // Assert that the number of NLB is at most t
    assert ((tmp_d.size() - std::unordered_set<int>(tmp_d.begin(), tmp_d.end()).size()) <= t);
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
        for (int i=0; i<h-1; i++) {
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
        for (int i=0; i<h-1; i++) {
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
        for (int i=0; i<h-1; i++) {
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
STRPMSolver::compare(int pindex, BitVec& other_b, LevVec& other_d)
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
    const int pindex = pl == 0 ? (h-1)-(pr+1)/2-1 : (h-1)-pr/2-1;

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
                logger << "to target " << label_vertex(target) << "(" << target << ")" << ":";
                stream_tmp(logger, h);
                logger << " =>";
            }
#endif
        if (pl == (pr&1)) prog_tmp(pindex, h);
        //else trunc_tmp(pindex);
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
        //else trunc_tmp(pindex);
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

struct Node
{
    int k;
    int t;
    int h;
    bool isU;
};


// Keep track of already computed sizes, this is cached beyond single calls of
// tree_size below
std::unordered_map<std::tuple<int, int, int>, unsigned,
                   boost::hash<std::tuple<int, int, int>>> treeU;
std::unordered_map<std::tuple<int, int, int>, unsigned,
                   boost::hash<std::tuple<int, int, int>>> treeV;

unsigned tree_size(int k, int t, int h) 
{
    std::stack<Node> stack;

    stack.push ({k, t, h, true});
    while (!stack.empty())
    {
        Node& tos = stack.top();
        if (tos.isU and tos.h == 1 and tos.k == 1)
        {
            treeU[std::make_tuple(tos.k, tos.t, tos.h)] = 1;
            stack.pop ();
        }
        else if (tos.isU and tos.h > 1 and tos.k == 1)
        {
            auto son = treeU.find(std::make_tuple(tos.k, tos.t, tos.h - 1)); 
            if (son != treeU.end())
            {
                treeU[std::make_tuple(tos.k, tos.t, tos.h)] = son->second;
                stack.pop ();
            }
            else stack.push ({tos.k, tos.t, tos.h - 1, true});
        }
        else if (tos.h >= tos.k and tos.k >= 2 and tos.t == 0)
        {
            auto son = treeU.find(std::make_tuple(tos.k - 1, tos.t, tos.h - 1));
            if (son != treeU.end())
            {
                if (tos.isU) treeU[std::make_tuple(tos.k, tos.t, tos.h)] = son->second;
                else treeV[std::make_tuple(tos.k, tos.t, tos.h)] = son->second;
                stack.pop ();
            }
            else stack.push ({tos.k - 1, tos.t, tos.h - 1, true});
        }
        else if (!tos.isU and tos.h >= tos.k and tos.k >= 2 and tos.t >= 1)
        {
            auto son1 = treeV.find(std::make_tuple(tos.k, tos.t - 1, tos. h));
            auto son2 = treeU.find(std::make_tuple(tos.k - 1, tos.t, tos.h - 1));
            if (son1 != treeV.end() and son2 != treeU.end())
            {
                treeV[std::make_tuple(tos.k, tos.t, tos.h)] = son1->second * 2 + son2->second;
                stack.pop ();
            }
            else
            {
                stack.push ({tos.k - 1, tos.t, tos.h - 1, true});
                stack.push ({tos.k, tos.t - 1, tos.h, false});
            }
        }
        else if (tos.isU and tos.h == tos.k and tos.k >= 2)
        {
            auto son = treeV.find(std::make_tuple(tos.k, tos.t, tos.h));
            if (son != treeV.end()) 
            {
                treeU[std::make_tuple(tos.k, tos.t, tos.h)] = son->second;
                stack.pop ();
            }
            else stack.push ({tos.k, tos.t, tos.h, false});
        }
        else if (tos.isU and tos.h > tos.k and tos.k >= 2)
        {
            auto son1 = treeV.find(std::make_tuple(tos.k, tos.t, tos.h));
            auto son2 = treeU.find(std::make_tuple(tos.k, tos.t, tos.h - 1));
            if (son1 != treeV.end() and son2 != treeU.end())
            {
                treeU[std::make_tuple(tos.k, tos.t, tos.h)] = son1->second * 2 + son2->second;
                stack.pop ();
            }
            else
            {
                stack.push ({tos.k, tos.t, tos.h - 1, true});
                stack.push ({tos.k, tos.t, tos.h, false});
            }
        }
        else assert(false); // We should never get here
    }
    
    return treeU[std::make_tuple(k, t, h)];
}

struct SizeCompare
{
    int h;

    bool operator()(const std::pair<int,int>& lhs,
                    const std::pair<int,int>& rhs) const 
    {
        auto lhs_size = tree_size(lhs.first, lhs.second, h);
        auto rhs_size = tree_size(rhs.first, rhs.second, h);
        
        return lhs_size > rhs_size;
    }
};


void
STRPMSolver::run(int t_val, int k_val, int depth, int player)
{
    // Marcin's word: think of h as the number of priorities of the
    // opponent... PLUS ONE!
    t = t_val;
    h = depth + 1;  // FIXME: This is Guillermo's hack, the +1
    k = k_val;  // Maybe possible: std::min(t + 2, h);

#ifndef NDEBUG
    logger << "Strahler-tree parameters for player " << player << ": k = " << k << ", t = " << t << ", h = " << h << std::endl;
#endif

    // initialize progress measures - Every node is set to the smallest leaf in the tree
    pm_b = std::vector<BitVec> (nodecount(), BitVec(k-1+t, 0));
    LevVec initial_d (k-1+t, 0);
    for (size_t i = t + 1; i < initial_d.size(); i++)
    {
        initial_d[i] = initial_d[i-1] + 1;
    }
    pm_d = std::vector<LevVec> (nodecount(), initial_d);

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
    int t_max = floor_log2(nodecount());
    int h0 = (max_prio/2)+1;
    int h1 = (max_prio+1)/2;

    int h_max = std::max(h0, h1);
    int k_max = std::min(t_max + 2, h_max);

    // create datastructures
    Q.resize(nodecount());
    dirty.resize(nodecount());
    unstable.resize(nodecount());

    // Create a priority queue for (k, t) pairs and push init with (1, 1)
    std::priority_queue<
        std::pair<int,int>,
        std::vector<std::pair<int,int>>,
        RatioCompare
        //ApproxSizeCompare
    > pq { };
    pq.push({1, 1});
    /*
    To use SizeCompare:
    std::priority_queue<
        std::pair<int,int>,
        std::vector<std::pair<int,int>>,
        SizeCompare
    > pq { SizeCompare { h_max } };
    */

    // Keep track of already tried combinations
    std::unordered_set<std::pair<int, int>, boost::hash<std::pair<int, int>>> already_tried;

#ifndef NDEBUG
    logger << "Max t: " << t_max << ", max k: " << k_max << std::endl;
#endif

#if ALWAYS_RESET
    bitset initial_disabled { disabled };
    bitset initial_solved { game.getSolved() };
#endif
    
    while (!pq.empty()) {
        // Step 1: Get values
        auto [k_val, t_val] = pq.top();
        pq.pop();

        // Step 2: Reset the game - we want to know whether this combination can solve the game on its own
        lift_count = 0, lift_attempt = 0;
        uint64_t c;

#if ALWAYS_RESET
        game.reset_to_initial(initial_solved);
        reset_to_initial(initial_disabled);
#endif

#ifndef NDEBUG
        logger << "Currently unsolved: " << game.count_unsolved() << std::endl;
#endif

        // Step 3: Actually do the solving
        if (ODDFIRST) {
            // run odd counters
            run(t_val, k_val, h1, 1);
            c = game.count_unsolved();
#ifndef NDEBUG
            logger << "after odd, " << std::setw(9) << lift_count << " lifts, " << std::setw(9) << lift_attempt << " lift attempts, " << c << " unsolved left." << std::endl;
#endif
            // if now solved, no need to run odd counters
            if (c != 0)
            {
                // run even counters
                run(t_val, k_val, h0, 0);
                c = game.count_unsolved();
#ifndef NDEBUG
                logger << "after even, " << std::setw(9) << lift_count << " lifts, " << std::setw(9) << lift_attempt << " lift attempts, " << c << " unsolved left." << std::endl;
#endif
            }
            
        } else {
            // run even counters
            run(t_val, k_val, h0, 0);
            c = game.count_unsolved();
#ifndef NDEBUG
            logger << "after even, " << std::setw(9) << lift_count << " lifts, " << std::setw(9) << lift_attempt << " lift attempts, " << c << " unsolved left." << std::endl;
#endif
            // if now solved, no need to run odd counters
            if (c != 0)
            {
                // run odd counters
                run(t_val, k_val, h1, 1);
                c = game.count_unsolved();
#ifndef NDEBUG
                logger << "after odd, " << std::setw(9) << lift_count << " lifts, " << std::setw(9) << lift_attempt << " lift attempts, " << c << " unsolved left." << std::endl;
#endif
            }
        }

        // Step 4: Check whether we solved the game
        if (c == 0)
        {
            // We can stop, everything is solved!
            logger << "Solved with k = " << k_val << ", t = " << t_val << std::endl;
            break;
        }
        else if (k_val < k_max or t_val < t_max)
        {
            std::pair<int, int> candidate {k_val + 1, t_val};
            if (k_val + 1 <= k_max and already_tried.find(candidate) == already_tried.end()) 
            {
                pq.push(candidate);
                already_tried.insert(candidate);
            }

            candidate = { k_val, t_val + 1 };
            if (t_val + 1 <= t_max and already_tried.find(candidate) == already_tried.end()) 
            {
                pq.push(candidate);
                already_tried.insert(candidate);
            }
        }
    }
    
}

}
