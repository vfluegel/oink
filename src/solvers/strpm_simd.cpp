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

#include "strpm_simd.hpp"

#define ODDFIRST 1
#define ALWAYS_RESET 0

namespace pg {

STRPM_SIMDSolver::STRPM_SIMDSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

STRPM_SIMDSolver::~STRPM_SIMDSolver()
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

/**
 * Load up to `count` successor PMs (given by node indices ids[0..count-1])
 * into the interleaved batch registers. PM j occupies lanes j*(k-1)+i for
 * string position i. Slots j >= count are zeroed (sentinel: nlanes = 0).
 */
void
STRPM_SIMDSolver::load_batch(const int* ids, int count)
{
    const int stride = std::max(k - 1, 1);
    alignas(32) uint8_t b[32] = {}, m[32] = {};
    for (int j = 0; j < batch_B; j++) {
        if (j < count) {
            int node = ids[j];
            batch_ids[j]    = node;
            batch_nlanes[j] = pm_nlanes[node];
            for (int i = 0; i < stride; i++) {
                int lane = j * stride + i;
                b[lane]            = pm_bits  [node * 8 + i];
                m[lane]            = pm_masks [node * 8 + i];
                batch_levels[lane] = pm_levels[node * 8 + i];
            }
        } else {
            batch_ids[j]    = -1;
            batch_nlanes[j] = 0;
            for (int i = 0; i < stride; i++)
                batch_levels[j * stride + i] = 0;
            // b[] and m[] already zero-initialised for these lanes
        }
    }
    batch_bits .copy_from(b, stdx::element_aligned);
    batch_masks.copy_from(m, stdx::element_aligned);
}

/**
 * Copy PM j from the interleaved batch registers into the tmp_* working registers.
 */
void
STRPM_SIMDSolver::extract_batch_to_tmp(int j)
{
    const int stride = std::max(k - 1, 1);
    alignas(8) uint8_t b[8] = {}, m[8] = {};
    for (int i = 0; i < k - 1; i++) {
        int l  = j * stride + i;
        b[i]   = static_cast<uint8_t>(batch_bits [l]);
        m[i]   = static_cast<uint8_t>(batch_masks[l]);
        tmp_levels[i] = batch_levels[l];
    }
    // k == 1: the loop above runs 0 times (k-1=0) but levels[0] must still be
    // propagated.  batch_prog_tmp sets batch_levels[j]=-1 for Top, and
    // load_batch copies pm_levels[succ*8] into batch_levels[j], so when
    // batch_nlanes[j]>0 the level is valid (e.g. -1 for Top).
    if (k == 1 && batch_nlanes[j] > 0) tmp_levels[0] = batch_levels[j];
    tmp_bits .copy_from(b, stdx::element_aligned);
    tmp_masks.copy_from(m, stdx::element_aligned);
    tmp_nlanes = batch_nlanes[j];
    fill_inactive_tmp();
}

/**
 * Copy PM j from the interleaved batch registers into the best_* working registers.
 */
void
STRPM_SIMDSolver::extract_batch_to_best(int j)
{
    const int stride = std::max(k - 1, 1);
    alignas(8) uint8_t b[8] = {}, m[8] = {};
    for (int i = 0; i < k - 1; i++) {
        int l  = j * stride + i;
        b[i]   = static_cast<uint8_t>(batch_bits [l]);
        m[i]   = static_cast<uint8_t>(batch_masks[l]);
        best_levels[i] = batch_levels[l];
    }
    // k == 1: same issue as extract_batch_to_tmp — levels[0] must be copied.
    if (k == 1 && batch_nlanes[j] > 0) best_levels[0] = batch_levels[j];
    best_bits .copy_from(b, stdx::element_aligned);
    best_masks.copy_from(m, stdx::element_aligned);
    best_nlanes = batch_nlanes[j];
}

/**
 * Compute prog for cnt PMs simultaneously using the interleaved layout.
 *
 * cnt = number of valid PMs currently loaded (≤ batch_B).  All scalar loops
 * are bounded by cnt so that unused slots (j ≥ cnt) are never iterated.
 *
 * Predication/blending idiom: all bit/mask mutations are accumulated into
 * working byte arrays wb[] / wm[] during the scalar per-PM loop; a single
 * 32-wide SIMD copy_from at the end applies every update in one pass.
 */
void
STRPM_SIMDSolver::batch_prog_tmp(int pindex, int h, int cnt)
{
    const int stride = std::max(k - 1, 1);

    // k = 1: each PM is either Bottom (nlanes=0 → Top) or already Top.
    if (k == 1) {
        for (int j = 0; j < cnt; j++) {
            if (batch_nlanes[j] == 0) { batch_nlanes[j] = 1; batch_levels[j] = -1; }
        }
        return;
    }

    // ── 32-wide SIMD NLB per-element popcount ────────────────────────────────
    simd_uint8_32 per_elem = simd_popcount32(batch_masks);
    where(batch_masks > simd_uint8_32(0), per_elem) -= simd_uint8_32(1);

    // =====================================================================
    // k == 2 fast path: stride = 1, so lane j == PM j.
    //
    // Saves the NLB prefix-sum (nlb_b=0 always, nlb_c=per_elem) and the
    // match-finding loop (merged into the case loop below).  Uses wb/wm
    // byte arrays — not SIMD element access — for scalar writes, then a
    // single copy_from pair at the end to apply all changes in bulk.
    // =====================================================================
    if (k == 2) {
        // stp / afa (scalar, O(cnt), no inner loop for k=2)
        alignas(32) uint8_t stp_arr[32] = {}, afa_arr[32] = {};
        for (int j = 0; j < cnt; j++) {
            if (batch_nlanes[j] == 0 || batch_levels[j] == -1) continue;
            stp_arr[j] = (batch_levels[j] <= pindex) ? 1 : 0;
            // nlanes=1, i=0: afa iff (h-1-level) == 1 iff level == h-2
            afa_arr[j] = (batch_levels[j] == h - 2) ? 1 : 0;
        }
        simd_uint8_32 stp_v, afa_v;
        stp_v.copy_from(stp_arr, stdx::element_aligned);
        afa_v.copy_from(afa_arr, stdx::element_aligned);

        // has_succ (32-wide SIMD)
        // nlb_b[j]=0 always → (nlb_b==t) fires only if t==0 (invalid); omit.
        simd_uint8_32_mask has_bits_m = (batch_masks > simd_uint8_32(0));
        simd_uint8_32_mask smaller_p  = has_bits_m & (stp_v > simd_uint8_32(0));
        simd_uint8_32      pat01      = batch_masks & simd_uint8_32(0xFE);
        simd_uint8_32_mask no_succ    = has_bits_m
            & (per_elem == simd_uint8_32(static_cast<uint8_t>(t)))
            & (((batch_bits == pat01) & (afa_v > simd_uint8_32(0)))
               | (batch_bits == batch_masks));
        simd_uint8_32_mask has_succ = smaller_p & !no_succ;

        // Working byte arrays: plain scalar writes, no SIMD element access.
        alignas(32) uint8_t wb[32], wm[32];
        batch_bits .copy_to(wb, stdx::element_aligned);
        batch_masks.copy_to(wm, stdx::element_aligned);

        const uint8_t tail_mask = static_cast<uint8_t>((1u << (t + 1)) - 1);

        for (int j = 0; j < cnt; j++) {
            if (batch_nlanes[j] == 0 || batch_levels[j] == -1) continue;
            const uint8_t bv = wb[j], mv = wm[j];

            if (has_succ[j]) {
                if (per_elem[j] >= static_cast<uint8_t>(t)) {
                    // Case A: all t NLB slots used
                    int rb = std::countl_one(static_cast<uint8_t>(bv | ~mv)) + 1;
                    uint8_t keep = (rb < 8) ? static_cast<uint8_t>((1u << (8 - rb)) - 1) : 0u;
                    wb[j] = bv & keep;
                    wm[j] = mv & keep;
                    if (rb == 8) {
                        // Full reset: level++, tail fill
                        batch_levels[j] += 1;
                        wm[j] = tail_mask;
                    }
                    // rb < 8: partial reset, mi → 1 == k-1, no tail fill
                } else {
                    // Case B: extend string + tail fill
                    wb[j] = bv | static_cast<uint8_t>(1u << std::popcount(mv));
                    wm[j] = tail_mask;
                }
            } else {
                // mi == -1: carry or overflow
                int new_level = std::min(batch_levels[j], pindex + 1) - 1;
                if (new_level < 0) {
                    wb[j] = 0; wm[j] = 0;
                    batch_nlanes[j] = 1;
                    batch_levels[j] = -1;
                } else {
                    // Carry: leading "1" + tail fill
                    wb[j]           = 1;
                    wm[j]           = tail_mask;
                    batch_levels[j] = new_level;
                }
            }
        }

        batch_bits .copy_from(wb, stdx::element_aligned);
        batch_masks.copy_from(wm, stdx::element_aligned);
        return;
    }

    // ── Scalar NLB prefix sum (O(cnt × stride)) ───────────────────────────────
    alignas(32) uint8_t nlb_b_arr[32] = {}, nlb_c_arr[32] = {};
    for (int j = 0; j < cnt; j++) {
        int run = 0;
        for (int i = 0; i < k - 1; i++) {
            int l = j * stride + i;
            nlb_b_arr[l] = static_cast<uint8_t>(run);
            run          += static_cast<int>(per_elem[l]);
            nlb_c_arr[l] = static_cast<uint8_t>(run);
        }
    }

    // ── Scalar stp / afa predicates (O(cnt × stride)) ────────────────────────
    alignas(32) uint8_t stp_arr[32] = {}, afa_arr[32] = {};
    for (int j = 0; j < cnt; j++) {
        for (int i = 0; i < static_cast<int>(batch_nlanes[j]); i++) {
            int l      = j * stride + i;
            stp_arr[l] = (batch_levels[l] <= pindex) ? 1 : 0;
            afa_arr[l] = ((h - 1 - batch_levels[l]) == (static_cast<int>(batch_nlanes[j]) - i)) ? 1 : 0;
        }
    }
    simd_uint8_32 stp_v, afa_v, nlb_b_v, nlb_c_v;
    stp_v  .copy_from(stp_arr,   stdx::element_aligned);
    afa_v  .copy_from(afa_arr,   stdx::element_aligned);
    nlb_b_v.copy_from(nlb_b_arr, stdx::element_aligned);
    nlb_c_v.copy_from(nlb_c_arr, stdx::element_aligned);

    // ── 32-wide SIMD has_succ mask ────────────────────────────────────────────
    simd_uint8_32_mask has_bits_m = (batch_masks > simd_uint8_32(0));
    simd_uint8_32_mask smaller_p  = has_bits_m & (stp_v > simd_uint8_32(0));
    simd_uint8_32      pat01      = batch_masks & simd_uint8_32(0xFE);
    simd_uint8_32_mask no_succ    = has_bits_m & (
        (nlb_b_v == simd_uint8_32(static_cast<uint8_t>(t))) |
        ((nlb_c_v == simd_uint8_32(static_cast<uint8_t>(t))) & (
            ((batch_bits == pat01) & (afa_v > simd_uint8_32(0))) |
            (batch_bits == batch_masks)
        ))
    );
    simd_uint8_32_mask has_succ = smaller_p & !no_succ;

    // ── Working byte arrays: all bit/mask mutations go here ───────────────────
    // Initialised from the current register state; a single copy_from at the
    // end of the function writes every change back in one 32-wide SIMD store.
    alignas(32) uint8_t wb[32], wm[32];
    batch_bits .copy_to(wb, stdx::element_aligned);
    batch_masks.copy_to(wm, stdx::element_aligned);

    // ── Per-PM: match finding + Cases A / B + tail fill  (O(cnt × stride)) ───
    std::fill(batch_match.begin(), batch_match.begin() + cnt, -1);

    for (int j = 0; j < cnt; j++) {
        if (batch_nlanes[j] == 0) continue;
        if (batch_levels[j * stride] == -1) continue;  // already Top

        // Rightmost has_succ lane within this PM's stride window
        int mi = -1;
        for (int i = static_cast<int>(batch_nlanes[j]) - 1; i >= 0; i--)
            if (has_succ[j * stride + i]) { mi = i; break; }

        if (mi == -1) {
            // ── mi == -1: no lane can be incremented ─────────────────────────
            int new_level = std::min(batch_levels[j * stride], pindex + 1) - 1;
            if (new_level < 0) {
                // Overflow to Top: zero all lanes for this PM in wb / wm
                batch_nlanes[j]          = 1;
                batch_levels[j * stride] = -1;
                for (int i = 0; i < stride; i++)
                    wb[j * stride + i] = wm[j * stride + i] = 0;
                continue;
            }
            // Carry: leading "1" at new_level; tail fill handles mask / nlanes
            wb[j * stride]           = 1;
            batch_levels[j * stride] = new_level;
            mi = 0;
            // fall through to tail fill

        } else {
            int l = j * stride + mi;

            if (nlb_c_arr[l] >= static_cast<uint8_t>(t)) {
                // ── Case A: all t NLB slots used; carry / reset ───────────────
                uint8_t b_val = wb[l], m_val = wm[l];
                int reset_bits = std::countl_one(static_cast<uint8_t>(b_val | ~m_val)) + 1;
                uint8_t keep   = (reset_bits < 8)
                                 ? static_cast<uint8_t>((1u << (8 - reset_bits)) - 1)
                                 : 0u;
                wb[l] = b_val & keep;
                wm[l] = m_val & keep;

                if (reset_bits < 8) {
                    // Partial reset: advance match, update NLB count for tail fill
                    int new_pop    = std::popcount(wm[l]);
                    nlb_c_arr[l]   = static_cast<uint8_t>(
                                         nlb_b_arr[l] + (new_pop > 0 ? new_pop - 1 : 0));
                    mi++;
                    l = j * stride + mi;
                    if (mi < k - 1) {
                        batch_levels[l] = batch_levels[l - 1] + 1;
                        wb[l] = 0;
                    }
                } else {
                    // Full reset: lane stays at mi, move to next available level
                    nlb_c_arr[l] = nlb_b_arr[l];  // lane contributed nothing
                    int mi1  = mi + 1;
                    bool gap = (mi1 == static_cast<int>(batch_nlanes[j]) && batch_levels[l] < h - 2) ||
                               (mi1 <  static_cast<int>(batch_nlanes[j]) &&
                                batch_levels[l] + 1 < batch_levels[j * stride + mi1]);
                    batch_levels[l] = gap
                        ? ((mi1 == static_cast<int>(batch_nlanes[j]))
                               ? batch_levels[l] + 1
                               : batch_levels[j * stride + mi1] - 1)
                        : batch_levels[l] + 1;
                }
            } else {
                // ── Case B: NLB available; extend or start a new string ───────
                int mi1 = mi + 1;
                if (mi1 < static_cast<int>(batch_nlanes[j]) &&
                    batch_levels[l] + 1 < batch_levels[j * stride + mi1] &&
                    batch_levels[l] != pindex)
                {
                    // Gap: start a new string at the next level
                    mi++;
                    l               = j * stride + mi;
                    wb[l]           = 1;
                    batch_levels[l] = std::min(pindex, batch_levels[l] - 1);
                } else {
                    // Extend: set the next unused bit position
                    int first_new = std::popcount(wm[l]);
                    wb[l] |= static_cast<uint8_t>(1u << first_new);
                }
            }
        }

        // ── Tail fill from mi onward ──────────────────────────────────────────
        if (mi < k - 1) {
            int l           = j * stride + mi;
            int bits_before = (mi > 0) ? static_cast<int>(nlb_c_arr[j * stride + mi - 1]) : 0;
            wm[l]           = static_cast<uint8_t>((1u << (t - bits_before + 1)) - 1);

            // Extend nlanes to fill consecutive levels after mi
            batch_nlanes[j] = static_cast<uint8_t>(mi + 1);
            int stl         = batch_levels[l] + 1;
            while (static_cast<int>(batch_nlanes[j]) < k - 1 && stl <= h - 2) {
                batch_levels[j * stride + batch_nlanes[j]] = stl;
                batch_nlanes[j]++;
                stl++;
            }
            // mask=1, bits=0 for positions (mi, nlanes)
            for (int i = mi + 1; i < static_cast<int>(batch_nlanes[j]); i++)
                { wm[j * stride + i] = 1; wb[j * stride + i] = 0; }
            // Zero positions [nlanes, stride)
            for (int i = static_cast<int>(batch_nlanes[j]); i < k - 1; i++)
                { wb[j * stride + i] = 0; wm[j * stride + i] = 0; }
        }
    }

    // ── Blend: apply all accumulated updates in one 32-wide SIMD store ───────
    batch_bits .copy_from(wb, stdx::element_aligned);
    batch_masks.copy_from(wm, stdx::element_aligned);
}

/**
 * Set tmp := min { m | m >_p tmp }
 */
void
STRPM_SIMDSolver::prog_tmp(int pindex, int h)
{
    if (k == 1 and tmp_nlanes == 0)
    {
        // We can immediately handle k = 1, it's just one branch and top
        tmp_levels[0] = -1;
        tmp_nlanes = 1;
        tmp_bits = 0;
        tmp_masks = 0;
        fill_inactive_tmp();
        return;
    }

    // Simple case 1: Top >_p Top
    if (tmp_levels[0] == -1) return; // already Top

    // --- NLB (Non-Leading-Bit) counting across all 8 bitstring lanes ---
    // Each bitstring is stored as bits/mask pair. The mask has a 1 for each
    // position in the string, and bits holds the actual bit values.
    // NLB count for a single string = popcount(mask) - 1 (subtract the leading bit).
    // For empty strings (mask==0), NLB contribution is 0.
    simd_uint8_mask has_bits = (tmp_masks > 0);
    simd_uint8 per_elem = simd_popcount8(tmp_masks);
    stdx::where(has_bits, per_elem) = per_elem - simd_uint8(1);
    // Build inclusive prefix sum: nlb_counts[i] = total NLB in lanes 0..i
    simd_uint8 nlb_counts;
    nlb_counts[0] = per_elem[0];
    for (size_t n = 1; n < 8; n++) {
        nlb_counts[n] = nlb_counts[n-1] + per_elem[n];
    }
    // nlb_before[i] = total NLB in lanes 0..i-1 (exclusive prefix sum)
    simd_uint8 nlb_before = nlb_counts - per_elem;

    // --- Per-lane predicates computed scalarly (levels are int, can exceed 255) ---
    // smaller_than_p: lane i's level is at or before the priority index pindex,
    //   meaning this string can potentially be incremented by prog.
    // all_filled_after: from lane i onward, every remaining level slot is occupied
    //   (no gaps), so no room to insert new strings after this point.
    alignas(8) uint8_t stp_arr[8] = {};
    alignas(8) uint8_t afa_arr[8] = {};
    for (int i = 0; i < tmp_nlanes; i++) {
        if (tmp_levels[i] <= pindex) stp_arr[i] = 1;
        if ((h - 1 - tmp_levels[i]) == (tmp_nlanes - i)) afa_arr[i] = 1;
    }
    simd_uint8 stp_v, afa_v;
    stp_v.copy_from(stp_arr, stdx::element_aligned);
    afa_v.copy_from(afa_arr, stdx::element_aligned);
    simd_uint8_mask smaller_than_p = has_bits and (stp_v > simd_uint8(0));

    // --- Determine which lanes have no valid successor (are "stuck") ---
    // clear_first_bit = 0xFE: mask to ignore bit position 0 (the leading bit).
    // pattern_zero_and_ones: what the bits look like if leading bit is 0 and
    //   all NLB positions are 1 (i.e., the string is 01^j).
    simd_uint8 clear_first_bit(~simd_uint8{0} & ~1);
    simd_uint8 pattern_zero_and_ones = tmp_masks & clear_first_bit;

    // A lane has no successor when:
    //   1) nlb_before == t: all t NLB slots are already used in earlier lanes, OR
    //   2) nlb_counts == t (all t NLB used up to and including this lane) AND either:
    //      a) bits == 01^j pattern AND all levels after this are filled (can't reset
    //         downward because there's nowhere to put new strings), OR
    //      b) bits == mask, i.e. all positions are 1 (string is 1^j, maximum value)
    simd_uint8_mask no_successor = has_bits and ((nlb_before == t) or ((nlb_counts == t) and (
             ((tmp_bits == pattern_zero_and_ones) and (afa_v > simd_uint8(0)))
          or (tmp_bits == tmp_masks)))
    );

    // A lane has a successor if it's before pindex, non-empty, and not stuck.
    simd_uint8_mask has_successor = smaller_than_p and has_bits and !no_successor;

    // --- Find the rightmost (highest-index) lane that can be incremented ---
    // We search from the right so that the resulting measure is the smallest
    // value strictly greater than the current one (lexicographic successor).
    int match = stdx::find_last_set(has_successor);
    if (match == -1)
    {
        // No lane can be incremented. Two sub-cases:
        if (tmp_levels[0] == 0)
        {
            // Already at the maximum for this level structure => overflow to Top
            tmp_nlanes = 1;
            tmp_levels[0] = -1;
            fill_inactive_tmp();
            return;
        }
        else
        {
            // Shift the first string to an earlier level and reset its bits to "1"
            // (the smallest non-empty bitstring with leading bit 1).
            // Mirrors scalar: set_index = min(levels[0], pindex+1) - 1
            int new_level = std::min(tmp_levels[0], pindex + 1) - 1;
            if (new_level < 0)
            {
                // Overflow to Top
                tmp_nlanes = 1;
                tmp_levels[0] = -1;
                fill_inactive_tmp();
                return;
            }
            tmp_bits[0] = 1;
            match = 0;
            tmp_levels[0] = new_level;
        }
    }
    else if (nlb_counts[match] == t)
    {
        // Case A: All t NLB slots are consumed up to this lane.
        // The string ends with a pattern like ...10...01^j (trailing ones after
        // the last zero). We must "carry": erase the trailing 01^j suffix.

        // countl_one on (bits | ~mask) counts consecutive 1s from the MSB of
        // the used portion — these are the trailing ones in the bitstring.
        // +1 includes the zero bit that precedes them (the "carry" position).
        int reset_bits = std::countl_one(static_cast<uint8_t>(tmp_bits[match] | ~tmp_masks[match])) + 1;
        int current_bits = std::popcount(static_cast<uint8_t>(tmp_masks[match]));
        // Clear the top reset_bits positions in both bits and mask:
        // (1 << (8 - reset_bits)) - 1 produces a mask keeping only the lower bits.
        tmp_bits[match] &= (1u << (8-reset_bits)) - 1;
        tmp_masks[match] &= (1u << (8-reset_bits)) - 1;
        if (reset_bits < 8)
        {
            // Partial reset: some bits remain in this lane. Update NLB count
            // and move to the next lane for the tail fill.
            nlb_counts[match] -= current_bits - std::popcount(static_cast<uint8_t>(tmp_masks[match]));
            match++;
            if (match < k-1)
            {
                tmp_levels[match] = tmp_levels[match-1] + 1;
                tmp_bits[match] = 0;
            }
        }
        else
        {
            // Full reset: entire string in this lane was erased.
            // Keep the lane but move it to the next available level.
            nlb_counts[match] -= (current_bits - 1);
            if ((match + 1 == tmp_nlanes and tmp_levels[match] < h-2) or
                (match + 1 < tmp_nlanes and tmp_levels[match] + 1 < tmp_levels[match + 1]) )
            {
                // There's a gap — place this empty string at the next level.
                tmp_levels[match] = (match + 1 == tmp_nlanes) ? tmp_levels[match] + 1 : tmp_levels[match + 1] - 1;
            }
            else tmp_levels[match] = tmp_levels[match] + 1;
        }
    }
    else
    {
        // Case B: There are still unused NLB slots available.
        // We can grow the measure by appending a bit.
        if ((match + 1 < tmp_nlanes and tmp_levels[match] + 1 < tmp_levels[match + 1] and tmp_levels[match] != pindex))
        {
            // There's an empty level between this lane and the next —
            // start a new string there with leading bit 1.
            match++;
            tmp_bits[match] = 1;
            tmp_levels[match] = std::min(pindex, tmp_levels[match] - 1);
        }
        else
        {
            // Extend the current string by setting the next bit position.
            // first_new = current string length = popcount(mask).
            int first_new = std::popcount(static_cast<uint8_t>(tmp_masks[match]));
            tmp_bits[match] |= (1u << first_new);
        }
    }
    if (match < k-1)
    {
        // --- Fill the tail: set all lanes after 'match' to minimum values ---
        // Compute how many NLB are left to distribute, then set this lane's
        // mask to use exactly (t - bits_before + 1) positions (leading bit + NLBs).
        int bits_before = match > 0 ? nlb_counts[match - 1] : 0;
        tmp_masks[match] = (1u << (t - bits_before + 1)) - 1;

        // Append new single-bit strings at consecutive levels after 'match'.
        tmp_nlanes = match + 1;
        int set_to_level = tmp_levels[match] + 1;
        while (tmp_nlanes < k-1 and set_to_level <= h-2)
        {
            tmp_levels[tmp_nlanes] = set_to_level;
            tmp_nlanes++;
            set_to_level++;
        }
        // SIMD bulk-set: for all lanes in (match, tmp_nlanes), set mask=1, bits=0.
        // Lanes beyond tmp_nlanes are zeroed by fill_inactive_tmp below.
        simd_uint8_mask after_match_mask   = LANE_INDICES > simd_uint8(static_cast<uint8_t>(match));
        simd_uint8_mask before_levels_mask = LANE_INDICES < simd_uint8(tmp_nlanes);
        stdx::where(before_levels_mask and after_match_mask, tmp_masks) = 1;
        stdx::where(after_match_mask, tmp_bits) = 0;

        fill_inactive_tmp();
    }
}

/**
 * Write pm to ostream.
 */
void
STRPM_SIMDSolver::stream_pm(std::ostream &out, int idx)
{
    uint8_t nlanes = pm_nlanes[idx];
    if (nlanes > 0 and pm_levels[idx*8] == -1) {
        out << " \033[1;33mTop\033[m";
    } else {
        out << " { ";
        int output_level = 0;
        for (int i = 0; i < nlanes; i++) {
            if (i>0) out << ", ";
            while (pm_levels[idx*8 + i] > output_level)
            {
                out << "ε, ";
                output_level++;
            }
            out << std::bitset<8>(pm_bits[idx*8 + i]) << "/" << std::bitset<8>(pm_masks[idx*8 + i]);
        }
        while (output_level < h-2)
        {
            out << ", ε";
            output_level++;
        }
        out << " }";
    }
}

/**
 * Write SIMD state to ostream.
 */
void
STRPM_SIMDSolver::stream_simd(std::ostream &out, simd_uint8& bits, simd_uint8& masks, int* levels, uint8_t nlanes)
{
    if (nlanes > 0 and levels[0] == -1) {
        out << " \033[1;33mTop\033[m";
    } else {
        out << " { ";
        int output_level = 0;
        for (int i = 0; i < nlanes; i++) {
            if (i>0) out << ", ";
            while (levels[i] > output_level)
            {
                out << "ε, ";
                output_level++;
            }
            out << std::bitset<8>(static_cast<uint8_t>(bits[i])) << "/" << std::bitset<8>(static_cast<uint8_t>(masks[i]));
        }
        while (output_level < h-2)
        {
            out << ", ε";
            output_level++;
        }
        out << " }";
    }
}

/**
 * Compare tmp and best
 * res := -1 :: tmp < best
 * res := 0  :: tmp = best
 * res := 1  :: tmp > best
 */
int
STRPM_SIMDSolver::compare(int pindex)
{
    if (k == 1)
    {
        // It is either empty or Top, so comparing sizes is enough
        if (tmp_nlanes == best_nlanes) return 0;
        else if (tmp_nlanes > 0) return 1;
        else if (best_nlanes > 0) return -1;

    }

    // cases involving Top
    if (tmp_levels[0] == -1 and best_levels[0] == -1) return 0;
    if (tmp_levels[0] == -1) return 1;
    if (best_levels[0] == -1) return -1;

    // --- SIMD precomputation: compare all 8 bitstrings in parallel ---
    // Each bitstring is encoded as (bits, mask) where mask indicates which
    // positions exist. A longer string (more mask bits) is compared by first
    // checking the shared prefix, then the extra bit decides ordering.
    //
    // shorter_string: intersection of masks = positions present in both strings.
    // diff: positions that exist in one string but not the other.
    // first_length_difference: isolate the lowest such position (where lengths diverge).
    //   diff & ~(diff << 1) clears all but the least-significant 1 in each
    //   contiguous run of 1s in diff, giving us the first position where one
    //   string is longer than the other.
    simd_uint8 shorter_string = tmp_masks & best_masks;
    simd_uint8 diff = tmp_masks ^ best_masks;
    simd_uint8 first_length_difference = diff & ~(diff << 1);
    simd_uint8 bit_xor = tmp_bits ^ best_bits;
    // combined = shared positions + the first extra position.
    // relevant_xor = bit differences within this relevant region.
    simd_uint8 combined     = shorter_string + first_length_difference;
    simd_uint8 relevant_xor = bit_xor & combined;
    // a_less: tmp's string is shorter (fewer mask bits) and the extra bit in
    //   best is different from tmp (relevant_xor matches the length diff), OR
    //   tmp is longer but all shared bits are identical (the extra 0-extension wins).
    simd_uint8_mask a_less = ((tmp_masks < best_masks) and (relevant_xor == first_length_difference)) or
                             ((tmp_masks > best_masks) and (relevant_xor == 0));
    simd_uint8_mask b_less = ((best_masks < tmp_masks) and (relevant_xor == first_length_difference)) or
                             ((best_masks > tmp_masks) and (relevant_xor == 0));

    // For strings of equal length (or within the shared prefix), find the
    // first bit position where they differ and check who has the 1.
    // different_bits: positions in the shared region where bits disagree.
    // Isolate the lowest such bit with x & (-x) (two's complement trick).
    simd_uint8 different_bits = (shorter_string & bit_xor);
    simd_uint8 first_bit_difference = different_bits & (simd_uint8(0) - different_bits);
    simd_uint8_mask a_greater = (different_bits > 0) and ((tmp_bits & first_bit_difference) > 0);
    simd_uint8_mask b_greater = (different_bits > 0) and ((best_bits & first_bit_difference) > 0);

    // --- Scalar level loop ---
    // Levels are int (can exceed 255), so we iterate scalarly over lanes.
    // The SIMD-precomputed a_greater/b_greater/a_less/b_less masks tell us
    // the bitstring comparison result for each lane at O(1) cost.
    uint8_t max_nlanes = std::max(tmp_nlanes, best_nlanes);
    for (uint8_t i = 0; i < max_nlanes; i++)
    {
        // One of the two has more nonempty strings but we were equal up until here
        if (i >= best_nlanes)
        {
            return (tmp_bits[i] & 1) == 0 ? -1 : 1;
        }
        else if (i >= tmp_nlanes)
        {
            return (best_bits[i] & 1) == 0 ? 1 : -1;
        }
        int tl = tmp_levels[i];
        int bl = best_levels[i];
        // We are past the considered index!
        if (tl > pindex and bl > pindex)
        {
            break;
        }
        // One of the two has a bitstring "earlier"
        else if ((tl <= pindex and bl > pindex) or (tl < bl))
        {
            return (tmp_bits[i] & 1) == 0 ? -1 : 1;
        }
        else if ((tl > pindex and bl <= pindex) or (tl > bl))
        {
            return (best_bits[i] & 1) == 0 ? 1 : -1;
        }
        // The two levels are equal... We have to compare strings
        else if (a_greater[i] or (!a_less[i] and b_less[i]))
        {
            return 1;
        }
        else if (b_greater[i] or (a_less[i] and !b_less[i]))
        {
            return -1;
        }
    }

    // The bitstrings are entirely equal
    return 0;
}

bool
STRPM_SIMDSolver::lift(int v, int target, int &str, int pl)
{
    // check if already Top
    if (pm_nlanes[v] > 0 and pm_levels[v*8] == -1) return false; // already Top

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
                stream_simd(logger, tmp_bits, tmp_masks, tmp_levels, tmp_nlanes);
                logger << " =>";
            }
#endif
        if (pl == (pr&1)) prog_tmp(pindex, h);
        //else trunc_tmp(pindex);
#ifndef NDEBUG
            if (trace >= 2) {
                stream_simd(logger, tmp_bits, tmp_masks, tmp_levels, tmp_nlanes);
                logger << std::endl;
            }
#endif
        to_best(v);
        if (compare(pindex) > 0) {
            from_tmp(v);
#ifndef NDEBUG
            if (trace >= 1) {
                logger << "\033[32;1mnew measure\033[m of \033[36;1m" << label_vertex(v) << "\033[m:";
                stream_simd(logger, tmp_bits, tmp_masks, tmp_levels, tmp_nlanes);
                logger << " (to " << label_vertex(target) << ")\n";
            }
#endif
            return true;
        } else {
            return false;
        }
    }

    // Pre-collect valid (non-disabled) successors to avoid branch
    // mispredictions from the disabled check inside the hot loop.
    // succs is a member vector, pre-reserved to nodecount() in run().
    succs.clear();
    for (auto curedge = outs(v); *curedge != -1; curedge++) {
        int to = *curedge;
        if (!disabled[to]) succs.push_back(to);
    }
    int nsuccs = static_cast<int>(succs.size());

    const bool do_prog = (pl == (pr&1));
    const bool want_max = (owner(v) == pl);

    bool first = true;

    // Use the interleaved batch path when batch_B >= 2 (and we have successors).
    // batch_prog_tmp computes prog for batch_B PMs simultaneously using AVX2-wide
    // SIMD; we then extract each result and compare against the running best.
    if (batch_B >= 2 && nsuccs > 0) {
        const int stride = std::max(k - 1, 1);
        for (int si = 0; si < nsuccs; si += batch_B) {
            int cnt = std::min(batch_B, nsuccs - si);
            load_batch(succs.data() + si, cnt);
            if (do_prog) batch_prog_tmp(pindex, h, cnt);

            for (int j = 0; j < cnt; j++) {
                // Early exit for want_max: Top is unsurpassable
                if (want_max && batch_nlanes[j] > 0 &&
                        batch_levels[j * stride] == -1) {
                    extract_batch_to_best(j);
                    str = batch_ids[j];
                    goto batch_done;
                }
                // Fast first-string filter: skip extract+compare if j clearly
                // cannot improve the running best.
                if (!first && batch_can_skip(j, want_max, pindex)) continue;
                extract_batch_to_tmp(j);
#ifndef NDEBUG
                if (trace >= 2) {
                    logger << "to successor " << label_vertex(batch_ids[j]) << " (batch) ";
                    stream_simd(logger, tmp_bits, tmp_masks, tmp_levels, tmp_nlanes);
                    logger << std::endl;
                }
#endif
                if (first) {
                    tmp_to_best();
                    str = batch_ids[j];
                    first = false;
                    if (want_max && best_nlanes > 0 && best_levels[0] == -1)
                        goto batch_done;
                } else if (want_max) {
                    if (compare(pindex) > 0) { tmp_to_best(); str = batch_ids[j]; }
                } else {
                    if (compare(pindex) < 0) { tmp_to_best(); str = batch_ids[j]; }
                }
            }
        }
        batch_done:;
    } else {
        // Scalar fallback: batch_B < 2 (k very large) or no successors
        for (int si = 0; si < nsuccs; si++) {
            int to = succs[si];
            to_tmp(to);
#ifndef NDEBUG
            if (trace >= 2) {
                logger << "to successor " << label_vertex(to) << " from";
                stream_simd(logger, tmp_bits, tmp_masks, tmp_levels, tmp_nlanes);
                logger << " =>";
            }
#endif
            if (do_prog) prog_tmp(pindex, h);
#ifndef NDEBUG
            if (trace >= 2) {
                stream_simd(logger, tmp_bits, tmp_masks, tmp_levels, tmp_nlanes);
                logger << std::endl;
            }
#endif
            if (first) {
                tmp_to_best();
                str = to;
                if (want_max and best_nlanes > 0 and best_levels[0] == -1) break;
            } else if (want_max) {
                if (tmp_nlanes > 0 and tmp_levels[0] == -1) {
                    tmp_to_best(); str = to; break;
                }
                if (compare(pindex) > 0) { tmp_to_best(); str = to; }
            } else {
                if (compare(pindex) < 0) { tmp_to_best(); str = to; }
            }
            first = false;
        }
    }

    // set best to pm if higher
    to_tmp(v);
    if (compare(pindex) < 0) {
#ifndef NDEBUG
        if (trace >= 1) {
            logger << "\033[1;32mnew measure\033[m of \033[36;1m" << label_vertex(v) << "\033[m:";
            stream_simd(logger, best_bits, best_masks, best_levels, best_nlanes);
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
                   boost::hash<std::tuple<int, int, int>>> treeU_simd;
std::unordered_map<std::tuple<int, int, int>, unsigned,
                   boost::hash<std::tuple<int, int, int>>> treeV_simd;

unsigned tree_size_simd(int k, int t, int h)
{
    std::stack<Node> stack;

    stack.push ({k, t, h, true});
    while (!stack.empty())
    {
        Node& tos = stack.top();
        if (tos.isU and tos.h == 1 and tos.k == 1)
        {
            treeU_simd[std::make_tuple(tos.k, tos.t, tos.h)] = 1;
            stack.pop ();
        }
        else if (tos.isU and tos.h > 1 and tos.k == 1)
        {
            auto son = treeU_simd.find(std::make_tuple(tos.k, tos.t, tos.h - 1));
            if (son != treeU_simd.end())
            {
                treeU_simd[std::make_tuple(tos.k, tos.t, tos.h)] = son->second;
                stack.pop ();
            }
            else stack.push ({tos.k, tos.t, tos.h - 1, true});
        }
        else if (tos.h >= tos.k and tos.k >= 2 and tos.t == 0)
        {
            auto son = treeU_simd.find(std::make_tuple(tos.k - 1, tos.t, tos.h - 1));
            if (son != treeU_simd.end())
            {
                if (tos.isU) treeU_simd[std::make_tuple(tos.k, tos.t, tos.h)] = son->second;
                else treeV_simd[std::make_tuple(tos.k, tos.t, tos.h)] = son->second;
                stack.pop ();
            }
            else stack.push ({tos.k - 1, tos.t, tos.h - 1, true});
        }
        else if (!tos.isU and tos.h >= tos.k and tos.k >= 2 and tos.t >= 1)
        {
            auto son1 = treeV_simd.find(std::make_tuple(tos.k, tos.t - 1, tos. h));
            auto son2 = treeU_simd.find(std::make_tuple(tos.k - 1, tos.t, tos.h - 1));
            if (son1 != treeV_simd.end() and son2 != treeU_simd.end())
            {
                treeV_simd[std::make_tuple(tos.k, tos.t, tos.h)] = son1->second * 2 + son2->second;
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
            auto son = treeV_simd.find(std::make_tuple(tos.k, tos.t, tos.h));
            if (son != treeV_simd.end())
            {
                treeU_simd[std::make_tuple(tos.k, tos.t, tos.h)] = son->second;
                stack.pop ();
            }
            else stack.push ({tos.k, tos.t, tos.h, false});
        }
        else if (tos.isU and tos.h > tos.k and tos.k >= 2)
        {
            auto son1 = treeV_simd.find(std::make_tuple(tos.k, tos.t, tos.h));
            auto son2 = treeU_simd.find(std::make_tuple(tos.k, tos.t, tos.h - 1));
            if (son1 != treeV_simd.end() and son2 != treeU_simd.end())
            {
                treeU_simd[std::make_tuple(tos.k, tos.t, tos.h)] = son1->second * 2 + son2->second;
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

    return treeU_simd[std::make_tuple(k, t, h)];
}

struct SizeCompare
{
    int h;

    bool operator()(const std::pair<int,int>& lhs,
                    const std::pair<int,int>& rhs) const
    {
        auto lhs_size = tree_size_simd(lhs.first, lhs.second, h);
        auto rhs_size = tree_size_simd(rhs.first, rhs.second, h);

        return lhs_size > rhs_size;
    }
};


void
STRPM_SIMDSolver::run(int t_val, int k_val, int depth, int player)
{
    // Marcin's word: think of h as the number of priorities of the
    // opponent... PLUS ONE!
    t = t_val;
    h = depth + 1;  // FIXME: This is Guillermo's hack, the +1
    k = k_val;  // Maybe possible: std::min(t + 2, h);

    // Initialize interleaved batch registers for successor processing.
    // batch_B PMs are packed into one 32-lane register; PM j occupies
    // lanes j*(k-1) .. j*(k-1)+(k-2).
    {
        int stride = std::max(k - 1, 1);
        batch_B = 32 / stride;
        batch_levels.assign(batch_B * stride, 0);
        batch_nlanes.assign(batch_B, 0);
        batch_ids   .assign(batch_B, -1);
        batch_match .assign(batch_B, -1);
    }

#ifndef NDEBUG
    logger << "Strahler-tree parameters for player " << player << ": k = " << k << ", t = " << t << ", h = " << h << std::endl;
#endif

    // Initialize progress measures using flat arrays.
    // Every node is set to the smallest leaf in the tree.
    const int nc = nodecount();
    pm_bits.assign(nc * 8, 0);
    pm_masks.assign(nc * 8, 0);
    pm_levels.assign(nc * 8, 0);
    pm_nlanes.assign(nc, static_cast<uint8_t>(k - 1));

    // Build the initial mask and levels pattern (once), then stamp it to every node.
    alignas(8) uint8_t initial_mask[8] = {};
    int initial_levels[8] = {};

    initial_mask[0] = (1 << (t+1)) - 1;
    initial_levels[0] = 0;
    for (int i = 1; i < k-1; i++)
    {
        initial_levels[i] = i;
        initial_mask[i] = 1;
    }
    for (int n = 0; n < nc; n++)
    {
        std::memcpy(&pm_masks[n * 8], initial_mask, 8);
        std::memcpy(&pm_levels[n * 8], initial_levels, 8 * sizeof(int));
    }

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
        if (pm_nlanes[v] == 0 or pm_levels[v*8] != -1) {
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

            if (pm_nlanes[v] == 0 or pm_levels[v*8] != -1) {
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
        if (pm_nlanes[v] == 0 or pm_levels[v*8] != -1) Solver::solve(v, 1-player, game.getStrategy(v));
    }

    Solver::flush();
}

void
STRPM_SIMDSolver::run()
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
    succs.reserve(nodecount());

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
