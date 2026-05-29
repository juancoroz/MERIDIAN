//  sobol_runner.cpp  --  Sobol sensitivity-index estimation
//
//  Implements the pick-freeze estimator (Saltelli 2010) and Jansen's
//  total-effect estimator on top of the existing rocket_6dof MC
//  infrastructure.  Reuses:
//    - json parsing & set_path mutation
//    - Distribution + quantile() for inverse-CDF sampling
//    - run_single_mission() for the per-run lifecycle
//
//  For k input variables and base sample size N, runs N*(k+2)
//  simulations.  Output: per-(output, input) S1 and ST estimates.
//
//  Estimators (with E_f0² = (1/N) Σ A_j * B_j):
//    V(Y) = (1/N) Σ A_j² - ((1/N) Σ A_j)²
//    S_i  = (1/N) Σ B_j (AB_i_j - A_j) / V(Y)         [Saltelli 2010]
//    ST_i = (1/(2N)) Σ (A_j - AB_i_j)^2 / V(Y)        [Jansen 1999]
//
//  S_i can come out slightly negative for small N when the input
//  has essentially no effect (the estimator has variance even at the
//  null).  We report the raw value -- negative S_i indicates
//  "indistinguishable from zero given this sample size".

#include "sobol_runner.h"
#include "json.h"
#include "distributions.h"
#include "sim_runner.h"
#include "parallel_runner.h"
#include "propulsion.h"

#include <fstream>
#include <regex>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <map>
#include <set>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

namespace rocket6dof {

namespace {

// ---- Per-output Sobol indices ----
struct Indices {
    std::vector<double> S1;    // first-order, length k
    std::vector<double> ST;    // total-effect, length k
    double mean;
    double var;
};

// ---- Per-(output, input) bootstrap CI ----
struct CI {
    std::vector<double> S1_lo, S1_hi;  // 2.5% / 97.5% percentiles, length k
    std::vector<double> ST_lo, ST_hi;
};

// Bootstrap CI for second-order indices S_ij.
// Indexed by pair position p (matches the pair_ij list in run_sobol).
struct CI_pairs {
    std::vector<double> Sij_lo, Sij_hi;  // length n_pairs
};

// Bootstrap CI for third-order indices S_ijk.
// Indexed by triplet position t (matches the trip_ijk list in run_sobol).
struct CI_triplets {
    std::vector<double> Sijk_lo, Sijk_hi;  // length n_triplets
};

// Bootstrap CI for group-pair indices S_{Ga,Gb}.
// Indexed by group-pair position p (matches the group_pair_ij list).
struct CI_gpairs {
    std::vector<double> SGij_lo, SGij_hi;  // length n_group_pairs
};

// Bootstrap CI for the per-input variance share
//   share[i] = S1[i] + ½·Σ_j S_ij[i,j] + ⅓·Σ_jk S_ijk[i,j,k]
// computed jointly inside each bootstrap resample.  This replaces the
// per-term independent-sum approximation (which overshoots because the
// S1/S_ij/S_ijk estimators share Y_A/Y_B samples and are positively
// correlated).
struct CI_share {
    std::vector<double> share_lo, share_hi;  // length k (per input)
};

// Saltelli/Jansen estimators given Y_A, Y_B, and Y_AB (one column per
// input).  All vectors have length N; Y_AB is k * N flattened
// (Y_AB[i*N + j] = j-th sample of the AB_i matrix).
Indices compute_indices(const std::vector<double>& YA,
                        const std::vector<double>& YB,
                        const std::vector<std::vector<double>>& YAB)
{
    const int N = static_cast<int>(YA.size());
    const int k = static_cast<int>(YAB.size());
    Indices out;
    out.S1.assign(k, 0.0);
    out.ST.assign(k, 0.0);

    // Combined-sample mean and variance: use all 2N samples for
    // robustness (recommended in Saltelli 2010 Appendix B).
    // NaN-robust: count only finite samples so failed simulations
    // don't propagate into the point estimates.  Same fix pattern
    // as bootstrap_indices below.
    double sumA = 0, sumB = 0;
    double sumA2 = 0;
    int n_valid_A = 0, n_valid_B = 0, n_valid_A2 = 0;
    for (int j = 0; j < N; ++j) {
        if (std::isfinite(YA[j])) { sumA += YA[j]; sumA2 += YA[j]*YA[j]; ++n_valid_A; ++n_valid_A2; }
        if (std::isfinite(YB[j])) { sumB += YB[j]; ++n_valid_B; }
    }
    if (n_valid_A < 2) {
        out.mean = 0.0;
        out.var  = 0.0;
        return out;
    }
    double mean = (sumA + sumB) / (double)(n_valid_A + n_valid_B);
    double meanA = sumA / n_valid_A;
    double varA = sumA2 / n_valid_A - meanA * meanA;  // variance estimator from A alone
    out.mean = mean;
    out.var  = varA;

    if (!(varA > 0) || !std::isfinite(varA)) {
        // No variance in output -- sensitivities undefined.  Return zeros.
        return out;
    }

    for (int i = 0; i < k; ++i) {
        const auto& YABi = YAB[i];
        double s_num = 0, t_num = 0;
        int n_pairs = 0;
        for (int j = 0; j < N; ++j) {
            // Skip if any of the three values are NaN; we need
            // all three to form a valid pair contribution.
            if (!std::isfinite(YA[j]) || !std::isfinite(YB[j])
                || !std::isfinite(YABi[j])) continue;
            // Saltelli 2010, eqn (b): S_i numerator
            s_num += YB[j] * (YABi[j] - YA[j]);
            // Jansen 1999: ST_i numerator
            double d = YA[j] - YABi[j];
            t_num += d * d;
            ++n_pairs;
        }
        if (n_pairs < 2) {
            out.S1[i] = 0.0;
            out.ST[i] = 0.0;
        } else {
            out.S1[i] = (s_num / n_pairs) / varA;
            out.ST[i] = (t_num / (2.0 * n_pairs)) / varA;
        }
    }
    return out;
}

// Bootstrap CI for a single output: resample (with replacement) the N
// sample indices B times, recompute indices on each resample, and
// return the percentile range across resamples.  Uses the same RNG
// stream as the main run so results are reproducible given the seed.
CI bootstrap_indices(const std::vector<double>& YA,
                     const std::vector<double>& YB,
                     const std::vector<std::vector<double>>& YAB,
                     int B, std::mt19937& rng,
                     double ci_level = 0.95)
{
    const int N = static_cast<int>(YA.size());
    const int k = static_cast<int>(YAB.size());

    // ---- Precompute finite masks ----
    // The inner bootstrap loop checks isfinite() on YA[q], YB[q], and
    // YAB[i][q] for every resampled index, three times.  This is hot:
    // ~B * k * N * 3 isfinite() calls per output.  Precompute a
    // per-input "all three are finite" mask so the inner loop becomes
    // a single byte read.  On B=500, k=10, N=128 (typical first-order),
    // this saves ~2M isfinite() calls per output.
    //
    // Also precompute the finite-A subset for the YA-variance step
    // (which doesn't involve YB or YAB[i]).
    std::vector<char> mask_A(N);
    for (int j = 0; j < N; ++j) mask_A[j] = std::isfinite(YA[j]) ? 1 : 0;
    // Per-input mask: 1 if YA[j] && YB[j] && YAB[i][j] all finite
    std::vector<std::vector<char>> mask(k, std::vector<char>(N));
    for (int i = 0; i < k; ++i) {
        const auto& YABi = YAB[i];
        for (int j = 0; j < N; ++j) {
            mask[i][j] = (mask_A[j] && std::isfinite(YB[j])
                          && std::isfinite(YABi[j])) ? 1 : 0;
        }
    }

    // Storage: B resamples of S1[i] and ST[i] for each input i
    std::vector<std::vector<double>> S1_samples(k), ST_samples(k);
    for (int i = 0; i < k; ++i) {
        S1_samples[i].reserve(B);
        ST_samples[i].reserve(B);
    }

    std::uniform_int_distribution<int> pick(0, N - 1);
    std::vector<int> idx(N);

    for (int b = 0; b < B; ++b) {
        // Resample indices
        for (int j = 0; j < N; ++j) idx[j] = pick(rng);

        // Recompute mean and variance on the resample, skipping
        // any NaN values in YA.  Uses precomputed mask_A.
        double sumA = 0, sumA2 = 0;
        int n_valid_A = 0;
        for (int j = 0; j < N; ++j) {
            int q = idx[j];
            if (!mask_A[q]) continue;
            double y = YA[q];
            sumA  += y;
            sumA2 += y * y;
            ++n_valid_A;
        }
        if (n_valid_A < 2) {
            for (int i = 0; i < k; ++i) {
                S1_samples[i].push_back(0.0);
                ST_samples[i].push_back(0.0);
            }
            continue;
        }
        double meanA = sumA / n_valid_A;
        double varA  = sumA2 / n_valid_A - meanA * meanA;
        if (!(varA > 0) || !std::isfinite(varA)) {
            // Degenerate bootstrap; record zeros for this resample
            for (int i = 0; i < k; ++i) {
                S1_samples[i].push_back(0.0);
                ST_samples[i].push_back(0.0);
            }
            continue;
        }

        for (int i = 0; i < k; ++i) {
            const auto& YABi = YAB[i];
            const auto& maski = mask[i];
            double s_num = 0, t_num = 0;
            int n_pairs = 0;
            for (int j = 0; j < N; ++j) {
                int q = idx[j];
                if (!maski[q]) continue;  // single byte check vs 3 isfinite()
                s_num += YB[q] * (YABi[q] - YA[q]);
                double d = YA[q] - YABi[q];
                t_num += d * d;
                ++n_pairs;
            }
            if (n_pairs < 2) {
                S1_samples[i].push_back(0.0);
                ST_samples[i].push_back(0.0);
            } else {
                S1_samples[i].push_back((s_num / n_pairs) / varA);
                ST_samples[i].push_back((t_num / (2.0 * n_pairs)) / varA);
            }
        }
    }

    // Compute 2.5% / 97.5% percentiles per input.  Filter NaN out of
    // each sample series before sorting (std::sort on NaN values is
    // technically undefined; nan_to_finite drop + finite count gives
    // a well-defined percentile on the clean portion).
    CI ci;
    ci.S1_lo.assign(k, 0); ci.S1_hi.assign(k, 0);
    ci.ST_lo.assign(k, 0); ci.ST_hi.assign(k, 0);
    auto pct = [](std::vector<double>& v, double p) {
        // Drop NaN/inf entries first
        v.erase(std::remove_if(v.begin(), v.end(),
                               [](double x){ return !std::isfinite(x); }),
                v.end());
        if (v.empty()) return std::nan("");
        std::sort(v.begin(), v.end());
        int n = static_cast<int>(v.size());
        double idx = p * (n - 1);
        int lo = static_cast<int>(idx);
        int hi = std::min(lo + 1, n - 1);
        double frac = idx - lo;
        return v[lo] * (1 - frac) + v[hi] * frac;
    };
    const double p_lo = 0.5 * (1.0 - ci_level);
    const double p_hi = 0.5 * (1.0 + ci_level);
    for (int i = 0; i < k; ++i) {
        ci.S1_lo[i] = pct(S1_samples[i], p_lo);
        ci.S1_hi[i] = pct(S1_samples[i], p_hi);
        ci.ST_lo[i] = pct(ST_samples[i], p_lo);
        ci.ST_hi[i] = pct(ST_samples[i], p_hi);
    }
    return ci;
}

// Bootstrap CI for second-order pair indices S_ij.
//   pair_ij[p] = (i, j) with i < j
//   YBA[p][n] = j-th output sample of XBA_p matrix (B with cols i,j from A)
// For each bootstrap resample we recompute on the *same* index set:
//   varA, S1[i], S1[j], and V_closed_ij = (1/N) Σ YA[q] * (YBA_p[q] - YB[q])
//   then S_ij_b = V_closed_ij / varA - S1[i] - S1[j]
// Uses a separate RNG offset so it's reproducible but independent of
// the first-order bootstrap stream.
CI_pairs bootstrap_pair_indices(const std::vector<double>& YA,
                                const std::vector<double>& YB,
                                const std::vector<std::vector<double>>& YAB,
                                const std::vector<std::vector<double>>& YBA,
                                const std::vector<std::pair<int,int>>& pair_ij,
                                int B, std::mt19937& rng,
                                double ci_level = 0.95)
{
    const int N      = static_cast<int>(YA.size());
    const int k      = static_cast<int>(YAB.size());
    const int npairs = static_cast<int>(YBA.size());

    // Precompute finite masks (see bootstrap_indices for rationale).
    // First-order mask per input i: YA && YB && YAB[i] finite at j
    // Pair mask per pair p: YA && YB && YBA[p] finite at j
    std::vector<char> mask_A(N);
    for (int j = 0; j < N; ++j) mask_A[j] = std::isfinite(YA[j]) ? 1 : 0;
    std::vector<std::vector<char>> mask_s1(k, std::vector<char>(N));
    for (int i = 0; i < k; ++i) {
        const auto& YABi = YAB[i];
        for (int j = 0; j < N; ++j) {
            mask_s1[i][j] = (mask_A[j] && std::isfinite(YB[j])
                             && std::isfinite(YABi[j])) ? 1 : 0;
        }
    }
    std::vector<std::vector<char>> mask_pair(npairs, std::vector<char>(N));
    for (int p = 0; p < npairs; ++p) {
        const auto& YBAp = YBA[p];
        for (int j = 0; j < N; ++j) {
            mask_pair[p][j] = (mask_A[j] && std::isfinite(YB[j])
                               && std::isfinite(YBAp[j])) ? 1 : 0;
        }
    }

    std::vector<std::vector<double>> Sij_samples(npairs);
    for (int p = 0; p < npairs; ++p) Sij_samples[p].reserve(B);

    std::uniform_int_distribution<int> pick(0, N - 1);
    std::vector<int> idx(N);
    std::vector<double> S1_local(k, 0.0);  // per-resample S1, reused

    for (int b = 0; b < B; ++b) {
        for (int j = 0; j < N; ++j) idx[j] = pick(rng);

        // Mean and variance on the resample, NaN-robust via mask_A.
        double sumA = 0, sumA2 = 0;
        int n_valid_A = 0;
        for (int j = 0; j < N; ++j) {
            int q = idx[j];
            if (!mask_A[q]) continue;
            double y = YA[q];
            sumA  += y;
            sumA2 += y * y;
            ++n_valid_A;
        }
        if (n_valid_A < 2) {
            for (int p = 0; p < npairs; ++p) Sij_samples[p].push_back(0.0);
            continue;
        }
        double meanA = sumA / n_valid_A;
        double varA  = sumA2 / n_valid_A - meanA * meanA;
        if (!(varA > 0) || !std::isfinite(varA)) {
            for (int p = 0; p < npairs; ++p) Sij_samples[p].push_back(0.0);
            continue;
        }

        // First-order S1 on the resample using mask_s1
        for (int i = 0; i < k; ++i) {
            const auto& YABi = YAB[i];
            const auto& maski = mask_s1[i];
            double s_num = 0;
            int n_pairs_i = 0;
            for (int j = 0; j < N; ++j) {
                int q = idx[j];
                if (!maski[q]) continue;
                s_num += YB[q] * (YABi[q] - YA[q]);
                ++n_pairs_i;
            }
            S1_local[i] = (n_pairs_i >= 2) ? ((s_num / n_pairs_i) / varA) : 0.0;
        }

        // S_ij per pair on the resample using mask_pair
        for (int p = 0; p < npairs; ++p) {
            int i = pair_ij[p].first;
            int j = pair_ij[p].second;
            const auto& YBAp = YBA[p];
            const auto& maskp = mask_pair[p];
            double c_num = 0;
            int n_pairs_p = 0;
            for (int n = 0; n < N; ++n) {
                int q = idx[n];
                if (!maskp[q]) continue;
                c_num += YA[q] * (YBAp[q] - YB[q]);
                ++n_pairs_p;
            }
            if (n_pairs_p < 2) {
                Sij_samples[p].push_back(0.0);
            } else {
                double V_closed = c_num / n_pairs_p;
                double Sij_b = V_closed / varA - S1_local[i] - S1_local[j];
                Sij_samples[p].push_back(Sij_b);
            }
        }
    }

    // Percentile across resamples for each pair, NaN-robust.
    CI_pairs out;
    out.Sij_lo.assign(npairs, 0);
    out.Sij_hi.assign(npairs, 0);
    auto pct = [](std::vector<double>& v, double p) {
        v.erase(std::remove_if(v.begin(), v.end(),
                               [](double x){ return !std::isfinite(x); }),
                v.end());
        if (v.empty()) return std::nan("");
        std::sort(v.begin(), v.end());
        int n = static_cast<int>(v.size());
        double idx = p * (n - 1);
        int lo = static_cast<int>(idx);
        int hi = std::min(lo + 1, n - 1);
        double frac = idx - lo;
        return v[lo] * (1 - frac) + v[hi] * frac;
    };
    const double p_lo = 0.5 * (1.0 - ci_level);
    const double p_hi = 0.5 * (1.0 + ci_level);
    for (int p = 0; p < npairs; ++p) {
        out.Sij_lo[p] = pct(Sij_samples[p], p_lo);
        out.Sij_hi[p] = pct(Sij_samples[p], p_hi);
    }
    return out;
}

// Bootstrap CI for third-order pair indices S_ijk.
//   trip_ijk[t] = {i, j, k_} with i<j<k_
//   YBA_3[t][n] = n-th output sample of XBA_ijk (XB with i,j,k_ cols from A)
// For each bootstrap resample we recompute on the *same* index set:
//   varA, S1[i], S1[j], S1[k_], S_ij, S_ik, S_jk, and V_closed_ijk
// Then:
//   S_ijk_b = V_closed/varA - S1[i] - S1[j] - S1[k_]
//                            - S_ij - S_ik - S_jk
// Uses a separate RNG offset so it's reproducible but independent of
// the first-order and pair bootstrap streams.
CI_triplets bootstrap_triplet_indices(
    const std::vector<double>& YA,
    const std::vector<double>& YB,
    const std::vector<std::vector<double>>& YAB,   // size k, each length N
    const std::vector<std::vector<double>>& YBA,   // size n_pairs, each length N
    const std::vector<std::vector<double>>& YBA_3, // size n_triplets, each length N
    const std::vector<std::pair<int,int>>& pair_ij,
    const std::vector<std::array<int,3>>& trip_ijk,
    int B, std::mt19937& rng,
    double ci_level = 0.95)
{
    const int N        = static_cast<int>(YA.size());
    const int k        = static_cast<int>(YAB.size());
    const int n_pairs  = static_cast<int>(YBA.size());
    const int n_triplets = static_cast<int>(YBA_3.size());

    // Precompute finite masks for the three NaN-skip paths.  Saves
    // ~9 isfinite() calls per (b, i, q) tuple -- significant on
    // triplets where the inner loop is B * n_triplets * N deep.
    std::vector<char> mask_A(N);
    for (int j = 0; j < N; ++j) mask_A[j] = std::isfinite(YA[j]) ? 1 : 0;
    std::vector<std::vector<char>> mask_s1(k, std::vector<char>(N));
    for (int i = 0; i < k; ++i) {
        const auto& YABi = YAB[i];
        for (int j = 0; j < N; ++j) {
            mask_s1[i][j] = (mask_A[j] && std::isfinite(YB[j])
                             && std::isfinite(YABi[j])) ? 1 : 0;
        }
    }
    std::vector<std::vector<char>> mask_pair(n_pairs, std::vector<char>(N));
    for (int p = 0; p < n_pairs; ++p) {
        const auto& YBAp = YBA[p];
        for (int j = 0; j < N; ++j) {
            mask_pair[p][j] = (mask_A[j] && std::isfinite(YB[j])
                               && std::isfinite(YBAp[j])) ? 1 : 0;
        }
    }
    std::vector<std::vector<char>> mask_trip(n_triplets, std::vector<char>(N));
    for (int t = 0; t < n_triplets; ++t) {
        const auto& YBA_t = YBA_3[t];
        for (int j = 0; j < N; ++j) {
            mask_trip[t][j] = (mask_A[j] && std::isfinite(YB[j])
                               && std::isfinite(YBA_t[j])) ? 1 : 0;
        }
    }

    // Pair lookup: pair_index_of(a,b) returns p such that pair_ij[p] = (min,max).
    auto pair_index_of = [&](int a, int b) -> int {
        if (a > b) std::swap(a, b);
        return a * k - a * (a + 1) / 2 + (b - a - 1);
    };

    std::vector<std::vector<double>> Sijk_samples(n_triplets);
    for (int t = 0; t < n_triplets; ++t) Sijk_samples[t].reserve(B);

    std::uniform_int_distribution<int> pick(0, N - 1);
    std::vector<int> idx(N);
    std::vector<double> S1_local(k, 0.0);
    std::vector<double> Sij_local(n_pairs, 0.0);

    for (int b = 0; b < B; ++b) {
        for (int j = 0; j < N; ++j) idx[j] = pick(rng);

        // Variance on the resample using mask_A.
        double sumA = 0, sumA2 = 0;
        int n_valid_A = 0;
        for (int j = 0; j < N; ++j) {
            int q = idx[j];
            if (!mask_A[q]) continue;
            double y = YA[q];
            sumA  += y;
            sumA2 += y * y;
            ++n_valid_A;
        }
        if (n_valid_A < 2) {
            for (int t = 0; t < n_triplets; ++t) Sijk_samples[t].push_back(0.0);
            continue;
        }
        double meanA = sumA / n_valid_A;
        double varA  = sumA2 / n_valid_A - meanA * meanA;
        if (!(varA > 0) || !std::isfinite(varA)) {
            for (int t = 0; t < n_triplets; ++t) Sijk_samples[t].push_back(0.0);
            continue;
        }

        // First-order S1 on the resample using mask_s1.
        for (int i = 0; i < k; ++i) {
            const auto& YABi = YAB[i];
            const auto& maski = mask_s1[i];
            double s_num = 0;
            int n_pairs_i = 0;
            for (int j = 0; j < N; ++j) {
                int q = idx[j];
                if (!maski[q]) continue;
                s_num += YB[q] * (YABi[q] - YA[q]);
                ++n_pairs_i;
            }
            S1_local[i] = (n_pairs_i >= 2) ? ((s_num / n_pairs_i) / varA) : 0.0;
        }

        // Second-order S_ij on the resample using mask_pair.
        for (int p = 0; p < n_pairs; ++p) {
            int i = pair_ij[p].first;
            int j = pair_ij[p].second;
            const auto& YBAp = YBA[p];
            const auto& maskp = mask_pair[p];
            double c_num = 0;
            int n_pairs_p = 0;
            for (int n = 0; n < N; ++n) {
                int q = idx[n];
                if (!maskp[q]) continue;
                c_num += YA[q] * (YBAp[q] - YB[q]);
                ++n_pairs_p;
            }
            if (n_pairs_p < 2) {
                Sij_local[p] = 0.0;
            } else {
                double V_closed = c_num / n_pairs_p;
                Sij_local[p] = V_closed / varA - S1_local[i] - S1_local[j];
            }
        }

        // Third-order S_ijk on the resample using mask_trip.
        for (int t = 0; t < n_triplets; ++t) {
            int i  = trip_ijk[t][0];
            int j  = trip_ijk[t][1];
            int k_ = trip_ijk[t][2];
            const auto& YBA_t = YBA_3[t];
            const auto& maskt = mask_trip[t];
            double c_num = 0;
            int n_pairs_t = 0;
            for (int n = 0; n < N; ++n) {
                int q = idx[n];
                if (!maskt[q]) continue;
                c_num += YA[q] * (YBA_t[q] - YB[q]);
                ++n_pairs_t;
            }
            if (n_pairs_t < 2) {
                Sijk_samples[t].push_back(0.0);
            } else {
                double V_closed = c_num / n_pairs_t;
                double Sijk_b = V_closed / varA
                              - S1_local[i] - S1_local[j] - S1_local[k_]
                              - Sij_local[pair_index_of(i, j)]
                              - Sij_local[pair_index_of(i, k_)]
                              - Sij_local[pair_index_of(j, k_)];
                Sijk_samples[t].push_back(Sijk_b);
            }
        }
    }

    CI_triplets out;
    out.Sijk_lo.assign(n_triplets, 0);
    out.Sijk_hi.assign(n_triplets, 0);
    auto pct = [](std::vector<double>& v, double p) {
        v.erase(std::remove_if(v.begin(), v.end(),
                               [](double x){ return !std::isfinite(x); }),
                v.end());
        if (v.empty()) return std::nan("");
        std::sort(v.begin(), v.end());
        int n = static_cast<int>(v.size());
        double i = p * (n - 1);
        int lo = static_cast<int>(i);
        int hi = std::min(lo + 1, n - 1);
        double frac = i - lo;
        return v[lo] * (1 - frac) + v[hi] * frac;
    };
    const double p_lo = 0.5 * (1.0 - ci_level);
    const double p_hi = 0.5 * (1.0 + ci_level);
    for (int t = 0; t < n_triplets; ++t) {
        out.Sijk_lo[t] = pct(Sijk_samples[t], p_lo);
        out.Sijk_hi[t] = pct(Sijk_samples[t], p_hi);
    }
    return out;
}

// Bootstrap CI for group-pair indices S_{Ga,Gb}.
//   group_pair_ij[p] = (ga, gb) with ga < gb
//   group_indices[g] = list of input indices belonging to group g
//   YAB_G[g][n]      = sample n of XABG_g matrix (XA with group g's
//                      columns from XB; for group S1 estimator)
//   YBA_G[p][n]      = sample n of XBA_GaGb matrix
// For each resample we recompute V, group S1[Ga]_b and S1[Gb]_b, then
// V_closed_{Ga,Gb}_b = (1/N) Σ Y_A · (Y_BA_GaGb - Y_B)
// S_{Ga,Gb}_b = V_closed/V - S1[Ga]_b - S1[Gb]_b
//
// Group S1 estimator on the resample (matches the run_sobol path):
// S1[Gg] = (1/N) Σ YB[q] * (YABG_g[q] - YA[q]) / V
CI_gpairs bootstrap_gpair_indices(
    const std::vector<double>& YA,
    const std::vector<double>& YB,
    const std::vector<std::vector<double>>& YAB_G, // size n_groups, each length N
    const std::vector<std::vector<double>>& YBA_G, // size n_group_pairs, each length N
    const std::vector<std::pair<int,int>>& group_pair_ij,
    int B, std::mt19937& rng,
    double ci_level = 0.95)
{
    const int N         = static_cast<int>(YA.size());
    const int n_groups  = static_cast<int>(YAB_G.size());
    const int n_gpairs  = static_cast<int>(YBA_G.size());

    // Precompute finite masks
    std::vector<char> mask_A(N);
    for (int j = 0; j < N; ++j) mask_A[j] = std::isfinite(YA[j]) ? 1 : 0;
    std::vector<std::vector<char>> mask_g(n_groups, std::vector<char>(N));
    for (int g = 0; g < n_groups; ++g) {
        const auto& YABg = YAB_G[g];
        for (int j = 0; j < N; ++j) {
            mask_g[g][j] = (mask_A[j] && std::isfinite(YB[j])
                            && std::isfinite(YABg[j])) ? 1 : 0;
        }
    }
    std::vector<std::vector<char>> mask_gp(n_gpairs, std::vector<char>(N));
    for (int p = 0; p < n_gpairs; ++p) {
        const auto& YBA_Gp = YBA_G[p];
        for (int j = 0; j < N; ++j) {
            mask_gp[p][j] = (mask_A[j] && std::isfinite(YB[j])
                             && std::isfinite(YBA_Gp[j])) ? 1 : 0;
        }
    }

    std::vector<std::vector<double>> SGij_samples(n_gpairs);
    for (int p = 0; p < n_gpairs; ++p) SGij_samples[p].reserve(B);

    std::uniform_int_distribution<int> pick(0, N - 1);
    std::vector<int> idx(N);
    std::vector<double> S1G_local(n_groups, 0.0);

    for (int b = 0; b < B; ++b) {
        for (int j = 0; j < N; ++j) idx[j] = pick(rng);

        // Variance on the resample using mask_A.
        double sumA = 0, sumA2 = 0;
        int n_valid_A = 0;
        for (int j = 0; j < N; ++j) {
            int q = idx[j];
            if (!mask_A[q]) continue;
            double y = YA[q];
            sumA  += y;
            sumA2 += y * y;
            ++n_valid_A;
        }
        if (n_valid_A < 2) {
            for (int p = 0; p < n_gpairs; ++p) SGij_samples[p].push_back(0.0);
            continue;
        }
        double meanA = sumA / n_valid_A;
        double varA  = sumA2 / n_valid_A - meanA * meanA;
        if (!(varA > 0) || !std::isfinite(varA)) {
            for (int p = 0; p < n_gpairs; ++p) SGij_samples[p].push_back(0.0);
            continue;
        }

        // Group first-order S1[Gg] on the resample using mask_g.
        for (int g = 0; g < n_groups; ++g) {
            const auto& YABg = YAB_G[g];
            const auto& maskg = mask_g[g];
            double s_num = 0;
            int n_pairs_g = 0;
            for (int j = 0; j < N; ++j) {
                int q = idx[j];
                if (!maskg[q]) continue;
                s_num += YB[q] * (YABg[q] - YA[q]);
                ++n_pairs_g;
            }
            S1G_local[g] = (n_pairs_g >= 2) ? ((s_num / n_pairs_g) / varA) : 0.0;
        }

        // Group-pair S_{Ga,Gb} on the resample using mask_gp.
        for (int p = 0; p < n_gpairs; ++p) {
            int ga = group_pair_ij[p].first;
            int gb = group_pair_ij[p].second;
            const auto& YBA_Gp = YBA_G[p];
            const auto& maskp = mask_gp[p];
            double c_num = 0;
            int n_pairs_p = 0;
            for (int n = 0; n < N; ++n) {
                int q = idx[n];
                if (!maskp[q]) continue;
                c_num += YA[q] * (YBA_Gp[q] - YB[q]);
                ++n_pairs_p;
            }
            if (n_pairs_p < 2) {
                SGij_samples[p].push_back(0.0);
            } else {
                double V_closed = c_num / n_pairs_p;
                double SGij_b = V_closed / varA - S1G_local[ga] - S1G_local[gb];
                SGij_samples[p].push_back(SGij_b);
            }
        }
    }

    CI_gpairs out;
    out.SGij_lo.assign(n_gpairs, 0);
    out.SGij_hi.assign(n_gpairs, 0);
    auto pct = [](std::vector<double>& v, double p) {
        v.erase(std::remove_if(v.begin(), v.end(),
                               [](double x){ return !std::isfinite(x); }),
                v.end());
        if (v.empty()) return std::nan("");
        std::sort(v.begin(), v.end());
        int n = static_cast<int>(v.size());
        double i = p * (n - 1);
        int lo = static_cast<int>(i);
        int hi = std::min(lo + 1, n - 1);
        double frac = i - lo;
        return v[lo] * (1 - frac) + v[hi] * frac;
    };
    const double p_lo = 0.5 * (1.0 - ci_level);
    const double p_hi = 0.5 * (1.0 + ci_level);
    for (int p = 0; p < n_gpairs; ++p) {
        out.SGij_lo[p] = pct(SGij_samples[p], p_lo);
        out.SGij_hi[p] = pct(SGij_samples[p], p_hi);
    }
    return out;
}

//  Joint bootstrap of clean_share[i] = S1[i] + ½·ΣS_ij + ⅓·ΣS_ijk
// Computes the per-input variance share INSIDE each bootstrap resample
// (rather than propagating per-term CIs after the fact).  This captures
// the joint variation of S1/S_ij/S_ijk that comes from sharing Y_A and
// Y_B samples and produces exact CI on clean_share[i] without needing
// any correlation shrinkage.
//
// Optional: pass empty pair_ij/YBA to skip the second-order term, and
// empty trip_ijk_arr/YBA_3 to skip the third-order term.
// Uses RNG offset distinct from the other bootstraps for reproducibility.
CI_share bootstrap_clean_share(
    const std::vector<double>& YA,
    const std::vector<double>& YB,
    const std::vector<std::vector<double>>& YAB,    // size k
    const std::vector<std::vector<double>>& YBA,    // size n_pairs (may be empty)
    const std::vector<std::vector<double>>& YBA_3,  // size n_triplets (may be empty)
    const std::vector<std::pair<int,int>>& pair_ij,
    const std::vector<std::array<int,3>>& trip_ijk,
    int B, std::mt19937& rng,
    double ci_level = 0.95)
{
    const int N        = static_cast<int>(YA.size());
    const int k        = static_cast<int>(YAB.size());
    const int n_pairs  = static_cast<int>(YBA.size());
    const int n_trip   = static_cast<int>(YBA_3.size());
    const bool do_pair = n_pairs > 0 && !pair_ij.empty();
    const bool do_trip = n_trip  > 0 && !trip_ijk.empty();

    // Precompute finite masks
    std::vector<char> mask_A(N);
    for (int j = 0; j < N; ++j) mask_A[j] = std::isfinite(YA[j]) ? 1 : 0;
    std::vector<std::vector<char>> mask_s1(k, std::vector<char>(N));
    for (int i = 0; i < k; ++i) {
        const auto& YABi = YAB[i];
        for (int j = 0; j < N; ++j) {
            mask_s1[i][j] = (mask_A[j] && std::isfinite(YB[j])
                             && std::isfinite(YABi[j])) ? 1 : 0;
        }
    }
    std::vector<std::vector<char>> mask_pair(do_pair ? n_pairs : 0,
                                              std::vector<char>(N));
    if (do_pair) {
        for (int p = 0; p < n_pairs; ++p) {
            const auto& YBAp = YBA[p];
            for (int j = 0; j < N; ++j) {
                mask_pair[p][j] = (mask_A[j] && std::isfinite(YB[j])
                                   && std::isfinite(YBAp[j])) ? 1 : 0;
            }
        }
    }
    std::vector<std::vector<char>> mask_trip(do_trip ? n_trip : 0,
                                              std::vector<char>(N));
    if (do_trip) {
        for (int t = 0; t < n_trip; ++t) {
            const auto& YBA_t = YBA_3[t];
            for (int j = 0; j < N; ++j) {
                mask_trip[t][j] = (mask_A[j] && std::isfinite(YB[j])
                                   && std::isfinite(YBA_t[j])) ? 1 : 0;
            }
        }
    }

    // Pair index lookup (matches existing convention in bootstrap_triplet_indices).
    auto pair_index_of = [&](int a, int b) -> int {
        if (a > b) std::swap(a, b);
        return a * k - a * (a + 1) / 2 + (b - a - 1);
    };

    std::vector<std::vector<double>> share_samples(k);
    for (int i = 0; i < k; ++i) share_samples[i].reserve(B);

    std::uniform_int_distribution<int> pick(0, N - 1);
    std::vector<int> idx(N);
    std::vector<double> S1_local(k, 0.0);
    std::vector<double> Sij_local(do_pair ? n_pairs : 0, 0.0);

    for (int b = 0; b < B; ++b) {
        // ONE index set for all three orders -- this is the crucial
        // difference from running three separate bootstraps and summing.
        for (int j = 0; j < N; ++j) idx[j] = pick(rng);

        // Variance on resample using mask_A.
        double sumA = 0, sumA2 = 0;
        int n_valid_A = 0;
        for (int j = 0; j < N; ++j) {
            int q = idx[j];
            if (!mask_A[q]) continue;
            double y = YA[q];
            sumA  += y;
            sumA2 += y * y;
            ++n_valid_A;
        }
        if (n_valid_A < 2) {
            for (int i = 0; i < k; ++i) share_samples[i].push_back(0.0);
            continue;
        }
        double meanA = sumA / n_valid_A;
        double varA  = sumA2 / n_valid_A - meanA * meanA;
        if (!(varA > 0) || !std::isfinite(varA)) {
            for (int i = 0; i < k; ++i) share_samples[i].push_back(0.0);
            continue;
        }

        // First-order S1 on this resample using mask_s1.
        for (int i = 0; i < k; ++i) {
            const auto& YABi = YAB[i];
            const auto& maski = mask_s1[i];
            double s_num = 0;
            int n_pairs_i = 0;
            for (int j = 0; j < N; ++j) {
                int q = idx[j];
                if (!maski[q]) continue;
                s_num += YB[q] * (YABi[q] - YA[q]);
                ++n_pairs_i;
            }
            S1_local[i] = (n_pairs_i >= 2) ? ((s_num / n_pairs_i) / varA) : 0.0;
        }

        // Second-order S_ij using mask_pair.
        if (do_pair) {
            for (int p = 0; p < n_pairs; ++p) {
                int i = pair_ij[p].first;
                int j = pair_ij[p].second;
                const auto& YBAp = YBA[p];
                const auto& maskp = mask_pair[p];
                double c_num = 0;
                int n_pairs_p = 0;
                for (int n = 0; n < N; ++n) {
                    int q = idx[n];
                    if (!maskp[q]) continue;
                    c_num += YA[q] * (YBAp[q] - YB[q]);
                    ++n_pairs_p;
                }
                if (n_pairs_p < 2) {
                    Sij_local[p] = 0.0;
                } else {
                    double V_closed = c_num / n_pairs_p;
                    Sij_local[p] = V_closed / varA - S1_local[i] - S1_local[j];
                }
            }
        }

        // Build per-input share for this resample.
        std::vector<double> share_b(k, 0.0);
        for (int i = 0; i < k; ++i) share_b[i] = S1_local[i];

        if (do_pair) {
            // Each pair contributes ½·S_ij to BOTH inputs.
            for (int p = 0; p < n_pairs; ++p) {
                int i = pair_ij[p].first;
                int j = pair_ij[p].second;
                double half = 0.5 * Sij_local[p];
                share_b[i] += half;
                share_b[j] += half;
            }
        }

        if (do_trip) {
            // Compute S_ijk on the resample using mask_trip.
            for (int t = 0; t < n_trip; ++t) {
                int i  = trip_ijk[t][0];
                int j  = trip_ijk[t][1];
                int k_ = trip_ijk[t][2];
                const auto& YBA_t = YBA_3[t];
                const auto& maskt = mask_trip[t];
                double c_num = 0;
                int n_pairs_t = 0;
                for (int n = 0; n < N; ++n) {
                    int q = idx[n];
                    if (!maskt[q]) continue;
                    c_num += YA[q] * (YBA_t[q] - YB[q]);
                    ++n_pairs_t;
                }
                if (n_pairs_t < 2) continue;  // skip this triplet contribution
                double V_closed = c_num / n_pairs_t;
                double Sijk_b = V_closed / varA
                              - S1_local[i] - S1_local[j] - S1_local[k_]
                              - Sij_local[pair_index_of(i, j)]
                              - Sij_local[pair_index_of(i, k_)]
                              - Sij_local[pair_index_of(j, k_)];
                double third = Sijk_b / 3.0;
                share_b[i]  += third;
                share_b[j]  += third;
                share_b[k_] += third;
            }
        }

        for (int i = 0; i < k; ++i) share_samples[i].push_back(share_b[i]);
    }

    CI_share out;
    out.share_lo.assign(k, 0);
    out.share_hi.assign(k, 0);
    auto pct = [](std::vector<double>& v, double p) {
        v.erase(std::remove_if(v.begin(), v.end(),
                               [](double x){ return !std::isfinite(x); }),
                v.end());
        if (v.empty()) return std::nan("");
        std::sort(v.begin(), v.end());
        int n = static_cast<int>(v.size());
        double i = p * (n - 1);
        int lo = static_cast<int>(i);
        int hi = std::min(lo + 1, n - 1);
        double frac = i - lo;
        return v[lo] * (1 - frac) + v[hi] * frac;
    };
    const double p_lo = 0.5 * (1.0 - ci_level);
    const double p_hi = 0.5 * (1.0 + ci_level);
    for (int i = 0; i < k; ++i) {
        out.share_lo[i] = pct(share_samples[i], p_lo);
        out.share_hi[i] = pct(share_samples[i], p_hi);
    }
    return out;
}

} // anon


//  Target CI width parsing
//  sobol.target_ci_width can be either:
//    - a scalar number       => applies to ST only (legacy form)
//    - an object             => may contain .ST and/or .S1 keys
//
//  Examples:
//    "target_ci_width": 0.05
//    "target_ci_width": { "ST": 0.10, "S1": 0.05 }
//    "target_ci_width": { "ST": 0.10 }
//
//  A target of 0 (or absent key in the object form) means "not
//  targeted -- ignore this estimator when computing the
//  recommendation".

struct TargetWidths {
    double ST;       // <= 0 means "not targeted" (global default)
    double S1;
    double S_ij;     // second-order target (forces pilot to enable second_order)
    double S_ijk;    // third-order target (forces pilot to enable third_order)
    // Per-output overrides.  If an output name appears in any of the
    // maps below, that target replaces the global default for that
    // output.  Outputs listed in `skip` are excluded from all order
    // scans entirely.
    std::map<std::string, double> per_out_ST;
    std::map<std::string, double> per_out_S1;
    std::map<std::string, double> per_out_S_ij;
    std::map<std::string, double> per_out_S_ijk;
    std::set<std::string> skip;
    bool has_ST()   const { return ST   > 0 || !per_out_ST.empty(); }
    bool has_S1()   const { return S1   > 0 || !per_out_S1.empty(); }
    bool has_S_ij() const { return S_ij > 0 || !per_out_S_ij.empty(); }
    bool has_S_ijk()const { return S_ijk> 0 || !per_out_S_ijk.empty(); }
    bool any() const { return has_ST() || has_S1() || has_S_ij() || has_S_ijk(); }
    // Per-output lookup with fallback to global.  Returns 0 if neither
    // is set OR if the output is in the skip list.
    double ST_for  (const std::string& o) const {
        if (skip.count(o)) return 0.0;
        auto it = per_out_ST.find(o);
        return (it != per_out_ST.end()) ? it->second : ST;
    }
    double S1_for  (const std::string& o) const {
        if (skip.count(o)) return 0.0;
        auto it = per_out_S1.find(o);
        return (it != per_out_S1.end()) ? it->second : S1;
    }
    double Sij_for (const std::string& o) const {
        if (skip.count(o)) return 0.0;
        auto it = per_out_S_ij.find(o);
        return (it != per_out_S_ij.end()) ? it->second : S_ij;
    }
    double Sijk_for(const std::string& o) const {
        if (skip.count(o)) return 0.0;
        auto it = per_out_S_ijk.find(o);
        return (it != per_out_S_ijk.end()) ? it->second : S_ijk;
    }
};

struct PilotResult {
    // For each estimator we track:
    //   worst_*_w : the worst CI half-width across kept rows (for display)
    //   p90_*_w   : the 90th-percentile CI half-width (for display)
    //   worst_*_N : the worst N_needed = pilot_N * (w/target)^2 across kept
    //               rows where the per-row target is > 0.  This is what
    //               drives the recommendation under "worst" strategy.
    //   p90_*_N   : the 90th-percentile of N_needed.  Same purpose, more
    //               robust to single noisy rows.
    //   driver_*_out / driver_*_inp / driver_*_tgt : the (output, input,
    //               target) tuple of the row that produced worst_*_N.
    //               With per-output targets this can differ from the
    //               worst-width row, so we track both.
    // Recommendation defaults to p90 (more robust to single pathological
    // rows); SOBOL_RECOMMEND_FROM=worst restores conservative behavior.
    double worst_ST_w;     std::string worst_ST_out,    worst_ST_inp;
    double worst_S1_w;     std::string worst_S1_out,    worst_S1_inp;
    double worst_Sij_w;    std::string worst_Sij_out,   worst_Sij_inp;
    double worst_Sijk_w;   std::string worst_Sijk_out,  worst_Sijk_inp;
    double p90_ST_w;
    double p90_S1_w;
    double p90_Sij_w;
    double p90_Sijk_w;
    double worst_ST_N, p90_ST_N, median_ST_N;
    double worst_S1_N, p90_S1_N, median_S1_N;
    double worst_Sij_N, p90_Sij_N, median_Sij_N;
    double worst_Sijk_N, p90_Sijk_N, median_Sijk_N;
    std::string driver_ST_out, driver_ST_inp;   double driver_ST_tgt;
    std::string driver_S1_out, driver_S1_inp;   double driver_S1_tgt;
    std::string driver_Sij_out, driver_Sij_inp; double driver_Sij_tgt;
    std::string driver_Sijk_out, driver_Sijk_inp; double driver_Sijk_tgt;
    int    n_rows;         int n_skipped;
    int    n_pair_rows;    int n_trip_rows;
};

static TargetWidths read_targets(const json::Value& sobol_block) {
    TargetWidths t;
    t.ST = 0.0; t.S1 = 0.0; t.S_ij = 0.0; t.S_ijk = 0.0;
    const json::Value& tcw = sobol_block["target_ci_width"];
    if (tcw.isNumber()) {
        // Scalar form -- legacy.  Applies to ST only.
        t.ST = tcw.asNumber(0.0);
    } else if (tcw.isObject()) {
        t.ST    = tcw["ST"   ].asNumber(0.0);
        t.S1    = tcw["S1"   ].asNumber(0.0);
        t.S_ij  = tcw["S_ij" ].asNumber(0.0);
        t.S_ijk = tcw["S_ijk"].asNumber(0.0);
        // Per-output overrides
        const json::Value& po = tcw["per_output"];
        if (po.isObject()) {
            for (const auto& kv : po.obj()) {
                const std::string& out_name = kv.first;
                const json::Value& ov = kv.second;
                if (!ov.isObject()) continue;
                if (ov["skip"].asBool(false)) {
                    t.skip.insert(out_name);
                    continue;
                }
                if (ov["ST"   ].isNumber()) t.per_out_ST   [out_name] = ov["ST"   ].asNumber(0.0);
                if (ov["S1"   ].isNumber()) t.per_out_S1   [out_name] = ov["S1"   ].asNumber(0.0);
                if (ov["S_ij" ].isNumber()) t.per_out_S_ij [out_name] = ov["S_ij" ].asNumber(0.0);
                if (ov["S_ijk"].isNumber()) t.per_out_S_ijk[out_name] = ov["S_ijk"].asNumber(0.0);
            }
        }
    }
    return t;
}

//  Pilot run helper -- shared by --suggest-n and embedded auto-scale
//  Executes one Sobol pass at pilot_N with bootstrap enabled, then
//  scans the resulting indices CSV to find the worst (largest) CI
//  half-width for ST and (separately) for S1.  Rows where the half-
//  width itself exceeds 1.0 are skipped (estimator in unstable regime
//  -- typically circular outputs like raw psi/theta).
//
//  Returns 0 on success.  The caller is responsible for the actual
//  N-scaling math and reporting.

// Linear-interpolated p-th percentile (0 < p < 1).
// Sorts the input vector in place.  Empty vector returns 0.
static double percentile_inplace(std::vector<double>& v, double p) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    double idx = p * (v.size() - 1);
    int lo = static_cast<int>(idx);
    int hi = std::min(lo + 1, (int)v.size() - 1);
    double frac = idx - lo;
    return v[lo] * (1 - frac) + v[hi] * frac;
}

static int run_pilot_and_measure(const std::string& config_path,
                                 int pilot_N,
                                 const TargetWidths& targets,
                                 PilotResult& out)
{
    out = PilotResult{};
    // Read user config raw for text-level patching.
    std::ifstream fin(config_path);
    if (!fin) {
        std::fprintf(stderr, "ERROR: cannot open %s\n", config_path.c_str());
        return 1;
    }
    std::string cfg_text((std::istreambuf_iterator<char>(fin)),
                          std::istreambuf_iterator<char>());
    fin.close();

    // Patch n_base, n_bootstrap, second_order, third_order, group_pairs.
    // Enable orders only if their CI target is set in `targets` --
    // saves sim cost when the user only cares about first-order.
    bool need_2nd = targets.has_S_ij() || targets.has_S_ijk();
    bool need_3rd = targets.has_S_ijk();
    auto patch_int = [](std::string& s, const std::string& key, int newval) {
        std::regex re(std::string("\"") + key + "\"\\s*:\\s*-?[0-9.eE+-]+");
        std::string repl = "\"" + key + "\": " + std::to_string(newval);
        if (std::regex_search(s, re)) {
            s = std::regex_replace(s, re, repl);
        } else {
            // Key not present -- inject into sobol block.
            std::regex sobol_re("\"sobol\"\\s*:\\s*\\{");
            std::string sobol_inj = "\"sobol\": {\n    " + repl + ",";
            s = std::regex_replace(s, sobol_re, sobol_inj,
                                   std::regex_constants::format_first_only);
        }
    };
    auto patch_bool = [](std::string& s, const std::string& key, bool newval) {
        std::regex re(std::string("\"") + key + "\"\\s*:\\s*(true|false)");
        std::string repl = "\"" + key + "\": " + (newval ? "true" : "false");
        if (std::regex_search(s, re)) {
            s = std::regex_replace(s, re, repl);
        } else {
            // Key not present -- inject into sobol block.
            std::regex sobol_re("\"sobol\"\\s*:\\s*\\{");
            std::string sobol_inj = "\"sobol\": {\n    " + repl + ",";
            s = std::regex_replace(s, sobol_re, sobol_inj,
                                   std::regex_constants::format_first_only);
        }
    };
    patch_int(cfg_text, "n_base", pilot_N);
    patch_int(cfg_text, "n_bootstrap", 500);
    patch_bool(cfg_text, "second_order", need_2nd);
    patch_bool(cfg_text, "third_order",  need_3rd);
    patch_bool(cfg_text, "group_pairs",  false);

    char tmp_cfg[256], tmp_csv[256], tmp_txt[256];
    std::snprintf(tmp_cfg, sizeof(tmp_cfg),
                  "/tmp/sobol_pilot_cfg_%d.json", (int)getpid());
    std::snprintf(tmp_csv, sizeof(tmp_csv),
                  "/tmp/sobol_pilot_%d.csv", (int)getpid());
    std::snprintf(tmp_txt, sizeof(tmp_txt),
                  "/tmp/sobol_pilot_%d.txt", (int)getpid());
    {
        std::ofstream fout(tmp_cfg);
        if (!fout) {
            std::fprintf(stderr, "ERROR: cannot write %s\n", tmp_cfg);
            return 1;
        }
        fout << cfg_text;
    }

    // Set the sentinel so the nested run_sobol skips its own embedded
    // auto-scale block.  Without this, --suggest-n triggers a full
    // production-N run inside the pilot.  Restore prior value on exit.
    const char* prev = std::getenv("SOBOL_IN_PILOT");
    setenv("SOBOL_IN_PILOT", "1", 1);
    int rc = run_sobol(tmp_cfg, tmp_csv, tmp_txt);
    if (prev) setenv("SOBOL_IN_PILOT", prev, 1);
    else      unsetenv("SOBOL_IN_PILOT");
    if (rc != 0) return rc;

    // Scan the indices CSV for worst widths.
    std::FILE* fp = std::fopen(tmp_csv, "r");
    if (!fp) {
        std::fprintf(stderr, "ERROR: cannot read pilot CSV %s\n", tmp_csv);
        return 1;
    }
    char line[2048];
    if (!std::fgets(line, sizeof(line), fp)) {
        std::fclose(fp); return 1;
    }
    std::vector<std::string> hdr;
    {
        std::string s = line;
        size_t p = 0;
        while (p < s.size()) {
            size_t e = s.find(',', p);
            if (e == std::string::npos) e = s.size();
            std::string tok = s.substr(p, e - p);
            while (!tok.empty() && (tok.back()=='\n' || tok.back()=='\r'))
                tok.pop_back();
            hdr.push_back(tok);
            p = e + 1;
        }
    }
    int col_out=-1, col_inp=-1, col_var=-1;
    int col_s1lo=-1, col_s1hi=-1, col_stlo=-1, col_sthi=-1;
    for (size_t i = 0; i < hdr.size(); ++i) {
        if (hdr[i] == "output") col_out = i;
        else if (hdr[i] == "input") col_inp = i;
        else if (hdr[i] == "S1_lo") col_s1lo = i;
        else if (hdr[i] == "S1_hi") col_s1hi = i;
        else if (hdr[i] == "ST_lo") col_stlo = i;
        else if (hdr[i] == "ST_hi") col_sthi = i;
        else if (hdr[i] == "output_variance") col_var = i;
    }
    if (col_stlo < 0 || col_sthi < 0 || col_s1lo < 0 || col_s1hi < 0 || col_var < 0) {
        std::fprintf(stderr,
            "ERROR: pilot CSV missing CI columns (bootstrap may have failed)\n");
        std::fclose(fp); return 1;
    }

    // First-order scan: collect ST and S1 half-widths AND per-row
    // N_needed = pilot_N * (w / target_for_this_output)^2.
    // The width vectors are for display; the N_needed vectors drive
    // the recommendation.
    std::vector<double> st_widths, s1_widths;
    std::vector<double> st_N, s1_N;
    while (std::fgets(line, sizeof(line), fp)) {
        std::vector<std::string> tok;
        std::string s = line;
        size_t p = 0;
        while (p < s.size()) {
            size_t e = s.find(',', p);
            if (e == std::string::npos) e = s.size();
            std::string t = s.substr(p, e - p);
            while (!t.empty() && (t.back()=='\n' || t.back()=='\r'))
                t.pop_back();
            tok.push_back(t);
            p = e + 1;
        }
        if ((int)tok.size() <= col_sthi) continue;
        const std::string& out_name = tok[col_out];
        double v = std::atof(tok[col_var].c_str());
        if (!(v > 0)) continue;
        // Skip outputs explicitly excluded by the user.
        if (targets.skip.count(out_name)) { out.n_skipped++; continue; }
        double st_lo = std::atof(tok[col_stlo].c_str());
        double st_hi = std::atof(tok[col_sthi].c_str());
        double s1_lo = std::atof(tok[col_s1lo].c_str());
        double s1_hi = std::atof(tok[col_s1hi].c_str());
        double st_w  = 0.5 * (st_hi - st_lo);
        double s1_w  = 0.5 * (s1_hi - s1_lo);
        // Skip truly degenerate rows (estimator unstable).
        if (st_w > 1.0 || s1_w > 1.0) { out.n_skipped++; continue; }
        out.n_rows++;
        st_widths.push_back(st_w);
        s1_widths.push_back(s1_w);
        // Per-row N_needed using per-output target.  Only contribute to
        // the recommendation if that output has a target for this
        // estimator (>0).  Also track the driving row (max N_needed)
        // separately from the worst-width row -- with per-output
        // targets they can be different.
        double st_tgt = targets.ST_for(out_name);
        double s1_tgt = targets.S1_for(out_name);
        if (st_tgt > 0) {
            double Nn = static_cast<double>(pilot_N) * (st_w/st_tgt) * (st_w/st_tgt);
            st_N.push_back(Nn);
            if (Nn > out.worst_ST_N) {
                out.worst_ST_N = Nn;
                out.driver_ST_out = out_name;
                out.driver_ST_inp = tok[col_inp];
                out.driver_ST_tgt = st_tgt;
            }
        }
        if (s1_tgt > 0) {
            double Nn = static_cast<double>(pilot_N) * (s1_w/s1_tgt) * (s1_w/s1_tgt);
            s1_N.push_back(Nn);
            if (Nn > out.worst_S1_N) {
                out.worst_S1_N = Nn;
                out.driver_S1_out = out_name;
                out.driver_S1_inp = tok[col_inp];
                out.driver_S1_tgt = s1_tgt;
            }
        }
        if (st_w > out.worst_ST_w) {
            out.worst_ST_w   = st_w;
            out.worst_ST_out = out_name;
            out.worst_ST_inp = tok[col_inp];
        }
        if (s1_w > out.worst_S1_w) {
            out.worst_S1_w   = s1_w;
            out.worst_S1_out = out_name;
            out.worst_S1_inp = tok[col_inp];
        }
    }
    std::fclose(fp);

    // Compute p90 of widths (display) and N_needed (recommendation).
    // worst_ST_N / worst_S1_N already populated during the scan above.
    // Note: percentile_inplace sorts the vector, so reuse works fine.
    out.p90_ST_w = percentile_inplace(st_widths, 0.90);
    out.p90_S1_w = percentile_inplace(s1_widths, 0.90);
    out.p90_ST_N    = percentile_inplace(st_N, 0.90);
    out.median_ST_N = percentile_inplace(st_N, 0.50);
    out.p90_S1_N    = percentile_inplace(s1_N, 0.90);
    out.median_S1_N = percentile_inplace(s1_N, 0.50);

    // Pair CSV scanning -- only if second_order was enabled in the pilot.
    if (need_2nd) {
        std::string pairs_path = tmp_csv;
        size_t dot = pairs_path.find_last_of('.');
        if (dot != std::string::npos)
            pairs_path = pairs_path.substr(0, dot) + "_pairs" + pairs_path.substr(dot);
        else
            pairs_path += "_pairs";
        std::FILE* fpp = std::fopen(pairs_path.c_str(), "r");
        if (fpp) {
            char pline[2048];
            if (std::fgets(pline, sizeof(pline), fpp)) {
                std::vector<std::string> phdr;
                {
                    std::string s = pline; size_t p = 0;
                    while (p < s.size()) {
                        size_t e = s.find(',', p);
                        if (e == std::string::npos) e = s.size();
                        std::string t = s.substr(p, e - p);
                        while (!t.empty() && (t.back()=='\n' || t.back()=='\r')) t.pop_back();
                        phdr.push_back(t);
                        p = e + 1;
                    }
                }
                int p_col_out=-1, p_col_i=-1, p_col_j=-1;
                int p_col_lo=-1, p_col_hi=-1, p_col_var=-1;
                for (size_t i = 0; i < phdr.size(); ++i) {
                    if (phdr[i] == "output") p_col_out = i;
                    else if (phdr[i] == "input_i") p_col_i = i;
                    else if (phdr[i] == "input_j") p_col_j = i;
                    else if (phdr[i] == "S_ij_lo") p_col_lo = i;
                    else if (phdr[i] == "S_ij_hi") p_col_hi = i;
                    else if (phdr[i] == "output_variance") p_col_var = i;
                }
                if (p_col_lo >= 0 && p_col_hi >= 0 && p_col_var >= 0) {
                    std::vector<double> sij_widths, sij_N;
                    while (std::fgets(pline, sizeof(pline), fpp)) {
                        std::vector<std::string> tok;
                        std::string s = pline; size_t p = 0;
                        while (p < s.size()) {
                            size_t e = s.find(',', p);
                            if (e == std::string::npos) e = s.size();
                            std::string t = s.substr(p, e - p);
                            while (!t.empty() && (t.back()=='\n' || t.back()=='\r')) t.pop_back();
                            tok.push_back(t);
                            p = e + 1;
                        }
                        if ((int)tok.size() <= p_col_hi) continue;
                        const std::string& out_name = tok[p_col_out];
                        double v = std::atof(tok[p_col_var].c_str());
                        if (!(v > 0)) continue;
                        if (targets.skip.count(out_name)) continue;
                        double lo = std::atof(tok[p_col_lo].c_str());
                        double hi = std::atof(tok[p_col_hi].c_str());
                        double w  = 0.5 * (hi - lo);
                        // Skip truly degenerate rows (estimator unstable).
                        // 10.0 is a safety net; recommendation uses p90.
                        if (w > 10.0) continue;
                        out.n_pair_rows++;
                        sij_widths.push_back(w);
                        double sij_tgt = targets.Sij_for(out_name);
                        if (sij_tgt > 0) {
                            double Nn = static_cast<double>(pilot_N) * (w/sij_tgt) * (w/sij_tgt);
                            sij_N.push_back(Nn);
                            if (Nn > out.worst_Sij_N) {
                                out.worst_Sij_N = Nn;
                                out.driver_Sij_out = out_name;
                                out.driver_Sij_inp = tok[p_col_i] + " × " + tok[p_col_j];
                                out.driver_Sij_tgt = sij_tgt;
                            }
                        }
                        if (w > out.worst_Sij_w) {
                            out.worst_Sij_w   = w;
                            out.worst_Sij_out = out_name;
                            out.worst_Sij_inp = tok[p_col_i] + " × " + tok[p_col_j];
                        }
                    }
                    out.p90_Sij_w    = percentile_inplace(sij_widths, 0.90);
                    out.p90_Sij_N    = percentile_inplace(sij_N, 0.90);
                    out.median_Sij_N = percentile_inplace(sij_N, 0.50);
                }
            }
            std::fclose(fpp);
        }
    }

    // Triplet CSV scanning -- only if third_order was enabled.
    if (need_3rd) {
        std::string trip_path = tmp_csv;
        size_t dot = trip_path.find_last_of('.');
        if (dot != std::string::npos)
            trip_path = trip_path.substr(0, dot) + "_triplets" + trip_path.substr(dot);
        else
            trip_path += "_triplets";
        std::FILE* fpt = std::fopen(trip_path.c_str(), "r");
        if (fpt) {
            char tline[2048];
            if (std::fgets(tline, sizeof(tline), fpt)) {
                std::vector<std::string> thdr;
                {
                    std::string s = tline; size_t p = 0;
                    while (p < s.size()) {
                        size_t e = s.find(',', p);
                        if (e == std::string::npos) e = s.size();
                        std::string t = s.substr(p, e - p);
                        while (!t.empty() && (t.back()=='\n' || t.back()=='\r')) t.pop_back();
                        thdr.push_back(t);
                        p = e + 1;
                    }
                }
                int t_col_out=-1, t_col_i=-1, t_col_j=-1, t_col_k=-1;
                int t_col_lo=-1, t_col_hi=-1, t_col_var=-1;
                for (size_t i = 0; i < thdr.size(); ++i) {
                    if (thdr[i] == "output") t_col_out = i;
                    else if (thdr[i] == "input_i") t_col_i = i;
                    else if (thdr[i] == "input_j") t_col_j = i;
                    else if (thdr[i] == "input_k") t_col_k = i;
                    else if (thdr[i] == "S_ijk_lo") t_col_lo = i;
                    else if (thdr[i] == "S_ijk_hi") t_col_hi = i;
                    else if (thdr[i] == "output_variance") t_col_var = i;
                }
                if (t_col_lo >= 0 && t_col_hi >= 0 && t_col_var >= 0) {
                    std::vector<double> sijk_widths, sijk_N;
                    while (std::fgets(tline, sizeof(tline), fpt)) {
                        std::vector<std::string> tok;
                        std::string s = tline; size_t p = 0;
                        while (p < s.size()) {
                            size_t e = s.find(',', p);
                            if (e == std::string::npos) e = s.size();
                            std::string t = s.substr(p, e - p);
                            while (!t.empty() && (t.back()=='\n' || t.back()=='\r')) t.pop_back();
                            tok.push_back(t);
                            p = e + 1;
                        }
                        if ((int)tok.size() <= t_col_hi) continue;
                        const std::string& out_name = tok[t_col_out];
                        double v = std::atof(tok[t_col_var].c_str());
                        if (!(v > 0)) continue;
                        if (targets.skip.count(out_name)) continue;
                        double lo = std::atof(tok[t_col_lo].c_str());
                        double hi = std::atof(tok[t_col_hi].c_str());
                        double w  = 0.5 * (hi - lo);
                        // Triplets have wider noise envelope.  20.0 is a
                        // safety net for truly degenerate rows;
                        // recommendation uses p90.
                        if (w > 20.0) continue;
                        out.n_trip_rows++;
                        sijk_widths.push_back(w);
                        double sijk_tgt = targets.Sijk_for(out_name);
                        if (sijk_tgt > 0) {
                            double Nn = static_cast<double>(pilot_N) * (w/sijk_tgt) * (w/sijk_tgt);
                            sijk_N.push_back(Nn);
                            if (Nn > out.worst_Sijk_N) {
                                out.worst_Sijk_N = Nn;
                                out.driver_Sijk_out = out_name;
                                out.driver_Sijk_inp = tok[t_col_i] + " × "
                                                    + tok[t_col_j] + " × "
                                                    + tok[t_col_k];
                                out.driver_Sijk_tgt = sijk_tgt;
                            }
                        }
                        if (w > out.worst_Sijk_w) {
                            out.worst_Sijk_w   = w;
                            out.worst_Sijk_out = out_name;
                            out.worst_Sijk_inp = tok[t_col_i] + " × "
                                                + tok[t_col_j] + " × "
                                                + tok[t_col_k];
                        }
                    }
                    out.p90_Sijk_w    = percentile_inplace(sijk_widths, 0.90);
                    out.p90_Sijk_N    = percentile_inplace(sijk_N, 0.90);
                    out.median_Sijk_N = percentile_inplace(sijk_N, 0.50);
                }
            }
            std::fclose(fpt);
        }
    }

    // Cleanup
    std::remove(tmp_cfg);
    std::remove(tmp_csv);
    std::remove(tmp_txt);
    {
        // Group/pair CSVs may have been emitted; remove if present.
        std::string base = tmp_csv;
        size_t dot = base.find_last_of('.');
        std::string stem = (dot != std::string::npos) ? base.substr(0, dot) : base;
        std::string ext  = (dot != std::string::npos) ? base.substr(dot) : "";
        std::remove((stem + "_groups" + ext).c_str());
        std::remove((stem + "_pairs"  + ext).c_str());
        std::remove((stem + "_triplets" + ext).c_str());
        std::remove((stem + "_group_pairs" + ext).c_str());
    }
    return 0;
}

int run_sobol(const std::string& config_path,
              const std::string& indices_csv,
              const std::string& summary_txt)
{
    // ---- Parse config ----
    json::Value cfg;
    try {
        cfg = json::parse_file(config_path);
    } catch (const json::ParseError& e) {
        std::fprintf(stderr, "ERROR parsing %s at line %d col %d: %s\n",
                     config_path.c_str(), e.line, e.col, e.what());
        return 1;
    }

    const json::Value& mc = cfg["monte_carlo"];
    if (!mc.isObject()) {
        std::fprintf(stderr, "ERROR: config has no 'monte_carlo' section\n");
        return 1;
    }

    // Read Sobol-specific override block, falling back to monte_carlo values
    const json::Value& so = cfg["sobol"];
    int          N_base = so["n_base"].asInt(mc["n_runs"].asInt(512));
    unsigned long seed   = static_cast<unsigned long>(
                              so["seed"].asNumber(mc["seed"].asNumber(42)));
    int          n_boot = so["n_bootstrap"].asInt(0);  // 0 disables bootstrap CI
    // Bootstrap confidence level (default 0.95).  Must be in (0, 1).
    // 0.90 -> 5%/95% percentiles; 0.99 -> 0.5%/99.5%.
    double       ci_level = so["bootstrap_ci_level"].asNumber(0.95);
    if (!(ci_level > 0.0 && ci_level < 1.0)) {
        std::fprintf(stderr, "WARN: sobol.bootstrap_ci_level must be in (0,1); "
                     "got %.4g; using 0.95.\n", ci_level);
        ci_level = 0.95;
    }
    // Default 0 means auto-detect (capped at 8); falls back to mc.n_workers
    // if sobol.n_workers is absent.  Explicit user values >= 1 are honored.
    int          raw_workers = so["n_workers"].asInt(mc["n_workers"].asInt(0));
    bool         do_second_order = so["second_order"].asBool(false);
    bool         do_third_order  = so["third_order"].asBool(false);
    bool         do_group_pairs  = so["group_pairs"].asBool(false);
    if (do_third_order && !do_second_order) {
        std::fprintf(stderr, "WARN: sobol.third_order requires second_order; "
                     "enabling second_order automatically.\n");
        do_second_order = true;
    }
    ResolvedWorkers rw = resolve_n_workers(raw_workers);
    int          n_workers = rw.n;
    if (N_base <= 0) {
        std::fprintf(stderr, "ERROR: sobol.n_base must be > 0\n"); return 1;
    }
    if (n_boot < 0) n_boot = 0;

    // ---- Embedded auto-scale mode ----
    // If sobol.target_ci_width is set (either scalar or object) AND
    // we're not already inside a pilot run, run a pilot first and
    // override N_base with the recommendation.  The pilot is gated
    // by an env-var sentinel to avoid infinite recursion: when
    // run_pilot_and_measure spawns run_sobol on the temp config, it
    // sets SOBOL_IN_PILOT=1 so the recursive call skips this block.
    if (!std::getenv("SOBOL_IN_PILOT")) {
        TargetWidths targets = read_targets(so);
        bool has_target = targets.any();
        if (has_target) {
            int pilot_N = 64;
            if (const char* env = std::getenv("SOBOL_PILOT_N")) {
                int v = std::atoi(env);
                if (v > 0) pilot_N = v;
            }
            std::printf("================================================================\n");
            std::printf("    Embedded auto-scale: running pilot at N=%d\n", pilot_N);
            if (targets.has_ST())    std::printf("    Target ST    half-width: %.4g\n", targets.ST);
            if (targets.has_S1())    std::printf("    Target S1    half-width: %.4g\n", targets.S1);
            if (targets.has_S_ij())  std::printf("    Target S_ij  half-width: %.4g  (pilot will enable second_order)\n", targets.S_ij);
            if (targets.has_S_ijk()) std::printf("    Target S_ijk half-width: %.4g  (pilot will enable third_order)\n", targets.S_ijk);
            std::printf("================================================================\n\n");
            setenv("SOBOL_IN_PILOT", "1", 1);
            PilotResult pr;
            int rc = run_pilot_and_measure(config_path, pilot_N, targets, pr);
            unsetenv("SOBOL_IN_PILOT");
            if (rc != 0) {
                std::fprintf(stderr, "ERROR: embedded pilot failed; "
                             "falling back to configured N=%d\n", N_base);
            } else {
                int N_rec = 1;
                auto bump_from_N = [&](double Nf) {
                    if (Nf <= 0) return;
                    int N = 1; while (N < Nf) N *= 2;
                    if (N < pilot_N) N = pilot_N;
                    if (N > N_rec) N_rec = N;
                };
                // Default to p90 (more robust); override with
                //   SOBOL_RECOMMEND_FROM=worst   conservative
                //   SOBOL_RECOMMEND_FROM=median  forgiving
                //   SOBOL_RECOMMEND_FROM=auto    adaptive (per-estimator
                //                                gap-based selection)
                enum class Strat { P90, WORST, MEDIAN, AUTO };
                Strat s = Strat::P90;
                if (const char* env = std::getenv("SOBOL_RECOMMEND_FROM")) {
                    std::string sv = env;
                    if      (sv == "worst")  s = Strat::WORST;
                    else if (sv == "median") s = Strat::MEDIAN;
                    else if (sv == "auto")   s = Strat::AUTO;
                }
                // pick_adaptive applies the same gap-based rule as
                // suggest_n: ≤3x → worst, ≤10x → p90, >10x → median.
                auto pick_adaptive = [](double w, double p, double m) {
                    if (w <= 0 && p <= 0) return 0.0;
                    double gap = (p > 0) ? (w / p) : 1.0;
                    if (gap <= 3.0)  return w;
                    if (gap <= 10.0) return p;
                    return m;
                };
                auto pick = [&](double w, double p, double m) {
                    return s == Strat::WORST  ? w
                         : s == Strat::MEDIAN ? m
                         : s == Strat::AUTO   ? pick_adaptive(w, p, m)
                         :                       p;
                };
                bump_from_N(pick(pr.worst_ST_N,   pr.p90_ST_N,   pr.median_ST_N));
                bump_from_N(pick(pr.worst_S1_N,   pr.p90_S1_N,   pr.median_S1_N));
                bump_from_N(pick(pr.worst_Sij_N,  pr.p90_Sij_N,  pr.median_Sij_N));
                bump_from_N(pick(pr.worst_Sijk_N, pr.p90_Sijk_N, pr.median_Sijk_N));
                const char* sname = s == Strat::WORST  ? "worst"
                                  : s == Strat::MEDIAN ? "median"
                                  : s == Strat::AUTO   ? "auto"   : "p90";
                std::printf("\nPilot recommends N=%d (was configured: N=%d) "
                            "[strategy: %s]\n", N_rec, N_base, sname);
                if (N_rec > N_base) {
                    std::printf("Embedded auto-scale: bumping N_base %d -> %d\n",
                                N_base, N_rec);
                    N_base = N_rec;
                } else {
                    std::printf("Embedded auto-scale: configured N (%d) already meets "
                                "target, keeping it.\n", N_base);
                }
                std::printf("================================================================\n\n");
            }
        }
    }

    // ---- Variables and distributions ----
    const json::Value& vars = mc["variables"];
    if (!vars.isArray() || vars.size() == 0) {
        std::fprintf(stderr, "ERROR: monte_carlo.variables must be non-empty\n");
        return 1;
    }
    struct VarSpec {
        std::string path;
        std::unique_ptr<Distribution> dist;
    };
    std::vector<VarSpec> spec_list;
    spec_list.reserve(vars.size());
    for (size_t i = 0; i < vars.size(); ++i) {
        VarSpec vs;
        vs.path = vars[i]["path"].asString();
        try {
            vs.dist = make_distribution(vars[i]);
        } catch (const std::exception& e) {
            std::fprintf(stderr, "ERROR variable '%s': %s\n", vs.path.c_str(), e.what());
            return 1;
        }
        spec_list.push_back(std::move(vs));
    }
    const int k = static_cast<int>(spec_list.size());

    // ---- Outputs ----
    const json::Value& outs = mc["outputs"];
    std::vector<std::string> output_names;
    if (outs.isArray() && outs.size() > 0) {
        for (size_t i = 0; i < outs.size(); ++i) {
            output_names.push_back(outs[i].asString());
        }
    } else {
        output_names = { "alt", "dvbi", "fpa", "theta", "psi", "t_end" };
    }
    const int n_out = static_cast<int>(output_names.size());

    // ---- Parse optional grouped-Sobol spec ----
    // Format:
    //   "sobol": {
    //     "groups": [
    //       { "name": "propulsion", "paths": ["propulsion.fmass0", ...] },
    //       { "name": "guidance",   "paths": ["launch.thtbdx0"] }
    //     ]
    //   }
    // For each group G, adds 2*N sims (XAB_G and XAB_~G swap matrices)
    // and computes S_G and ST_G.  Cost is O(g * N) where g = #groups.
    struct Group {
        std::string name;
        std::vector<int> indices;  // positions in spec_list
    };
    std::vector<Group> groups;
    const auto& groups_v = so["groups"];
    if (groups_v.isArray()) {
        // Build path -> spec_list index map
        std::unordered_map<std::string, int> path_to_idx;
        for (int i = 0; i < k; ++i) path_to_idx[spec_list[i].path] = i;
        for (size_t gi = 0; gi < groups_v.size(); ++gi) {
            Group g;
            g.name = groups_v[gi]["name"].asString("group_" + std::to_string(gi));
            const auto& paths_v = groups_v[gi]["paths"];
            if (!paths_v.isArray()) continue;
            for (size_t pi = 0; pi < paths_v.size(); ++pi) {
                std::string p = paths_v[pi].asString();
                auto it = path_to_idx.find(p);
                if (it == path_to_idx.end()) {
                    std::fprintf(stderr, "WARN: group '%s' references unknown path '%s' (skipping)\n",
                                 g.name.c_str(), p.c_str());
                    continue;
                }
                g.indices.push_back(it->second);
            }
            if (!g.indices.empty()) groups.push_back(g);
        }
    }
    const int n_groups = static_cast<int>(groups.size());

    // ---- Banner ----
    // Total sims = base (N*(k+2)) + groups (N*n_groups) + 2nd-order (N*k*(k-1)/2)
    //            + group pairs (N*ng*(ng-1)/2 if group_pairs enabled)
    int n_pairs = do_second_order ? (k * (k - 1) / 2) : 0;
    int n_triplets = do_third_order ? (k * (k - 1) * (k - 2) / 6) : 0;
    int n_group_pairs = (do_group_pairs && n_groups >= 2)
                       ? (n_groups * (n_groups - 1) / 2) : 0;
    int total_runs = N_base * (k + 2) + N_base * n_groups + N_base * n_pairs
                     + N_base * n_triplets + N_base * n_group_pairs;
    std::printf("================================================================\n");
    std::printf("    rocket_6dof  --  Sobol sensitivity analysis\n");
    std::printf("================================================================\n");
    std::printf("Config:        %s\n", config_path.c_str());
    std::printf("Build:         MAX_STAGES=%d\n", Propulsion::MAX_STAGES);
    std::printf("Seed:          %lu\n", seed);
    std::printf("Base N:        %d\n", N_base);
    std::printf("Variables:     %d\n", k);
    std::printf("Outputs:       %d\n", n_out);
    if (n_groups > 0) {
        std::printf("Groups:        %d (+ N sims per group)\n", n_groups);
        for (const auto& g : groups) {
            std::printf("  - %s: %zu vars\n", g.name.c_str(), g.indices.size());
        }
    }
    if (do_second_order) {
        std::printf("Second-order:  %d pairs (+ N sims per pair)\n", n_pairs);
    }
    if (do_third_order) {
        std::printf("Third-order:   %d triplets (+ N sims per triplet)\n", n_triplets);
    }
    if (n_group_pairs > 0) {
        std::printf("Group pairs:   %d (+ N sims per group pair)\n", n_group_pairs);
    }
    std::printf("Total sims:    %d\n", total_runs);
    std::printf("Workers:       %d  [%s]%s\n", n_workers, rw.source.c_str(),
                n_workers > 1 ? "  (subprocess parallel)" : "  (sequential)");
    if (n_boot > 0) {
        std::printf("Bootstrap:     %d resamples (95%% CI)\n", n_boot);
    } else {
        std::printf("Bootstrap:     disabled (set sobol.n_bootstrap > 0)\n");
    }
    std::printf("================================================================\n\n");

    // ---- Generate A, B uniform matrices [0,1)^(N x k) ----
    // and physical samples X_A, X_B by mapping through each variable's quantile().
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U(0.0, 1.0);

    auto draw_matrix = [&](int N) {
        std::vector<std::vector<double>> M(k, std::vector<double>(N));
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < N; ++j) M[i][j] = U(rng);
        }
        return M;
    };
    auto map_to_physical = [&](const std::vector<std::vector<double>>& M) {
        std::vector<std::vector<double>> X(k, std::vector<double>(M[0].size()));
        for (int i = 0; i < k; ++i) {
            for (size_t j = 0; j < M[i].size(); ++j) {
                X[i][j] = spec_list[i].dist->quantile(M[i][j]);
            }
        }
        return X;
    };

    auto A = draw_matrix(N_base);
    auto B = draw_matrix(N_base);
    auto XA = map_to_physical(A);
    auto XB = map_to_physical(B);

    // ---- Per-run helper ----
    auto run_at = [&](const std::vector<std::vector<double>>& X, int j,
                      std::vector<double>& Y_out) {
        json::Value run_cfg = cfg;
        for (int i = 0; i < k; ++i) {
            run_cfg.set_path(spec_list[i].path, X[i][j]);
        }
        SingleRunResult res = run_single_mission(run_cfg, output_names);
        Y_out.assign(res.outputs.begin(), res.outputs.end());
    };

    // ---- Storage for all outputs (one Y vector of length N per output) ----
    std::vector<std::vector<double>> YA(n_out, std::vector<double>(N_base));
    std::vector<std::vector<double>> YB(n_out, std::vector<double>(N_base));
    // YAB[o][i][j]: output o, input i, sample j
    std::vector<std::vector<std::vector<double>>> YAB(
        n_out, std::vector<std::vector<double>>(k, std::vector<double>(N_base)));
    // YAB_G[o][g][j]: output o, group g, sample j (empty if no groups)
    std::vector<std::vector<std::vector<double>>> YAB_G(
        n_out, std::vector<std::vector<double>>(n_groups, std::vector<double>(N_base)));
    // YBA[o][p][j]: output o, pair-index p, sample j (empty unless do_second_order)
    // Pair p enumerates (i,j) with i<j in lex order: (0,1),(0,2),...,(0,k-1),(1,2),...
    std::vector<std::pair<int,int>> pair_ij;
    pair_ij.reserve(n_pairs);
    for (int i = 0; i < k; ++i) {
        for (int j = i + 1; j < k; ++j) pair_ij.emplace_back(i, j);
    }
    std::vector<std::vector<std::vector<double>>> YBA(
        n_out, std::vector<std::vector<double>>(n_pairs, std::vector<double>(N_base)));

    // Group pairs: enumerate (gi, gj) with gi < gj
    std::vector<std::pair<int,int>> group_pair_ij;
    group_pair_ij.reserve(n_group_pairs);
    for (int gi = 0; gi < n_groups; ++gi) {
        for (int gj = gi + 1; gj < n_groups; ++gj) {
            group_pair_ij.emplace_back(gi, gj);
        }
    }
    std::vector<std::vector<std::vector<double>>> YBA_G(
        n_out, std::vector<std::vector<double>>(n_group_pairs, std::vector<double>(N_base)));

    // Triplets: enumerate (i,j,k_) with i<j<k_, lex order.
    // YBA_3[o][t][n]: output o, triplet-index t, sample n (empty unless third_order).
    struct Trip { int i, j, k_; };
    std::vector<Trip> trip_ijk;
    trip_ijk.reserve(n_triplets);
    for (int i = 0; i < k; ++i) {
        for (int j = i + 1; j < k; ++j) {
            for (int k_ = j + 1; k_ < k; ++k_) {
                trip_ijk.push_back({i, j, k_});
            }
        }
    }
    std::vector<std::vector<std::vector<double>>> YBA_3(
        n_out, std::vector<std::vector<double>>(n_triplets, std::vector<double>(N_base)));

    int sim_done = 0;
    auto t0 = std::chrono::steady_clock::now();
    auto progress = [&](const char* phase, int idx, int total) {
        sim_done++;
        if (sim_done % std::max(1, total_runs / 40) == 0 || sim_done == total_runs) {
            auto t1 = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(t1 - t0).count();
            double per     = elapsed / sim_done;
            double eta     = per * (total_runs - sim_done);
            std::printf("  %-7s %4d/%-4d   sim %5d/%-5d   %.3f s/run   ETA %.0f s\n",
                        phase, idx, total, sim_done, total_runs, per, eta);
        }
    };

    // Build XAB_i and XAB_G matrices.  We materialize them all up
    // front because (a) the parallel path needs a flat row-major sample
    // matrix covering every phase, and (b) it lets us share one
    // codepath between sequential and parallel execution.
    std::vector<std::vector<std::vector<double>>> XAB(k);   // [i][col][row]
    for (int i = 0; i < k; ++i) {
        XAB[i] = XA;
        XAB[i][i] = XB[i];  // replace column i with B's column
    }
    std::vector<std::vector<std::vector<double>>> XABG(n_groups);
    for (int gi = 0; gi < n_groups; ++gi) {
        XABG[gi] = XA;
        for (int idx : groups[gi].indices) XABG[gi][idx] = XB[idx];
    }
    // XBA[p]: XB with columns i and j (from pair_ij[p]) replaced with A's columns.
    // For pair (i,j) the estimator V(E[Y|Xi,Xj]) = (1/N) Σ Y_A · (Y_BA_ij - Y_B).
    std::vector<std::vector<std::vector<double>>> XBA(n_pairs);
    for (int p = 0; p < n_pairs; ++p) {
        XBA[p] = XB;
        XBA[p][pair_ij[p].first]  = XA[pair_ij[p].first];
        XBA[p][pair_ij[p].second] = XA[pair_ij[p].second];
    }
    // XBA_G[p]: XB with ALL columns from group_pair_ij[p].first AND
    //           group_pair_ij[p].second replaced with A's columns.
    // Estimator: V(E[Y|X_Ga, X_Gb]) ≈ (1/N) Σ Y_A · (Y_BA_Ga_Gb - Y_B)
    std::vector<std::vector<std::vector<double>>> XBA_G(n_group_pairs);
    for (int p = 0; p < n_group_pairs; ++p) {
        XBA_G[p] = XB;
        int ga = group_pair_ij[p].first;
        int gb = group_pair_ij[p].second;
        for (int idx : groups[ga].indices) XBA_G[p][idx] = XA[idx];
        for (int idx : groups[gb].indices) XBA_G[p][idx] = XA[idx];
    }

    // XBA_3[t]: XB with columns i, j, k_ from XA (for triplet t).
    // Estimator: V(E[Y|X_i, X_j, X_k]) ≈ (1/N) Σ Y_A · (Y_BA_ijk - Y_B)
    std::vector<std::vector<std::vector<double>>> XBA_3(n_triplets);
    for (int t = 0; t < n_triplets; ++t) {
        XBA_3[t] = XB;
        XBA_3[t][trip_ijk[t].i]  = XA[trip_ijk[t].i];
        XBA_3[t][trip_ijk[t].j]  = XA[trip_ijk[t].j];
        XBA_3[t][trip_ijk[t].k_] = XA[trip_ijk[t].k_];
    }

    if (n_workers >= 2) {
        // ---- Parallel path ----
        // Concatenate all phase matrices into one flat row-major sample
        // matrix.  Total rows = total_runs.
        // Row layout:
        //   rows [0, N):                                          XA
        //   rows [N, 2N):                                         XB
        //   rows [2N + i*N, 2N + (i+1)*N) for i in [0, k):        XAB[i]
        //   rows [2N + k*N + g*N, ...)    for g in [0, n_groups): XABG[g]
        //   rows [2N + k*N + n_groups*N + p*N, ...) for p in [0, n_pairs): XBA[p]
        std::vector<std::string> input_paths(k);
        for (int i = 0; i < k; ++i) input_paths[i] = spec_list[i].path;

        std::vector<std::vector<double>> sample_matrix(total_runs,
            std::vector<double>(k));

        auto fill_rows = [&](int row_base,
                             const std::vector<std::vector<double>>& X) {
            for (int j = 0; j < N_base; ++j) {
                for (int i = 0; i < k; ++i) {
                    sample_matrix[row_base + j][i] = X[i][j];
                }
            }
        };

        fill_rows(0,         XA);
        fill_rows(N_base,    XB);
        for (int i = 0; i < k; ++i) {
            fill_rows(2*N_base + i*N_base, XAB[i]);
        }
        const int group_base = 2*N_base + k*N_base;
        for (int gi = 0; gi < n_groups; ++gi) {
            fill_rows(group_base + gi*N_base, XABG[gi]);
        }
        const int pair_base = group_base + n_groups*N_base;
        for (int p = 0; p < n_pairs; ++p) {
            fill_rows(pair_base + p*N_base, XBA[p]);
        }
        const int gpair_base = pair_base + n_pairs*N_base;
        for (int p = 0; p < n_group_pairs; ++p) {
            fill_rows(gpair_base + p*N_base, XBA_G[p]);
        }
        const int trip_base = gpair_base + n_group_pairs*N_base;
        for (int t = 0; t < n_triplets; ++t) {
            fill_rows(trip_base + t*N_base, XBA_3[t]);
        }

        ParallelResult pr = run_chunked_parallel(
            config_path, input_paths, sample_matrix, output_names,
            n_workers, "./mission_worker", "/tmp", "Sobol");
        if (!pr.ok) {
            std::fprintf(stderr, "ERROR: parallel Sobol dispatch failed\n");
            return 1;
        }
        sim_done = total_runs;

        // Unpack results back into YA / YB / YAB / YAB_G / YBA
        for (int j = 0; j < N_base; ++j) {
            for (int o = 0; o < n_out; ++o) YA[o][j] = pr.row_outputs[j][o];
        }
        for (int j = 0; j < N_base; ++j) {
            for (int o = 0; o < n_out; ++o) YB[o][j] = pr.row_outputs[N_base + j][o];
        }
        for (int i = 0; i < k; ++i) {
            int rb = 2*N_base + i*N_base;
            for (int j = 0; j < N_base; ++j) {
                for (int o = 0; o < n_out; ++o)
                    YAB[o][i][j] = pr.row_outputs[rb + j][o];
            }
        }
        for (int gi = 0; gi < n_groups; ++gi) {
            int rb = group_base + gi*N_base;
            for (int j = 0; j < N_base; ++j) {
                for (int o = 0; o < n_out; ++o)
                    YAB_G[o][gi][j] = pr.row_outputs[rb + j][o];
            }
        }
        for (int p = 0; p < n_pairs; ++p) {
            int rb = pair_base + p*N_base;
            for (int j = 0; j < N_base; ++j) {
                for (int o = 0; o < n_out; ++o)
                    YBA[o][p][j] = pr.row_outputs[rb + j][o];
            }
        }
        for (int p = 0; p < n_group_pairs; ++p) {
            int rb = gpair_base + p*N_base;
            for (int j = 0; j < N_base; ++j) {
                for (int o = 0; o < n_out; ++o)
                    YBA_G[o][p][j] = pr.row_outputs[rb + j][o];
            }
        }
        for (int t = 0; t < n_triplets; ++t) {
            int rb = trip_base + t*N_base;
            for (int j = 0; j < N_base; ++j) {
                for (int o = 0; o < n_out; ++o)
                    YBA_3[o][t][j] = pr.row_outputs[rb + j][o];
            }
        }
    } else {
        // ---- Sequential path ----
        std::vector<double> Y;

        // Phase 1: A
        for (int j = 0; j < N_base; ++j) {
            run_at(XA, j, Y);
            for (int o = 0; o < n_out; ++o) YA[o][j] = Y[o];
            progress("A", j+1, N_base);
        }
        // Phase 2: B
        for (int j = 0; j < N_base; ++j) {
            run_at(XB, j, Y);
            for (int o = 0; o < n_out; ++o) YB[o][j] = Y[o];
            progress("B", j+1, N_base);
        }
        // Phase 3: AB_i for each i
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < N_base; ++j) {
                run_at(XAB[i], j, Y);
                for (int o = 0; o < n_out; ++o) YAB[o][i][j] = Y[o];
                progress(("AB_" + std::to_string(i)).c_str(), j+1, N_base);
            }
        }
        // Phase 4: AB_G for each group
        for (int gi = 0; gi < n_groups; ++gi) {
            for (int j = 0; j < N_base; ++j) {
                run_at(XABG[gi], j, Y);
                for (int o = 0; o < n_out; ++o) YAB_G[o][gi][j] = Y[o];
                progress(("AB_G_" + groups[gi].name).c_str(), j+1, N_base);
            }
        }
        // Phase 5: BA_ij for each pair (i<j)  -- only if second_order enabled
        for (int p = 0; p < n_pairs; ++p) {
            for (int j = 0; j < N_base; ++j) {
                run_at(XBA[p], j, Y);
                for (int o = 0; o < n_out; ++o) YBA[o][p][j] = Y[o];
                char tag[16];
                std::snprintf(tag, sizeof(tag), "BA_%d_%d",
                              pair_ij[p].first, pair_ij[p].second);
                progress(tag, j+1, N_base);
            }
        }
        // Phase 6: BA_Ga_Gb for each group pair  -- only if group_pairs enabled
        for (int p = 0; p < n_group_pairs; ++p) {
            for (int j = 0; j < N_base; ++j) {
                run_at(XBA_G[p], j, Y);
                for (int o = 0; o < n_out; ++o) YBA_G[o][p][j] = Y[o];
                std::string tag = "BA_G_" + groups[group_pair_ij[p].first].name
                                + "_" + groups[group_pair_ij[p].second].name;
                progress(tag.c_str(), j+1, N_base);
            }
        }
        // Phase 7: BA_ijk for each triplet  -- only if third_order enabled
        for (int t = 0; t < n_triplets; ++t) {
            for (int j = 0; j < N_base; ++j) {
                run_at(XBA_3[t], j, Y);
                for (int o = 0; o < n_out; ++o) YBA_3[o][t][j] = Y[o];
                char tag[24];
                std::snprintf(tag, sizeof(tag), "BA_%d_%d_%d",
                              trip_ijk[t].i, trip_ijk[t].j, trip_ijk[t].k_);
                progress(tag, j+1, N_base);
            }
        }
    }

    // ---- Compute indices ----
    std::vector<Indices> per_out(n_out);
    for (int o = 0; o < n_out; ++o) {
        per_out[o] = compute_indices(YA[o], YB[o], YAB[o]);
    }

    // Group indices per output: both S_G and ST_G use YAB_G
    struct GroupIndices {
        std::vector<double> S1;  // length n_groups
        std::vector<double> ST;
    };
    std::vector<GroupIndices> per_out_group(n_out);
    for (int o = 0; o < n_out; ++o) {
        per_out_group[o].S1.assign(n_groups, 0.0);
        per_out_group[o].ST.assign(n_groups, 0.0);
        double varA = per_out[o].var;
        if (!(varA > 0) || !std::isfinite(varA)) continue;
        const auto& ya = YA[o]; const auto& yb = YB[o];
        for (int gi = 0; gi < n_groups; ++gi) {
            const auto& yabg = YAB_G[o][gi];
            double s_num = 0, t_num = 0;
            int n_pairs = 0;
            for (int j = 0; j < N_base; ++j) {
                if (!std::isfinite(ya[j]) || !std::isfinite(yb[j])
                    || !std::isfinite(yabg[j])) continue;
                s_num += yb[j] * (yabg[j] - ya[j]);
                double d = ya[j] - yabg[j];
                t_num += d * d;
                ++n_pairs;
            }
            if (n_pairs < 2) {
                per_out_group[o].S1[gi] = 0.0;
                per_out_group[o].ST[gi] = 0.0;
            } else {
                per_out_group[o].S1[gi] = (s_num / n_pairs) / varA;
                per_out_group[o].ST[gi] = (t_num / (2.0 * n_pairs)) / varA;
            }
        }
    }

    // Second-order indices per output, per pair (i,j).
    //   V(E[Y|X_i,X_j]) ≈ (1/N) Σ Y_A · (Y_BA_ij - Y_B)
    //   S_ij = (V_closed_ij - V_i - V_j) / V(Y)
    // where V_i = S1[i] * V (already in per_out[o].S1).
    struct PairIndices {
        std::vector<double> Sij;  // length n_pairs
    };
    std::vector<PairIndices> per_out_pairs(n_out);
    for (int o = 0; o < n_out; ++o) {
        per_out_pairs[o].Sij.assign(n_pairs, 0.0);
        double varA = per_out[o].var;
        if (!(varA > 0) || !std::isfinite(varA)) continue;
        const auto& ya = YA[o]; const auto& yb = YB[o];
        for (int p = 0; p < n_pairs; ++p) {
            int i = pair_ij[p].first;
            int j = pair_ij[p].second;
            const auto& yba = YBA[o][p];
            double c_num = 0;
            int n_pairs_p = 0;
            for (int n = 0; n < N_base; ++n) {
                if (!std::isfinite(ya[n]) || !std::isfinite(yb[n])
                    || !std::isfinite(yba[n])) continue;
                c_num += ya[n] * (yba[n] - yb[n]);
                ++n_pairs_p;
            }
            if (n_pairs_p < 2) {
                per_out_pairs[o].Sij[p] = 0.0;
            } else {
                double V_closed = c_num / n_pairs_p;  // ≈ V(E[Y|X_i,X_j])
                double V_i = per_out[o].S1[i] * varA;
                double V_j = per_out[o].S1[j] * varA;
                per_out_pairs[o].Sij[p] = (V_closed - V_i - V_j) / varA;
            }
        }
    }

    // Second-order group-pair indices S_{Ga,Gb}.
    //   V(E[Y|X_Ga, X_Gb]) ≈ (1/N) Σ Y_A · (Y_BA_Ga_Gb - Y_B)
    //   S_{Ga,Gb} = (V_closed - V_Ga - V_Gb) / V
    // where V_Ga = per_out_group.S1[Ga] * varA.
    struct GroupPairIndices {
        std::vector<double> SGij;  // length n_group_pairs
    };
    std::vector<GroupPairIndices> per_out_gpair(n_out);
    for (int o = 0; o < n_out; ++o) {
        per_out_gpair[o].SGij.assign(n_group_pairs, 0.0);
        double varA = per_out[o].var;
        if (!(varA > 0) || !std::isfinite(varA)) continue;
        const auto& ya = YA[o]; const auto& yb = YB[o];
        for (int p = 0; p < n_group_pairs; ++p) {
            int ga = group_pair_ij[p].first;
            int gb = group_pair_ij[p].second;
            const auto& yba = YBA_G[o][p];
            double c_num = 0;
            int n_pairs_p = 0;
            for (int n = 0; n < N_base; ++n) {
                if (!std::isfinite(ya[n]) || !std::isfinite(yb[n])
                    || !std::isfinite(yba[n])) continue;
                c_num += ya[n] * (yba[n] - yb[n]);
                ++n_pairs_p;
            }
            if (n_pairs_p < 2) {
                per_out_gpair[o].SGij[p] = 0.0;
            } else {
                double V_closed = c_num / n_pairs_p;
                double V_Ga = per_out_group[o].S1[ga] * varA;
                double V_Gb = per_out_group[o].S1[gb] * varA;
                per_out_gpair[o].SGij[p] = (V_closed - V_Ga - V_Gb) / varA;
            }
        }
    }

    // Third-order indices S_ijk per output, per triplet (i<j<k_).
    //   V(E[Y|X_i,X_j,X_k]) ≈ (1/N) Σ Y_A · (Y_BA_ijk - Y_B)
    //   S_ijk = (V_closed_ijk - V_i - V_j - V_k - V_ij - V_ik - V_jk) / V
    // where V_p = per_out.S1[p] * V and V_pq = per_out_pairs.Sij[p,q] * V
    // (after looking up the pair index for the corresponding (i,j), etc.)
    struct TripletIndices {
        std::vector<double> Sijk;  // length n_triplets
    };
    std::vector<TripletIndices> per_out_trip(n_out);
    // Build a quick lookup: pair_index_of(i,j) -> p
    auto pair_index_of = [&](int a, int b) -> int {
        if (a > b) std::swap(a, b);
        // Pairs enumerated as (0,1),(0,2),...,(0,k-1),(1,2),...,(k-2,k-1).
        // For (i,j) with i<j: index = i*k - i*(i+1)/2 + (j - i - 1).
        return a * k - a * (a + 1) / 2 + (b - a - 1);
    };
    for (int o = 0; o < n_out; ++o) {
        per_out_trip[o].Sijk.assign(n_triplets, 0.0);
        double varA = per_out[o].var;
        if (!(varA > 0) || !std::isfinite(varA)) continue;
        const auto& ya = YA[o]; const auto& yb = YB[o];
        for (int t = 0; t < n_triplets; ++t) {
            int i = trip_ijk[t].i;
            int j = trip_ijk[t].j;
            int k_ = trip_ijk[t].k_;
            const auto& yba = YBA_3[o][t];
            double c_num = 0;
            int n_pairs_t = 0;
            for (int n = 0; n < N_base; ++n) {
                if (!std::isfinite(ya[n]) || !std::isfinite(yb[n])
                    || !std::isfinite(yba[n])) continue;
                c_num += ya[n] * (yba[n] - yb[n]);
                ++n_pairs_t;
            }
            if (n_pairs_t < 2) {
                per_out_trip[o].Sijk[t] = 0.0;
            } else {
                double V_closed = c_num / n_pairs_t;
                double V_i = per_out[o].S1[i] * varA;
                double V_j = per_out[o].S1[j] * varA;
                double V_k = per_out[o].S1[k_] * varA;
                double V_ij = per_out_pairs[o].Sij[pair_index_of(i, j)] * varA;
                double V_ik = per_out_pairs[o].Sij[pair_index_of(i, k_)] * varA;
                double V_jk = per_out_pairs[o].Sij[pair_index_of(j, k_)] * varA;
                per_out_trip[o].Sijk[t] = (V_closed - V_i - V_j - V_k
                                           - V_ij - V_ik - V_jk) / varA;
            }
        }
    }

    // ---- Bootstrap CI (optional) ----
    std::vector<CI> per_out_ci(n_out);
    std::vector<CI_pairs> per_out_ci_pairs(n_out);
    std::vector<CI_triplets> per_out_ci_trip(n_out);
    std::vector<CI_gpairs> per_out_ci_gpair(n_out);
    std::vector<CI_share> per_out_ci_share(n_out);

    // Convert trip_ijk to std::array<int,3> form (the bootstrap API
    // doesn't see the locally-defined Trip struct).
    std::vector<std::array<int,3>> trip_ijk_arr(n_triplets);
    for (int t = 0; t < n_triplets; ++t) {
        trip_ijk_arr[t] = {trip_ijk[t].i, trip_ijk[t].j, trip_ijk[t].k_};
    }

    if (n_boot > 0) {
        std::printf("\nBootstrap: computing %.0f%% CIs from %d resamples per output...\n",
                    100.0 * ci_level, n_boot);
        std::mt19937 boot_rng(seed + 0xDEADBEEFu);
        std::mt19937 pair_boot_rng(seed + 0xBEEFDEADu);
        std::mt19937 trip_boot_rng(seed + 0xFEEDFACEu);
        std::mt19937 gpair_boot_rng(seed + 0xCAFEBABEu);
        std::mt19937 share_boot_rng(seed + 0xD15EA5Eu);  // share-CI stream
        auto tb0 = std::chrono::steady_clock::now();
        for (int o = 0; o < n_out; ++o) {
            per_out_ci[o] = bootstrap_indices(YA[o], YB[o], YAB[o],
                                              n_boot, boot_rng, ci_level);
            if (do_second_order && n_pairs > 0) {
                per_out_ci_pairs[o] = bootstrap_pair_indices(
                    YA[o], YB[o], YAB[o], YBA[o], pair_ij,
                    n_boot, pair_boot_rng, ci_level);
            } else {
                per_out_ci_pairs[o].Sij_lo.assign(n_pairs, 0);
                per_out_ci_pairs[o].Sij_hi.assign(n_pairs, 0);
            }
            if (do_third_order && n_triplets > 0) {
                per_out_ci_trip[o] = bootstrap_triplet_indices(
                    YA[o], YB[o], YAB[o], YBA[o], YBA_3[o],
                    pair_ij, trip_ijk_arr, n_boot, trip_boot_rng, ci_level);
            } else {
                per_out_ci_trip[o].Sijk_lo.assign(n_triplets, 0);
                per_out_ci_trip[o].Sijk_hi.assign(n_triplets, 0);
            }
            if (n_group_pairs > 0) {
                per_out_ci_gpair[o] = bootstrap_gpair_indices(
                    YA[o], YB[o], YAB_G[o], YBA_G[o], group_pair_ij,
                    n_boot, gpair_boot_rng, ci_level);
            } else {
                per_out_ci_gpair[o].SGij_lo.assign(n_group_pairs, 0);
                per_out_ci_gpair[o].SGij_hi.assign(n_group_pairs, 0);
            }
            // Joint clean_share[i] bootstrap.  Uses a separate index
            // stream from the per-order bootstraps so they remain
            // statistically independent of each other, while the three
            // orders WITHIN this bootstrap share an index set per
            // resample (which is the whole point).
            // We only need this when pair data is available; without it
            // share is just S1 and its CI matches the existing first-
            // order CI.
            static const std::vector<std::vector<double>> EMPTY_YBA{};
            static const std::vector<std::pair<int,int>>  EMPTY_PAIR{};
            static const std::vector<std::array<int,3>>   EMPTY_TRIP{};
            const auto& YBA_use   = (do_second_order && n_pairs > 0)
                                  ? YBA[o]    : EMPTY_YBA;
            const auto& YBA_3_use = (do_third_order && n_triplets > 0)
                                  ? YBA_3[o]  : EMPTY_YBA;
            const auto& pair_use  = (do_second_order && n_pairs > 0)
                                  ? pair_ij   : EMPTY_PAIR;
            const auto& trip_use  = (do_third_order && n_triplets > 0)
                                  ? trip_ijk_arr : EMPTY_TRIP;
            per_out_ci_share[o] = bootstrap_clean_share(
                YA[o], YB[o], YAB[o], YBA_use, YBA_3_use,
                pair_use, trip_use, n_boot, share_boot_rng, ci_level);
        }
        auto tb1 = std::chrono::steady_clock::now();
        double boot_s = std::chrono::duration<double>(tb1 - tb0).count();
        std::printf("Bootstrap complete: %.2f s\n\n", boot_s);
    } else {
        for (int o = 0; o < n_out; ++o) {
            per_out_ci[o].S1_lo.assign(k, 0); per_out_ci[o].S1_hi.assign(k, 0);
            per_out_ci[o].ST_lo.assign(k, 0); per_out_ci[o].ST_hi.assign(k, 0);
            per_out_ci_pairs[o].Sij_lo.assign(n_pairs, 0);
            per_out_ci_pairs[o].Sij_hi.assign(n_pairs, 0);
            per_out_ci_trip[o].Sijk_lo.assign(n_triplets, 0);
            per_out_ci_trip[o].Sijk_hi.assign(n_triplets, 0);
            per_out_ci_gpair[o].SGij_lo.assign(n_group_pairs, 0);
            per_out_ci_gpair[o].SGij_hi.assign(n_group_pairs, 0);
            per_out_ci_share[o].share_lo.assign(k, 0);
            per_out_ci_share[o].share_hi.assign(k, 0);
        }
    }

    // ---- Write CSV ----
    std::FILE* fi = std::fopen(indices_csv.c_str(), "w");
    if (!fi) {
        std::fprintf(stderr, "ERROR: cannot open %s\n", indices_csv.c_str());
        return 1;
    }
    // Helper: point estimate of clean_share[i] from the production
    // arrays.  Mirrors the math used inside bootstrap_clean_share.
    auto compute_share = [&](int o, int i) -> double {
        double s = per_out[o].S1[i];
        if (do_second_order && n_pairs > 0) {
            for (int p = 0; p < n_pairs; ++p) {
                if (pair_ij[p].first == i || pair_ij[p].second == i) {
                    s += 0.5 * per_out_pairs[o].Sij[p];
                }
            }
        }
        if (do_third_order && n_triplets > 0) {
            for (int t = 0; t < n_triplets; ++t) {
                int ti = trip_ijk[t].i, tj = trip_ijk[t].j, tk = trip_ijk[t].k_;
                if (ti == i || tj == i || tk == i) {
                    s += per_out_trip[o].Sijk[t] / 3.0;
                }
            }
        }
        return s;
    };
    if (n_boot > 0) {
        std::fprintf(fi, "output,input,S1,ST,S1_lo,S1_hi,ST_lo,ST_hi,"
                          "clean_share,clean_share_lo,clean_share_hi,"
                          "output_mean,output_variance\n");
        for (int o = 0; o < n_out; ++o) {
            for (int i = 0; i < k; ++i) {
                double share = compute_share(o, i);
                std::fprintf(fi, "%s,%s,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,"
                                 "%.9g,%.9g,%.9g,%.9g,%.9g\n",
                             output_names[o].c_str(), spec_list[i].path.c_str(),
                             per_out[o].S1[i], per_out[o].ST[i],
                             per_out_ci[o].S1_lo[i], per_out_ci[o].S1_hi[i],
                             per_out_ci[o].ST_lo[i], per_out_ci[o].ST_hi[i],
                             share,
                             per_out_ci_share[o].share_lo[i],
                             per_out_ci_share[o].share_hi[i],
                             per_out[o].mean, per_out[o].var);
            }
        }
    } else {
        // Original schema for backward compatibility with downstream readers
        std::fprintf(fi, "output,input,S1,ST,output_mean,output_variance\n");
        for (int o = 0; o < n_out; ++o) {
            for (int i = 0; i < k; ++i) {
                std::fprintf(fi, "%s,%s,%.9g,%.9g,%.9g,%.9g\n",
                             output_names[o].c_str(), spec_list[i].path.c_str(),
                             per_out[o].S1[i], per_out[o].ST[i],
                             per_out[o].mean, per_out[o].var);
            }
        }
    }
    std::fclose(fi);

    // ---- Write targets sidecar (optional) ----
    // If the config has sobol.target_ci_width, emit a small JSON
    // sidecar listing the global defaults and per-output overrides.
    // The viewer can embed this to render a "targets in effect" header.
    {
        TargetWidths sidecar_t = read_targets(so);
        if (sidecar_t.any() || !sidecar_t.skip.empty()) {
            std::string targets_path = indices_csv;
            size_t dot = targets_path.find_last_of('.');
            if (dot != std::string::npos) {
                targets_path = targets_path.substr(0, dot) + "_targets.json";
            } else {
                targets_path += "_targets.json";
            }
            std::FILE* ft = std::fopen(targets_path.c_str(), "w");
            if (ft) {
                std::fprintf(ft, "{\n");
                // Globals
                std::fprintf(ft, "  \"global\": {");
                bool first = true;
                auto kv_num = [&](const char* key, double v) {
                    if (v > 0) {
                        std::fprintf(ft, "%s\"%s\": %.6g", first ? "" : ", ", key, v);
                        first = false;
                    }
                };
                kv_num("ST",    sidecar_t.ST);
                kv_num("S1",    sidecar_t.S1);
                kv_num("S_ij",  sidecar_t.S_ij);
                kv_num("S_ijk", sidecar_t.S_ijk);
                std::fprintf(ft, "},\n");
                // Per-output
                std::fprintf(ft, "  \"per_output\": {\n");
                std::set<std::string> all_outs;
                for (const auto& kv : sidecar_t.per_out_ST)    all_outs.insert(kv.first);
                for (const auto& kv : sidecar_t.per_out_S1)    all_outs.insert(kv.first);
                for (const auto& kv : sidecar_t.per_out_S_ij)  all_outs.insert(kv.first);
                for (const auto& kv : sidecar_t.per_out_S_ijk) all_outs.insert(kv.first);
                for (const auto& s : sidecar_t.skip)           all_outs.insert(s);
                bool out_first = true;
                for (const auto& out_name : all_outs) {
                    if (!out_first) std::fprintf(ft, ",\n");
                    out_first = false;
                    std::fprintf(ft, "    \"%s\": {", out_name.c_str());
                    bool ko_first = true;
                    auto pkv = [&](const char* key, const std::map<std::string,double>& m) {
                        auto it = m.find(out_name);
                        if (it != m.end()) {
                            std::fprintf(ft, "%s\"%s\": %.6g",
                                         ko_first ? "" : ", ", key, it->second);
                            ko_first = false;
                        }
                    };
                    pkv("ST",    sidecar_t.per_out_ST);
                    pkv("S1",    sidecar_t.per_out_S1);
                    pkv("S_ij",  sidecar_t.per_out_S_ij);
                    pkv("S_ijk", sidecar_t.per_out_S_ijk);
                    if (sidecar_t.skip.count(out_name)) {
                        std::fprintf(ft, "%s\"skip\": true", ko_first ? "" : ", ");
                        ko_first = false;
                    }
                    std::fprintf(ft, "}");
                }
                std::fprintf(ft, "\n  }\n}\n");
                std::fclose(ft);
                std::printf("Targets sidecar: %s\n", targets_path.c_str());
            }
        }
    }

    // ---- Write summary ----
    std::FILE* fs = std::fopen(summary_txt.c_str(), "w");
    if (fs) {
        std::fprintf(fs, "rocket_6dof Sobol sensitivity summary\n");
        std::fprintf(fs, "=====================================\n");
        std::fprintf(fs, "Config: %s   seed=%lu   N=%d   total sims=%d\n",
                     config_path.c_str(), seed, N_base, total_runs);
        if (n_boot > 0) {
            std::fprintf(fs, "Bootstrap: %d resamples, %.0f%% CI "
                         "(brackets below denote [lo, hi] percentile range)\n",
                         n_boot, 100.0 * ci_level);
        }
        std::fprintf(fs, "\n");

        // When the CI level is non-default (anything other than 95%),
        // each bracket inline-prefixes the level (e.g. "[99%: +0.31,
        // +0.50]") so that a single quoted row is self-describing.  At
        // the default 95% we leave brackets bare to preserve byte-
        // identical summary output for existing diff workflows.
        const bool show_level_inline = (n_boot > 0)
                                    && (std::fabs(ci_level - 0.95) > 1e-9);
        char lvl_prefix[16] = "";
        if (show_level_inline) {
            std::snprintf(lvl_prefix, sizeof(lvl_prefix), "%.0f%%: ",
                          100.0 * ci_level);
        }

        for (int o = 0; o < n_out; ++o) {
            const auto& I = per_out[o];
            std::fprintf(fs, "Output: %s   mean=%.6g   var=%.6g\n",
                         output_names[o].c_str(), I.mean, I.var);
            if (n_boot > 0) {
                char s1_lbl[32], st_lbl[32];
                std::snprintf(s1_lbl, sizeof(s1_lbl), "S1 [%.0f%% CI]", 100.0 * ci_level);
                std::snprintf(st_lbl, sizeof(st_lbl), "ST [%.0f%% CI]", 100.0 * ci_level);
                std::fprintf(fs, "  %-32s  %18s  %18s\n", "input", s1_lbl, st_lbl);
            } else {
                std::fprintf(fs, "  %-32s  %9s  %9s\n", "input", "S1", "ST");
            }
            // Sort by ST descending for readability
            std::vector<int> order(k);
            for (int i = 0; i < k; ++i) order[i] = i;
            std::sort(order.begin(), order.end(),
                      [&](int a, int b) { return I.ST[a] > I.ST[b]; });
            for (int idx : order) {
                if (n_boot > 0) {
                    std::fprintf(fs, "  %-32s  %5.2f [%s%+.2f,%+.2f]   %5.2f [%s%+.2f,%+.2f]\n",
                                 spec_list[idx].path.c_str(),
                                 I.S1[idx], lvl_prefix, per_out_ci[o].S1_lo[idx], per_out_ci[o].S1_hi[idx],
                                 I.ST[idx], lvl_prefix, per_out_ci[o].ST_lo[idx], per_out_ci[o].ST_hi[idx]);
                } else {
                    std::fprintf(fs, "  %-32s  %9.4f  %9.4f\n",
                                 spec_list[idx].path.c_str(), I.S1[idx], I.ST[idx]);
                }
            }
            // Append group indices
            if (n_groups > 0) {
                std::fprintf(fs, "  -- groups --\n");
                std::vector<int> g_order(n_groups);
                for (int gi = 0; gi < n_groups; ++gi) g_order[gi] = gi;
                std::sort(g_order.begin(), g_order.end(),
                          [&](int a, int b) {
                              return per_out_group[o].ST[a] > per_out_group[o].ST[b];
                          });
                for (int gi : g_order) {
                    std::string lbl = "[group] " + groups[gi].name;
                    std::fprintf(fs, "  %-32s  %9.4f  %9.4f\n",
                                 lbl.c_str(),
                                 per_out_group[o].S1[gi], per_out_group[o].ST[gi]);
                }
            }
            // Append second-order pair indices (S_ij), sorted by |S_ij|
            // descending, top 5 only for readability.
            if (do_second_order && n_pairs > 0) {
                std::fprintf(fs, "  -- second-order pairs (top by |S_ij|) --\n");
                std::vector<int> p_order(n_pairs);
                for (int p = 0; p < n_pairs; ++p) p_order[p] = p;
                std::sort(p_order.begin(), p_order.end(),
                          [&](int a, int b) {
                              return std::fabs(per_out_pairs[o].Sij[a])
                                   > std::fabs(per_out_pairs[o].Sij[b]);
                          });
                int show = std::min(5, n_pairs);
                for (int idx = 0; idx < show; ++idx) {
                    int p = p_order[idx];
                    int i = pair_ij[p].first;
                    int j = pair_ij[p].second;
                    std::string lbl = "[pair] " + spec_list[i].path + " * " + spec_list[j].path;
                    if (n_boot > 0) {
                        std::fprintf(fs, "  %-60s  %+8.3f [%s%+7.3f,%+7.3f]\n",
                                     lbl.c_str(),
                                     per_out_pairs[o].Sij[p],
                                     lvl_prefix,
                                     per_out_ci_pairs[o].Sij_lo[p],
                                     per_out_ci_pairs[o].Sij_hi[p]);
                    } else {
                        std::fprintf(fs, "  %-60s  %9.4f\n",
                                     lbl.c_str(), per_out_pairs[o].Sij[p]);
                    }
                }
            }
            // Append group-pair indices, sorted by |S_GaGb| descending.
            if (n_group_pairs > 0) {
                std::fprintf(fs, "  -- group pairs (by |S_GaGb|) --\n");
                std::vector<int> gp_order(n_group_pairs);
                for (int p = 0; p < n_group_pairs; ++p) gp_order[p] = p;
                std::sort(gp_order.begin(), gp_order.end(),
                          [&](int a, int b) {
                              return std::fabs(per_out_gpair[o].SGij[a])
                                   > std::fabs(per_out_gpair[o].SGij[b]);
                          });
                for (int idx : gp_order) {
                    int ga = group_pair_ij[idx].first;
                    int gb = group_pair_ij[idx].second;
                    std::string lbl = "[gpair] " + groups[ga].name + " * " + groups[gb].name;
                    if (n_boot > 0) {
                        std::fprintf(fs, "  %-40s  %+8.3f [%s%+7.3f,%+7.3f]\n",
                                     lbl.c_str(),
                                     per_out_gpair[o].SGij[idx],
                                     lvl_prefix,
                                     per_out_ci_gpair[o].SGij_lo[idx],
                                     per_out_ci_gpair[o].SGij_hi[idx]);
                    } else {
                        std::fprintf(fs, "  %-40s  %9.4f\n",
                                     lbl.c_str(), per_out_gpair[o].SGij[idx]);
                    }
                }
            }
            // Append third-order triplet indices, sorted by |S_ijk|, top 5.
            if (do_third_order && n_triplets > 0) {
                std::fprintf(fs, "  -- third-order triplets (top by |S_ijk|) --\n");
                std::vector<int> t_order(n_triplets);
                for (int t = 0; t < n_triplets; ++t) t_order[t] = t;
                std::sort(t_order.begin(), t_order.end(),
                          [&](int a, int b) {
                              return std::fabs(per_out_trip[o].Sijk[a])
                                   > std::fabs(per_out_trip[o].Sijk[b]);
                          });
                int show = std::min(5, n_triplets);
                for (int idx = 0; idx < show; ++idx) {
                    int t = t_order[idx];
                    std::string lbl = "[trip] "
                        + spec_list[trip_ijk[t].i].path + " * "
                        + spec_list[trip_ijk[t].j].path + " * "
                        + spec_list[trip_ijk[t].k_].path;
                    if (n_boot > 0) {
                        std::fprintf(fs, "  %-80s  %+8.3f [%s%+7.3f,%+7.3f]\n",
                                     lbl.c_str(),
                                     per_out_trip[o].Sijk[t],
                                     lvl_prefix,
                                     per_out_ci_trip[o].Sijk_lo[t],
                                     per_out_ci_trip[o].Sijk_hi[t]);
                    } else {
                        std::fprintf(fs, "  %-80s  %9.4f\n",
                                     lbl.c_str(), per_out_trip[o].Sijk[t]);
                    }
                }
            }
            std::fprintf(fs, "\n");
        }
        std::fclose(fs);
    }

    // ---- Write group CSV if any groups were defined ----
    if (n_groups > 0) {
        std::string group_csv = indices_csv;
        size_t dot = group_csv.find_last_of('.');
        if (dot != std::string::npos) {
            group_csv = group_csv.substr(0, dot) + "_groups" + group_csv.substr(dot);
        } else {
            group_csv += "_groups";
        }
        std::FILE* fg = std::fopen(group_csv.c_str(), "w");
        if (fg) {
            std::fprintf(fg, "output,group,S1,ST,output_mean,output_variance\n");
            for (int o = 0; o < n_out; ++o) {
                for (int gi = 0; gi < n_groups; ++gi) {
                    std::fprintf(fg, "%s,%s,%.9g,%.9g,%.9g,%.9g\n",
                                 output_names[o].c_str(), groups[gi].name.c_str(),
                                 per_out_group[o].S1[gi], per_out_group[o].ST[gi],
                                 per_out[o].mean, per_out[o].var);
                }
            }
            std::fclose(fg);
            std::printf("Group CSV:     %s\n", group_csv.c_str());
        }
    }

    // ---- Write second-order CSV if second_order enabled ----
    if (do_second_order && n_pairs > 0) {
        std::string pair_csv = indices_csv;
        size_t dot = pair_csv.find_last_of('.');
        if (dot != std::string::npos) {
            pair_csv = pair_csv.substr(0, dot) + "_pairs" + pair_csv.substr(dot);
        } else {
            pair_csv += "_pairs";
        }
        std::FILE* fp = std::fopen(pair_csv.c_str(), "w");
        if (fp) {
            // When bootstrap is enabled, emit CI columns; otherwise keep
            // the legacy 6-column schema for backward compatibility.
            bool with_ci = (n_boot > 0);
            if (with_ci) {
                std::fprintf(fp, "output,input_i,input_j,S_ij,S_ij_lo,S_ij_hi,output_mean,output_variance\n");
            } else {
                std::fprintf(fp, "output,input_i,input_j,S_ij,output_mean,output_variance\n");
            }
            for (int o = 0; o < n_out; ++o) {
                for (int p = 0; p < n_pairs; ++p) {
                    if (with_ci) {
                        std::fprintf(fp, "%s,%s,%s,%.9g,%.9g,%.9g,%.9g,%.9g\n",
                                     output_names[o].c_str(),
                                     spec_list[pair_ij[p].first].path.c_str(),
                                     spec_list[pair_ij[p].second].path.c_str(),
                                     per_out_pairs[o].Sij[p],
                                     per_out_ci_pairs[o].Sij_lo[p],
                                     per_out_ci_pairs[o].Sij_hi[p],
                                     per_out[o].mean, per_out[o].var);
                    } else {
                        std::fprintf(fp, "%s,%s,%s,%.9g,%.9g,%.9g\n",
                                     output_names[o].c_str(),
                                     spec_list[pair_ij[p].first].path.c_str(),
                                     spec_list[pair_ij[p].second].path.c_str(),
                                     per_out_pairs[o].Sij[p],
                                     per_out[o].mean, per_out[o].var);
                    }
                }
            }
            std::fclose(fp);
            std::printf("Pairs CSV:     %s\n", pair_csv.c_str());
        }
    }

    // ---- Write group-pair CSV if group_pairs enabled ----
    if (n_group_pairs > 0) {
        std::string gpair_csv = indices_csv;
        size_t dot = gpair_csv.find_last_of('.');
        if (dot != std::string::npos) {
            gpair_csv = gpair_csv.substr(0, dot) + "_group_pairs" + gpair_csv.substr(dot);
        } else {
            gpair_csv += "_group_pairs";
        }
        std::FILE* fp = std::fopen(gpair_csv.c_str(), "w");
        if (fp) {
            bool with_ci = (n_boot > 0);
            if (with_ci) {
                std::fprintf(fp, "output,group_a,group_b,S_GaGb,S_GaGb_lo,S_GaGb_hi,output_mean,output_variance\n");
            } else {
                std::fprintf(fp, "output,group_a,group_b,S_GaGb,output_mean,output_variance\n");
            }
            for (int o = 0; o < n_out; ++o) {
                for (int p = 0; p < n_group_pairs; ++p) {
                    if (with_ci) {
                        std::fprintf(fp, "%s,%s,%s,%.9g,%.9g,%.9g,%.9g,%.9g\n",
                                     output_names[o].c_str(),
                                     groups[group_pair_ij[p].first].name.c_str(),
                                     groups[group_pair_ij[p].second].name.c_str(),
                                     per_out_gpair[o].SGij[p],
                                     per_out_ci_gpair[o].SGij_lo[p],
                                     per_out_ci_gpair[o].SGij_hi[p],
                                     per_out[o].mean, per_out[o].var);
                    } else {
                        std::fprintf(fp, "%s,%s,%s,%.9g,%.9g,%.9g\n",
                                     output_names[o].c_str(),
                                     groups[group_pair_ij[p].first].name.c_str(),
                                     groups[group_pair_ij[p].second].name.c_str(),
                                     per_out_gpair[o].SGij[p],
                                     per_out[o].mean, per_out[o].var);
                    }
                }
            }
            std::fclose(fp);
            std::printf("Group pairs CSV: %s\n", gpair_csv.c_str());
        }
    }

    // ---- Write triplet CSV if third_order enabled ----
    if (do_third_order && n_triplets > 0) {
        std::string trip_csv = indices_csv;
        size_t dot = trip_csv.find_last_of('.');
        if (dot != std::string::npos) {
            trip_csv = trip_csv.substr(0, dot) + "_triplets" + trip_csv.substr(dot);
        } else {
            trip_csv += "_triplets";
        }
        std::FILE* fp = std::fopen(trip_csv.c_str(), "w");
        if (fp) {
            bool with_ci = (n_boot > 0);
            if (with_ci) {
                std::fprintf(fp, "output,input_i,input_j,input_k,S_ijk,S_ijk_lo,S_ijk_hi,output_mean,output_variance\n");
            } else {
                std::fprintf(fp, "output,input_i,input_j,input_k,S_ijk,output_mean,output_variance\n");
            }
            for (int o = 0; o < n_out; ++o) {
                for (int t = 0; t < n_triplets; ++t) {
                    if (with_ci) {
                        std::fprintf(fp, "%s,%s,%s,%s,%.9g,%.9g,%.9g,%.9g,%.9g\n",
                                     output_names[o].c_str(),
                                     spec_list[trip_ijk[t].i].path.c_str(),
                                     spec_list[trip_ijk[t].j].path.c_str(),
                                     spec_list[trip_ijk[t].k_].path.c_str(),
                                     per_out_trip[o].Sijk[t],
                                     per_out_ci_trip[o].Sijk_lo[t],
                                     per_out_ci_trip[o].Sijk_hi[t],
                                     per_out[o].mean, per_out[o].var);
                    } else {
                        std::fprintf(fp, "%s,%s,%s,%s,%.9g,%.9g,%.9g\n",
                                     output_names[o].c_str(),
                                     spec_list[trip_ijk[t].i].path.c_str(),
                                     spec_list[trip_ijk[t].j].path.c_str(),
                                     spec_list[trip_ijk[t].k_].path.c_str(),
                                     per_out_trip[o].Sijk[t],
                                     per_out[o].mean, per_out[o].var);
                    }
                }
            }
            std::fclose(fp);
            std::printf("Triplets CSV:  %s\n", trip_csv.c_str());
        }
    }

    // ---- Print summary to stdout ----
    std::printf("\n================================================================\n");
    std::printf("SOBOL ANALYSIS COMPLETE\n");
    std::printf("================================================================\n");
    std::printf("Indices CSV:   %s\n", indices_csv.c_str());
    std::printf("Summary:       %s\n", summary_txt.c_str());
    std::printf("\nKey sensitivities (ST > 0.05 only):\n");
    for (int o = 0; o < n_out; ++o) {
        const auto& I = per_out[o];
        // Print one line per output showing the dominant inputs
        std::printf("  %-15s  var=%.4g  ", output_names[o].c_str(), I.var);
        std::vector<int> order(k);
        for (int i = 0; i < k; ++i) order[i] = i;
        std::sort(order.begin(), order.end(),
                  [&](int a, int b) { return I.ST[a] > I.ST[b]; });
        for (int idx : order) {
            if (I.ST[idx] > 0.05) {
                std::printf("%s=%.2f ", spec_list[idx].path.c_str(), I.ST[idx]);
            }
        }
        std::printf("\n");
    }
    if (n_groups > 0) {
        std::printf("\nGroup ST values (S1, ST):\n");
        for (int o = 0; o < n_out; ++o) {
            if (!(per_out[o].var > 0)) continue;
            std::printf("  %-15s  ", output_names[o].c_str());
            for (int gi = 0; gi < n_groups; ++gi) {
                std::printf("%s: S1=%.2f ST=%.2f   ",
                            groups[gi].name.c_str(),
                            per_out_group[o].S1[gi],
                            per_out_group[o].ST[gi]);
            }
            std::printf("\n");
        }
    }
    std::printf("================================================================\n");

    return 0;
}
//  suggest_n: pilot-N advisory mode
//  Strategy:
//    1. Read user config.  Find sobol.target_ci_width (default 0.05)
//       and SOBOL_PILOT_N env var (default 64).
//    2. Write a temp config with sobol.n_base = pilot_N and
//       sobol.n_bootstrap = 500 (force-enable so CIs come out).
//       Disable second_order to save time (we're tuning N from S1/ST).
//    3. Call run_sobol on the temp config, writing to scratch CSV.
//    4. Parse the indices CSV.  For every (output, input) where
//       output_variance > 0, compute the ST CI half-width.  Track
//       the worst (largest) value.
//    5. Recommend N_prod = pilot_N * (worst_half_width / target)^2.
//
//  Bootstrap CI half-widths scale as 1/sqrt(N).  So if pilot at N0
//  gives half-width w0 and you want target w, then:
//    w_N = w0 * sqrt(N0/N)  =>  N = N0 * (w0/w)^2
//
//  This is an estimate, not a guarantee; the constant in 1/sqrt(N)
//  varies with output and input.  Worst-case across all (output,
//  input) pairs is what we report so the recommendation is
//  conservative.
int suggest_n(const std::string& config_path) {
    // Parse to read sobol.target_ci_width (object or scalar form).
    json::Value cfg_check = json::parse_file(config_path);
    TargetWidths targets = read_targets(cfg_check["sobol"]);
    // For --suggest-n: if user didn't set any target, default to ST=0.05
    // (legacy behavior).  Embedded mode is stricter and requires an
    // explicit target.
    if (!targets.any()) {
        targets.ST = 0.05;
        std::fprintf(stderr,
            "Note: no sobol.target_ci_width specified, defaulting to ST=0.05\n\n");
    }

    int pilot_N = 64;
    if (const char* env = std::getenv("SOBOL_PILOT_N")) {
        int v = std::atoi(env);
        if (v > 0) pilot_N = v;
    }

    std::printf("================================================================\n");
    std::printf("    Sobol N suggestion (pilot mode)\n");
    std::printf("================================================================\n");
    std::printf("Pilot N:       %d  (override with SOBOL_PILOT_N env var)\n", pilot_N);
    if (targets.has_ST())    std::printf("Target ST    half-width:  %.4g\n", targets.ST);
    if (targets.has_S1())    std::printf("Target S1    half-width:  %.4g\n", targets.S1);
    if (targets.has_S_ij())  std::printf("Target S_ij  half-width:  %.4g  (pilot will enable second_order)\n", targets.S_ij);
    if (targets.has_S_ijk()) std::printf("Target S_ijk half-width:  %.4g  (pilot will enable third_order)\n", targets.S_ijk);
    std::printf("================================================================\n\n");

    PilotResult pr;
    int rc = run_pilot_and_measure(config_path, pilot_N, targets, pr);
    if (rc != 0) return rc;
    if (pr.n_rows == 0) {
        std::fprintf(stderr, "ERROR: no usable rows in pilot CSV\n");
        return 1;
    }

    // Recommendation strategy.
    //   SOBOL_RECOMMEND_FROM=p90      (default; balanced)
    //   SOBOL_RECOMMEND_FROM=worst    (conservative; max over all rows)
    //   SOBOL_RECOMMEND_FROM=median   (forgiving; covers half the rows)
    //   SOBOL_RECOMMEND_FROM=auto     (adaptive: picks based on worst/p90 gap
    //                                  per estimator -- ≤3x: worst, ≤10x: p90,
    //                                  >10x: median.  May produce a different
    //                                  strategy per estimator.)
    enum class Strategy { P90, WORST, MEDIAN, AUTO };
    Strategy strat = Strategy::P90;
    if (const char* env = std::getenv("SOBOL_RECOMMEND_FROM")) {
        std::string s = env;
        if      (s == "worst")  strat = Strategy::WORST;
        else if (s == "median") strat = Strategy::MEDIAN;
        else if (s == "auto")   strat = Strategy::AUTO;
        else if (s == "p90")    strat = Strategy::P90;
    }
    const char* strategy =
        strat == Strategy::WORST  ? "worst"  :
        strat == Strategy::MEDIAN ? "median" :
        strat == Strategy::AUTO   ? "auto"   : "p90";
    // Helper: when strategy is auto, picks worst/p90/median per estimator
    // based on the worst/p90 gap.  Returns the chosen N and a label.
    auto pick_adaptive = [](double worst_N, double p90_N, double median_N)
        -> std::pair<double, const char*>
    {
        if (worst_N <= 0 && p90_N <= 0) return {0, "n/a"};
        double gap = (p90_N > 0) ? (worst_N / p90_N) : 1.0;
        if (gap <= 3.0)  return {worst_N,  "worst (gap≤3x)"};
        if (gap <= 10.0) return {p90_N,    "p90 (gap≤10x)"};
        return {median_N, "median (gap>10x)"};
    };

    std::printf("\n================================================================\n");
    std::printf("    PILOT RESULTS\n");
    std::printf("================================================================\n");
    std::printf("Rows analyzed:  %d first-order (output, input) pairs\n", pr.n_rows);
    if (pr.n_skipped > 0) {
        std::printf("Rows skipped:   %d  (estimator unstable or output in skip list)\n",
                    pr.n_skipped);
    }
    if (targets.has_S_ij()) {
        std::printf("Pair rows analyzed:     %d\n", pr.n_pair_rows);
    }
    if (targets.has_S_ijk()) {
        std::printf("Triplet rows analyzed:  %d\n", pr.n_trip_rows);
    }
    if (!targets.skip.empty()) {
        std::printf("Skipped outputs (per_output.skip): ");
        bool first = true;
        for (const auto& s : targets.skip) {
            std::printf("%s%s", first ? "" : ", ", s.c_str());
            first = false;
        }
        std::printf("\n");
    }
    if (!targets.per_out_ST.empty() || !targets.per_out_S1.empty()
        || !targets.per_out_S_ij.empty() || !targets.per_out_S_ijk.empty()) {
        std::printf("Per-output target overrides:  active\n");
    }
    std::printf("\nRecommendation strategy: %s "
                "(override with SOBOL_RECOMMEND_FROM=worst|p90|median|auto)\n\n", strategy);

    int N_rec_overall = 1;
    auto pow2_up = [](double Nf) {
        int N = 1; while (N < Nf) N *= 2; return N;
    };
    // Floor recommendation at pilot_N (no point recommending less than
    // what the pilot already used).  Median can otherwise drop to 1
    // when the median row trivially meets the target.
    auto pow2_floor = [pilot_N, &pow2_up](double Nf) {
        int N = pow2_up(Nf);
        if (N < pilot_N) N = pilot_N;
        return N;
    };

    // For each order, show: worst-width row (for context), median/p90/worst
    // N_needed values, the chosen one, the driving row (which row produced
    // worst_*_N), and a gap diagnostic if worst and p90 disagree by >10x.
    auto report_order = [&](const char* tag, double worst_w, double p90_w,
                             double worst_N, double p90_N, double median_N,
                             const std::string& worst_out,
                             const std::string& worst_inp,
                             const std::string& driver_out,
                             const std::string& driver_inp,
                             double driver_tgt)
    {
        double use_N;
        const char* strat_label;
        if (strat == Strategy::AUTO) {
            auto adp = pick_adaptive(worst_N, p90_N, median_N);
            use_N = adp.first;
            strat_label = adp.second;
        } else {
            use_N = (strat == Strategy::WORST)  ? worst_N
                  : (strat == Strategy::MEDIAN) ? median_N
                  :                                p90_N;
            strat_label = strategy;
        }
        if (use_N <= 0) {
            std::printf("%-6s worst_w=%.4g (on %s / %s), p90_w=%.4g  -- no targeted rows\n",
                        tag, worst_w, worst_out.c_str(), worst_inp.c_str(), p90_w);
            return;
        }
        int N_rec = pow2_floor(use_N);
        std::printf("%-6s worst_w=%.4g (on %s / %s), p90_w=%.4g\n",
                    tag, worst_w, worst_out.c_str(), worst_inp.c_str(), p90_w);
        std::printf("  N_needed:  median=%.0f  p90=%.0f  worst=%.0f\n",
                    median_N, p90_N, worst_N);
        if (!driver_out.empty()) {
            std::printf("  worst N driven by %s / %s (target=%.4g)\n",
                        driver_out.c_str(), driver_inp.c_str(), driver_tgt);
        }
        // Gap diagnostic: flag if worst and p90 disagree sharply.
        if (p90_N > 0 && worst_N / p90_N > 10) {
            std::printf("  ⚠  worst/p90 ratio = %.1fx -- single row dominates; consider per_output.skip\n",
                        worst_N / p90_N);
        }
        std::printf("  using %s N_needed = %.0f  ->  %d (next pow2)\n",
                    strat_label, use_N, N_rec);
        if (N_rec > N_rec_overall) N_rec_overall = N_rec;
    };

    if (targets.has_ST()) {
        report_order("ST:",   pr.worst_ST_w,   pr.p90_ST_w,
                     pr.worst_ST_N,  pr.p90_ST_N,  pr.median_ST_N,
                     pr.worst_ST_out, pr.worst_ST_inp,
                     pr.driver_ST_out, pr.driver_ST_inp, pr.driver_ST_tgt);
    }
    if (targets.has_S1()) {
        report_order("S1:",   pr.worst_S1_w,   pr.p90_S1_w,
                     pr.worst_S1_N,  pr.p90_S1_N,  pr.median_S1_N,
                     pr.worst_S1_out, pr.worst_S1_inp,
                     pr.driver_S1_out, pr.driver_S1_inp, pr.driver_S1_tgt);
    }
    if (targets.has_S_ij() && pr.worst_Sij_w > 0) {
        report_order("S_ij:", pr.worst_Sij_w,  pr.p90_Sij_w,
                     pr.worst_Sij_N, pr.p90_Sij_N, pr.median_Sij_N,
                     pr.worst_Sij_out, pr.worst_Sij_inp,
                     pr.driver_Sij_out, pr.driver_Sij_inp, pr.driver_Sij_tgt);
    }
    if (targets.has_S_ijk() && pr.worst_Sijk_w > 0) {
        report_order("S_ijk:",pr.worst_Sijk_w, pr.p90_Sijk_w,
                     pr.worst_Sijk_N,pr.p90_Sijk_N, pr.median_Sijk_N,
                     pr.worst_Sijk_out, pr.worst_Sijk_inp,
                     pr.driver_Sijk_out, pr.driver_Sijk_inp, pr.driver_Sijk_tgt);
    }

    std::printf("\nRECOMMENDED N (max across targets, next pow2): %d\n", N_rec_overall);
    std::printf("================================================================\n");
    return 0;
}
} // namespace rocket6dof
