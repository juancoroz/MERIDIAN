//  distributions_test.cpp  --  Verify the MC distribution machinery
//
//  Statistical tests sample N points and check that the empirical
//  moments (mean, std deviation) match analytical expectations within
//  a tolerance that scales like 1/sqrt(N).  N = 50,000 gives a
//  tolerance of about 0.5%; we use 1% for safety margin against
//  spurious failures.

#include "distributions.h"
#include "json.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using namespace rocket6dof;

namespace {

int fails = 0;
bool check(bool cond, const char* desc) {
    std::printf("  %s %s\n", cond ? "OK  " : "FAIL", desc);
    if (!cond) fails++;
    return cond;
}

bool close(double a, double b, double tol) {
    return std::fabs(a - b) <= tol;
}

void empirical_moments(const std::vector<double>& xs, double& mean, double& stddev) {
    double s = 0;
    for (double v : xs) s += v;
    mean = s / xs.size();
    double s2 = 0;
    for (double v : xs) s2 += (v - mean) * (v - mean);
    stddev = std::sqrt(s2 / xs.size());
}

std::vector<double> sample_n(const Distribution& d, int n, unsigned seed = 12345) {
    std::mt19937 rng(seed);
    std::vector<double> out(n);
    for (int i = 0; i < n; ++i) out[i] = d.sample(rng);
    return out;
}

int test_normal() {
    std::printf("\n=== Test 1: Normal distribution ===\n");
    NormalDist d(10.0, 2.0, -1e300, 1e300);
    auto xs = sample_n(d, 50000);
    double m, s;
    empirical_moments(xs, m, s);
    std::printf("  empirical mean=%.4f (expect 10.0), stddev=%.4f (expect 2.0)\n", m, s);
    check(close(m, 10.0, 0.05), "mean within 0.05 of 10.0");
    check(close(s,  2.0, 0.05), "stddev within 0.05 of 2.0");

    // Truncated
    NormalDist t(10.0, 2.0, 8.0, 12.0);
    auto ts = sample_n(t, 10000);
    double tmin = *std::min_element(ts.begin(), ts.end());
    double tmax = *std::max_element(ts.begin(), ts.end());
    check(tmin >= 8.0,  "truncated: min in [low, high]");
    check(tmax <= 12.0, "truncated: max in [low, high]");

    // Determinism: same seed, same sample
    std::mt19937 r1(42), r2(42);
    check(d.sample(r1) == d.sample(r2), "same seed -> same sample");

    // zero std-dev returns mean
    NormalDist z(7.0, 0.0, -1e300, 1e300);
    std::mt19937 zr(99);
    check(z.sample(zr) == 7.0, "stddev=0 -> always returns mean");
    return 0;
}

int test_lognormal() {
    std::printf("\n=== Test 2: Lognormal distribution ===\n");
    // mu_log=0, sigma_log=0.5
    // Theoretical mean of X = exp(mu + sigma^2/2) = exp(0.125) = 1.1331
    // Theoretical variance = (exp(sigma^2) - 1) * exp(2*mu + sigma^2)
    //                      = (exp(0.25) - 1) * exp(0.25) = 0.2840 * 1.2840 = 0.3647
    // stddev = sqrt(0.3647) = 0.6039
    LognormalDist d(0.0, 0.5, -1e300, 1e300);
    auto xs = sample_n(d, 100000);
    double m, s;
    empirical_moments(xs, m, s);
    std::printf("  empirical mean=%.4f (expect 1.1331), stddev=%.4f (expect 0.6039)\n", m, s);
    check(close(m, 1.1331, 0.02), "mean ~ exp(mu + sigma^2/2)");
    check(close(s, 0.6039, 0.03), "stddev matches analytical");

    // Always positive
    int negative = 0;
    for (double v : xs) if (v < 0) negative++;
    check(negative == 0, "all samples >= 0");

    return 0;
}

int test_uniform() {
    std::printf("\n=== Test 3: Uniform distribution ===\n");
    UniformDist d(5.0, 15.0);
    auto xs = sample_n(d, 50000);
    double m, s;
    empirical_moments(xs, m, s);
    // mean = (low+high)/2 = 10.0
    // stddev = (high-low)/sqrt(12) = 10/sqrt(12) = 2.8868
    std::printf("  empirical mean=%.4f (expect 10.0), stddev=%.4f (expect 2.887)\n", m, s);
    check(close(m, 10.0, 0.05),   "mean = (low+high)/2");
    check(close(s,  2.887, 0.05), "stddev = (high-low)/sqrt(12)");

    // Bounds
    double mn = *std::min_element(xs.begin(), xs.end());
    double mx = *std::max_element(xs.begin(), xs.end());
    check(mn >= 5.0,  "min in range");
    check(mx <= 15.0, "max in range");
    return 0;
}

int test_weibull() {
    std::printf("\n=== Test 4: Weibull distribution ===\n");
    // scale=1, shape=2 -> Rayleigh.
    // mean = scale * Gamma(1 + 1/shape) = Gamma(1.5) = 0.8862
    // var  = scale^2 * (Gamma(1+2/k) - (Gamma(1+1/k))^2)
    //      = 1 * (Gamma(2) - 0.8862^2) = 1 - 0.7855 = 0.2146
    // stddev = 0.4633
    WeibullDist d(1.0, 2.0, 1e300);
    auto xs = sample_n(d, 50000);
    double m, s;
    empirical_moments(xs, m, s);
    std::printf("  Rayleigh: empirical mean=%.4f (expect 0.886), stddev=%.4f (expect 0.463)\n", m, s);
    check(close(m, 0.8862, 0.02), "Rayleigh mean");
    check(close(s, 0.4633, 0.02), "Rayleigh stddev");

    // Shape=1 -> exponential with rate 1/scale
    WeibullDist d2(2.0, 1.0, 1e300);
    auto xs2 = sample_n(d2, 50000);
    empirical_moments(xs2, m, s);
    // mean = scale (for shape=1, this is exponential's mean = 1/rate = scale)
    std::printf("  shape=1: mean=%.4f (expect 2.0), stddev=%.4f (expect 2.0)\n", m, s);
    check(close(m, 2.0, 0.05), "shape=1 -> exponential mean = scale");
    check(close(s, 2.0, 0.05), "shape=1 -> exponential stddev = scale");

    return 0;
}

int test_triangular() {
    std::printf("\n=== Test 5: Triangular distribution ===\n");
    // low=0, mode=5, high=10
    // mean = (low + mode + high) / 3 = 5.0
    // var  = (low^2 + mode^2 + high^2 - low*mode - low*high - mode*high) / 18
    //      = (0 + 25 + 100 - 0 - 0 - 50) / 18 = 75 / 18 = 4.1667
    // stddev = 2.041
    TriangularDist d(0.0, 5.0, 10.0);
    auto xs = sample_n(d, 50000);
    double m, s;
    empirical_moments(xs, m, s);
    std::printf("  empirical mean=%.4f (expect 5.0), stddev=%.4f (expect 2.041)\n", m, s);
    check(close(m, 5.0, 0.05),   "mean = (low+mode+high)/3");
    check(close(s, 2.041, 0.05), "stddev matches analytical");

    // Asymmetric: low=0, mode=2, high=10 -> mean = 4.0
    TriangularDist d2(0.0, 2.0, 10.0);
    auto xs2 = sample_n(d2, 50000);
    empirical_moments(xs2, m, s);
    std::printf("  asymm:    mean=%.4f (expect 4.0)\n", m);
    check(close(m, 4.0, 0.05), "asymmetric: mean = (0+2+10)/3");

    // Bounds
    double mn = *std::min_element(xs.begin(), xs.end());
    double mx = *std::max_element(xs.begin(), xs.end());
    check(mn >= 0.0 && mx <= 10.0, "samples in [low, high]");

    return 0;
}

int test_exponential() {
    std::printf("\n=== Test 6: Exponential distribution ===\n");
    // rate = 0.5 -> mean = 1/rate = 2.0, stddev = 2.0
    ExponentialDist d(0.5);
    auto xs = sample_n(d, 50000);
    double m, s;
    empirical_moments(xs, m, s);
    std::printf("  empirical mean=%.4f (expect 2.0), stddev=%.4f (expect 2.0)\n", m, s);
    check(close(m, 2.0, 0.05), "mean = 1/rate");
    check(close(s, 2.0, 0.05), "stddev = 1/rate");
    return 0;
}

int test_constant() {
    std::printf("\n=== Test 7: Constant distribution ===\n");
    ConstantDist d(42.0);
    std::mt19937 rng(0);
    for (int i = 0; i < 100; ++i) {
        if (d.sample(rng) != 42.0) {
            check(false, "constant always returns the same value");
            return 0;
        }
    }
    check(true, "constant always returns 42.0");
    return 0;
}

int test_factory() {
    std::printf("\n=== Test 8: JSON factory ===\n");

    auto v = json::parse_string(R"({
      "distribution": "normal",
      "mean": 100.0,
      "stddev": 5.0
    })");
    auto d = make_distribution(v);
    check(d != nullptr && d->kind() == "normal", "normal factory");
    std::printf("    -> %s\n", d->describe().c_str());

    v = json::parse_string(R"({
      "distribution": "uniform", "low": -1.0, "high": 1.0
    })");
    d = make_distribution(v);
    check(d->kind() == "uniform", "uniform factory");

    v = json::parse_string(R"({
      "distribution": "lognormal", "mean": 0.0, "stddev": 0.5
    })");
    d = make_distribution(v);
    check(d->kind() == "lognormal", "lognormal factory");

    v = json::parse_string(R"({
      "distribution": "triangular", "low": 0, "mode": 1, "high": 3
    })");
    d = make_distribution(v);
    check(d->kind() == "triangular", "triangular factory");

    v = json::parse_string(R"({ "distribution": "weibull", "scale": 1, "shape": 2 })");
    d = make_distribution(v);
    check(d->kind() == "weibull", "weibull factory");

    v = json::parse_string(R"({ "distribution": "exponential", "rate": 0.1 })");
    d = make_distribution(v);
    check(d->kind() == "exponential", "exponential factory");

    v = json::parse_string(R"({ "distribution": "constant", "value": 3.14 })");
    d = make_distribution(v);
    check(d->kind() == "constant", "constant factory");
    {
        std::mt19937 rng(1);
        check(d->sample(rng) == 3.14, "constant sample == 3.14");
    }

    // Error cases
    bool threw = false;
    try {
        auto bad = json::parse_string(R"({ "distribution": "normal", "mean": 0 })");
        make_distribution(bad);
    } catch (const std::exception& e) {
        threw = true;
        std::printf("    -> caught: %s\n", e.what());
    }
    check(threw, "missing required field throws");

    threw = false;
    try {
        auto bad = json::parse_string(R"({ "distribution": "weird" })");
        make_distribution(bad);
    } catch (const std::exception& e) { threw = true; }
    check(threw, "unknown distribution throws");

    threw = false;
    try {
        UniformDist bad(5.0, 5.0);    // low == high
        (void)bad;
    } catch (const std::exception&) { threw = true; }
    check(threw, "uniform low >= high throws");

    return 0;
}

int test_reproducibility() {
    std::printf("\n=== Test 9: Reproducibility ===\n");

    NormalDist d(0.0, 1.0, -1e300, 1e300);
    auto a = sample_n(d, 1000, 42);
    auto b = sample_n(d, 1000, 42);
    bool match = true;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) { match = false; break; }
    }
    check(match, "same seed -> identical samples");

    auto c = sample_n(d, 1000, 43);
    bool differ = false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != c[i]) { differ = true; break; }
    }
    check(differ, "different seed -> different samples");

    return 0;
}

int test_quantile() {
    std::printf("\n=== Test 10: Quantile (inverse CDF) ===\n");

    // For each distribution, draw N uniforms u_i, transform through
    // quantile(), and check that the resulting empirical moments match
    // the analytical mean/stddev of the distribution within ~1%.
    auto stratify = [](int n, std::mt19937& rng) {
        // Stratified uniforms in [0,1]: better convergence than plain uniform
        std::uniform_real_distribution<double> jit(0.0, 1.0);
        std::vector<double> us(n);
        for (int i = 0; i < n; ++i) us[i] = (i + jit(rng)) / n;
        return us;
    };

    std::mt19937 rng(123);
    const int N = 20000;
    auto us = stratify(N, rng);
    std::vector<double> xs(N);

    // Normal(10, 2)
    NormalDist d(10.0, 2.0, -1e300, 1e300);
    for (int i = 0; i < N; ++i) xs[i] = d.quantile(us[i]);
    double m, s; empirical_moments(xs, m, s);
    std::printf("  Normal(10,2):  q-mean=%.4f, q-stddev=%.4f\n", m, s);
    check(close(m, 10.0, 0.05), "normal quantile mean");
    check(close(s,  2.0, 0.05), "normal quantile stddev");

    // Uniform(5, 15) -- exact reproduction of (5 + u*10)
    UniformDist du(5.0, 15.0);
    check(close(du.quantile(0.0),  5.0,  1e-10), "uniform q(0)  = low");
    check(close(du.quantile(0.5), 10.0,  1e-10), "uniform q(0.5)= mid");
    check(close(du.quantile(1.0), 15.0,  1e-10), "uniform q(1)  = high");

    // Lognormal: median is exp(mu)
    LognormalDist dl(0.0, 0.5, -1e300, 1e300);
    check(close(dl.quantile(0.5), 1.0,  1e-3), "lognormal q(0.5) = exp(mu) = 1.0");
    for (int i = 0; i < N; ++i) xs[i] = dl.quantile(us[i]);
    empirical_moments(xs, m, s);
    std::printf("  Lognormal(0,0.5):  q-mean=%.4f (expect 1.1331)\n", m);
    check(close(m, 1.1331, 0.02), "lognormal quantile mean");

    // Weibull(scale=1, shape=2) - Rayleigh
    WeibullDist dw(1.0, 2.0, 1e300);
    for (int i = 0; i < N; ++i) xs[i] = dw.quantile(us[i]);
    empirical_moments(xs, m, s);
    std::printf("  Weibull(1,2)/Rayleigh:  q-mean=%.4f (expect 0.886)\n", m);
    check(close(m, 0.8862, 0.02), "weibull quantile mean");

    // Triangular(0, 5, 10) -- boundary check via clamp_u: q(0) and
    // q(1) come within sqrt(1e-12)*range ≈ 1e-5 of the actual boundary
    TriangularDist dt(0.0, 5.0, 10.0);
    check(close(dt.quantile(0.0),  0.0, 1e-4), "triangular q(0)  ≈ low");
    check(close(dt.quantile(1.0), 10.0, 1e-4), "triangular q(1)  ≈ high");
    for (int i = 0; i < N; ++i) xs[i] = dt.quantile(us[i]);
    empirical_moments(xs, m, s);
    check(close(m, 5.0, 0.05), "triangular quantile mean = 5.0");

    // Exponential(rate=0.5)
    ExponentialDist de(0.5);
    for (int i = 0; i < N; ++i) xs[i] = de.quantile(us[i]);
    empirical_moments(xs, m, s);
    check(close(m, 2.0, 0.05), "exponential quantile mean = 1/rate");

    // Constant
    ConstantDist dc(7.0);
    check(dc.quantile(0.3) == 7.0, "constant q(any) = value");
    check(dc.quantile(0.9) == 7.0, "constant q(any) = value (2nd check)");

    // Roundtrip: u=0.5 of any symmetric distribution should be the mean
    NormalDist dn2(50.0, 5.0, -1e300, 1e300);
    check(close(dn2.quantile(0.5), 50.0, 1e-6), "Normal q(0.5) = mean (symmetric)");

    return 0;
}

} // anon

int main() {
    test_normal();
    test_lognormal();
    test_uniform();
    test_weibull();
    test_triangular();
    test_exponential();
    test_constant();
    test_factory();
    test_reproducibility();
    test_quantile();
    std::printf("\n=== Total failures: %d ===\n", fails);
    return fails == 0 ? 0 : 1;
}
