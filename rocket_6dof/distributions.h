//  distributions.h  --  Probability distributions for Monte Carlo
//
//  Implements the standard set used by Ansys ModelCenter's
//  Probabilistic Analysis Tool:
//
//    Normal       Gaussian with mean mu, std deviation sigma
//                 Optional truncation: bounds [low, high] reject + resample
//
//    Lognormal    log(X) is Normal(mu_log, sigma_log)
//                 Parameters specified in LOG space (matches ModelCenter)
//                 Optional truncation: bounds [low, high] in real space
//
//    Uniform      flat distribution on [low, high]
//
//    Weibull      scale lambda, shape k.  PDF: (k/lambda) (x/lambda)^(k-1)
//                 exp(-(x/lambda)^k) for x >= 0.  Special cases:
//                 k=1 -> exponential; k=2 -> Rayleigh.
//                 Optional truncation: upper bound only (x >= 0)
//
//    Triangular   piecewise-linear PDF defined by (low, mode, high)
//                 where mode is the peak / most likely value
//
//    Exponential  rate parameter lambda. PDF: lambda * exp(-lambda * x)
//                 Special case of Weibull with shape=1.
//
//    Constant     deterministic; same value every sample.  Useful for
//                 disabling a variable without removing it from the
//                 schema.
//
//  All distributions share a common `sample(rng)` interface that
//  returns a double.  The RNG is std::mt19937 seeded externally.
//
//  Truncated distributions use rejection sampling: draw repeatedly
//  until the sample is in [low, high].  After 1000 rejected draws we
//  give up and return the closest boundary (extremely unlikely unless
//  the bounds are pathologically tight).
//
//  Distribution can be constructed from a parsed JSON object:
//
//    {
//      "distribution": "normal",
//      "mean": 500.0,
//      "stddev": 10.0,
//      "low": 480.0,   // optional truncation
//      "high": 520.0
//    }

#ifndef ROCKET6DOF_DISTRIBUTIONS_H
#define ROCKET6DOF_DISTRIBUTIONS_H

#include "json.h"

#include <memory>
#include <random>
#include <string>

namespace rocket6dof {

class Distribution {
public:
    virtual ~Distribution() = default;
    virtual double sample(std::mt19937& rng) const = 0;
    // Inverse CDF: map u in (0,1) to a physical value.  Required for
    // quasi-random / Latin-hypercube / Sobol-sequence sampling where
    // we generate uniform samples in [0,1] and need to transform them
    // through the variable's distribution.  u=0 maps to the lower
    // support boundary, u=1 to the upper boundary.  Implementations
    // should clamp inputs to (eps, 1-eps) to avoid infinities at the
    // tails.
    virtual double quantile(double u) const = 0;
    virtual std::string kind() const = 0;
    virtual std::string describe() const = 0;
};

// Construct a Distribution from a JSON spec.  Throws std::runtime_error
// with a descriptive message if the spec is missing fields or has an
// unknown distribution type.  On success, the unique_ptr owns the
// distribution and the caller can keep it for the lifetime of the MC
// run.
std::unique_ptr<Distribution> make_distribution(const json::Value& spec);

// ---- Concrete distributions, exposed for direct construction and
// unit tests.  In normal use, prefer make_distribution(). ----

class NormalDist : public Distribution {
public:
    NormalDist(double mean, double stddev, double low, double high);
    double sample(std::mt19937& rng) const override;
    double quantile(double u) const override;
    std::string kind()     const override { return "normal"; }
    std::string describe() const override;
private:
    double mean_, stddev_, low_, high_;
    bool   truncated_;
};

class LognormalDist : public Distribution {
public:
    // mu_log, sigma_log are mean and std-dev of LOG of the variable.
    // low/high are in REAL space; pass -inf/+inf for un-truncated.
    LognormalDist(double mu_log, double sigma_log, double low, double high);
    double sample(std::mt19937& rng) const override;
    double quantile(double u) const override;
    std::string kind()     const override { return "lognormal"; }
    std::string describe() const override;
private:
    double mu_log_, sigma_log_, low_, high_;
    bool   truncated_;
};

class UniformDist : public Distribution {
public:
    UniformDist(double low, double high);
    double sample(std::mt19937& rng) const override;
    double quantile(double u) const override;
    std::string kind()     const override { return "uniform"; }
    std::string describe() const override;
private:
    double low_, high_;
};

class WeibullDist : public Distribution {
public:
    // scale = lambda > 0, shape = k > 0.  high is optional upper-truncation.
    WeibullDist(double scale, double shape, double high);
    double sample(std::mt19937& rng) const override;
    double quantile(double u) const override;
    std::string kind()     const override { return "weibull"; }
    std::string describe() const override;
private:
    double scale_, shape_, high_;
    bool   truncated_;
};

class TriangularDist : public Distribution {
public:
    TriangularDist(double low, double mode, double high);
    double sample(std::mt19937& rng) const override;
    double quantile(double u) const override;
    std::string kind()     const override { return "triangular"; }
    std::string describe() const override;
private:
    double low_, mode_, high_;
};

class ExponentialDist : public Distribution {
public:
    explicit ExponentialDist(double rate);
    double sample(std::mt19937& rng) const override;
    double quantile(double u) const override;
    std::string kind()     const override { return "exponential"; }
    std::string describe() const override;
private:
    double rate_;
};

class ConstantDist : public Distribution {
public:
    explicit ConstantDist(double value);
    double sample(std::mt19937& rng) const override;
    double quantile(double u) const override;
    std::string kind()     const override { return "constant"; }
    std::string describe() const override;
private:
    double value_;
};

} // namespace rocket6dof

#endif
