//  distributions.cpp  --  Probability distribution implementations
//
//  Each distribution wraps a std::random distribution adapter with our
//  uniform Distribution interface.  Truncation is handled by rejection
//  sampling with a safety cap to avoid infinite loops in pathological
//  cases (very tight bounds far from the mode).

#include "distributions.h"

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace rocket6dof {

namespace {
const int   REJECT_MAX  = 1000;
const double NEG_INF    = -std::numeric_limits<double>::infinity();
const double POS_INF    =  std::numeric_limits<double>::infinity();

// Clamp u into (eps, 1-eps) so quantile() never returns ±inf for the
// tail-heavy distributions (normal, lognormal, exponential).
double clamp_u(double u) {
    const double eps = 1e-12;
    if (u < eps)        return eps;
    if (u > 1.0 - eps)  return 1.0 - eps;
    return u;
}

// Standard normal inverse CDF (Acklam's rational approximation, 2003).
// Accurate to ~1.15e-9 absolute over the central region.  Outside
// p in (1e-15, 1 - 1e-15) the approximation degrades but our caller
// clamps to (1e-12, 1 - 1e-12) so we stay well within usable range.
double inv_normal_cdf(double p) {
    static const double a[] = {
        -3.969683028665376e+01,  2.209460984245205e+02,
        -2.759285104469687e+02,  1.383577518672690e+02,
        -3.066479806614716e+01,  2.506628277459239e+00
    };
    static const double b[] = {
        -5.447609879822406e+01,  1.615858368580409e+02,
        -1.556989798598866e+02,  6.680131188771972e+01,
        -1.328068155288572e+01
    };
    static const double c[] = {
        -7.784894002430293e-03, -3.223964580411365e-01,
        -2.400758277161838e+00, -2.549732539343734e+00,
         4.374664141464968e+00,  2.938163982698783e+00
    };
    static const double d[] = {
         7.784695709041462e-03,  3.224671290700398e-01,
         2.445134137142996e+00,  3.754408661907416e+00
    };
    const double plow  = 0.02425;
    const double phigh = 1.0 - plow;
    double q, r;
    if (p < plow) {
        q = std::sqrt(-2.0 * std::log(p));
        return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
               ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
    }
    if (p <= phigh) {
        q = p - 0.5;
        r = q * q;
        return (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q /
               (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1.0);
    }
    q = std::sqrt(-2.0 * std::log(1.0 - p));
    return -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
            ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
}
} // anon

//  Normal
NormalDist::NormalDist(double mean, double stddev, double low, double high)
    : mean_(mean), stddev_(stddev), low_(low), high_(high),
      truncated_(low > NEG_INF || high < POS_INF)
{
    if (stddev < 0) throw std::runtime_error("normal: stddev must be >= 0");
    if (truncated_ && low >= high)
        throw std::runtime_error("normal: low must be < high when truncated");
}

double NormalDist::sample(std::mt19937& rng) const {
    if (stddev_ == 0.0) return mean_;
    std::normal_distribution<double> d(mean_, stddev_);
    if (!truncated_) return d(rng);
    for (int i = 0; i < REJECT_MAX; i++) {
        double v = d(rng);
        if (v >= low_ && v <= high_) return v;
    }
    // Hit the rejection cap; clamp to the nearer bound
    return (mean_ < low_) ? low_ : (mean_ > high_ ? high_ : mean_);
}

std::string NormalDist::describe() const {
    std::ostringstream o;
    o << "normal(mean=" << mean_ << ", stddev=" << stddev_;
    if (truncated_) o << ", low=" << low_ << ", high=" << high_;
    o << ")";
    return o.str();
}

double NormalDist::quantile(double u) const {
    if (stddev_ == 0.0) return mean_;
    double v = mean_ + stddev_ * inv_normal_cdf(clamp_u(u));
    if (truncated_) {
        if (v < low_)  v = low_;
        if (v > high_) v = high_;
    }
    return v;
}

//  Lognormal
LognormalDist::LognormalDist(double mu_log, double sigma_log,
                             double low, double high)
    : mu_log_(mu_log), sigma_log_(sigma_log), low_(low), high_(high),
      truncated_(low > NEG_INF || high < POS_INF)
{
    if (sigma_log < 0) throw std::runtime_error("lognormal: stddev must be >= 0");
    if (truncated_) {
        if (low_ < 0) low_ = 0;   // lognormal is non-negative
        if (low >= high)
            throw std::runtime_error("lognormal: low must be < high when truncated");
    }
}

double LognormalDist::sample(std::mt19937& rng) const {
    if (sigma_log_ == 0.0) return std::exp(mu_log_);
    std::lognormal_distribution<double> d(mu_log_, sigma_log_);
    if (!truncated_) return d(rng);
    for (int i = 0; i < REJECT_MAX; i++) {
        double v = d(rng);
        if (v >= low_ && v <= high_) return v;
    }
    double median = std::exp(mu_log_);
    return (median < low_) ? low_ : (median > high_ ? high_ : median);
}

std::string LognormalDist::describe() const {
    std::ostringstream o;
    o << "lognormal(mu_log=" << mu_log_ << ", sigma_log=" << sigma_log_;
    if (truncated_) o << ", low=" << low_ << ", high=" << high_;
    o << ")";
    return o.str();
}

double LognormalDist::quantile(double u) const {
    if (sigma_log_ == 0.0) return std::exp(mu_log_);
    double v = std::exp(mu_log_ + sigma_log_ * inv_normal_cdf(clamp_u(u)));
    if (truncated_) {
        if (v < low_)  v = low_;
        if (v > high_) v = high_;
    }
    return v;
}

//  Uniform
UniformDist::UniformDist(double low, double high) : low_(low), high_(high) {
    if (!(low < high))
        throw std::runtime_error("uniform: low must be < high");
}

double UniformDist::sample(std::mt19937& rng) const {
    std::uniform_real_distribution<double> d(low_, high_);
    return d(rng);
}

std::string UniformDist::describe() const {
    std::ostringstream o;
    o << "uniform(low=" << low_ << ", high=" << high_ << ")";
    return o.str();
}

double UniformDist::quantile(double u) const {
    return low_ + clamp_u(u) * (high_ - low_);
}

//  Weibull
WeibullDist::WeibullDist(double scale, double shape, double high)
    : scale_(scale), shape_(shape), high_(high),
      truncated_(high < POS_INF)
{
    if (!(scale > 0)) throw std::runtime_error("weibull: scale must be > 0");
    if (!(shape > 0)) throw std::runtime_error("weibull: shape must be > 0");
    if (truncated_ && !(high > 0))
        throw std::runtime_error("weibull: high must be > 0 when truncated");
}

double WeibullDist::sample(std::mt19937& rng) const {
    std::weibull_distribution<double> d(shape_, scale_);
    if (!truncated_) return d(rng);
    for (int i = 0; i < REJECT_MAX; i++) {
        double v = d(rng);
        if (v >= 0 && v <= high_) return v;
    }
    return high_;
}

std::string WeibullDist::describe() const {
    std::ostringstream o;
    o << "weibull(scale=" << scale_ << ", shape=" << shape_;
    if (truncated_) o << ", high=" << high_;
    o << ")";
    return o.str();
}

double WeibullDist::quantile(double u) const {
    // CDF: F(x) = 1 - exp(-(x/scale)^shape).  Inverse: x = scale * (-ln(1-u))^(1/shape)
    double v = scale_ * std::pow(-std::log(1.0 - clamp_u(u)), 1.0 / shape_);
    if (truncated_ && v > high_) v = high_;
    return v;
}

//  Triangular
TriangularDist::TriangularDist(double low, double mode, double high)
    : low_(low), mode_(mode), high_(high)
{
    if (!(low < high))
        throw std::runtime_error("triangular: low must be < high");
    if (!(mode >= low && mode <= high))
        throw std::runtime_error("triangular: mode must be in [low, high]");
}

double TriangularDist::sample(std::mt19937& rng) const {
    // Inverse CDF method.  Let F_c = (mode - low) / (high - low).  Draw
    // U ~ Uniform(0, 1):
    //   if U < F_c:  X = low + sqrt(U * (high-low) * (mode-low))
    //   else:        X = high - sqrt((1-U) * (high-low) * (high-mode))
    std::uniform_real_distribution<double> u(0.0, 1.0);
    double U = u(rng);
    double range = high_ - low_;
    double F_c = (mode_ - low_) / range;
    if (U < F_c) {
        return low_ + std::sqrt(U * range * (mode_ - low_));
    } else {
        return high_ - std::sqrt((1.0 - U) * range * (high_ - mode_));
    }
}

std::string TriangularDist::describe() const {
    std::ostringstream o;
    o << "triangular(low=" << low_ << ", mode=" << mode_ << ", high=" << high_ << ")";
    return o.str();
}

double TriangularDist::quantile(double u) const {
    // Inverse CDF: with F_c = (mode-low)/(high-low):
    //   u < F_c:  low  + sqrt(u * (high-low) * (mode-low))
    //   u >= F_c: high - sqrt((1-u) * (high-low) * (high-mode))
    double U = clamp_u(u);
    double range = high_ - low_;
    double F_c   = (mode_ - low_) / range;
    if (U < F_c) {
        return low_ + std::sqrt(U * range * (mode_ - low_));
    } else {
        return high_ - std::sqrt((1.0 - U) * range * (high_ - mode_));
    }
}

//  Exponential
ExponentialDist::ExponentialDist(double rate) : rate_(rate) {
    if (!(rate > 0)) throw std::runtime_error("exponential: rate must be > 0");
}

double ExponentialDist::sample(std::mt19937& rng) const {
    std::exponential_distribution<double> d(rate_);
    return d(rng);
}

std::string ExponentialDist::describe() const {
    std::ostringstream o;
    o << "exponential(rate=" << rate_ << ")";
    return o.str();
}

double ExponentialDist::quantile(double u) const {
    return -std::log(1.0 - clamp_u(u)) / rate_;
}

//  Constant
ConstantDist::ConstantDist(double value) : value_(value) {}

double ConstantDist::sample(std::mt19937& /*rng*/) const { return value_; }

std::string ConstantDist::describe() const {
    std::ostringstream o;
    o << "constant(" << value_ << ")";
    return o.str();
}

double ConstantDist::quantile(double /*u*/) const {
    return value_;
}

//  Factory: parse a JSON object into a Distribution
std::unique_ptr<Distribution> make_distribution(const json::Value& spec) {
    if (!spec.isObject())
        throw std::runtime_error("distribution: expected object");

    std::string kind = spec["distribution"].asString("");
    if (kind.empty())
        throw std::runtime_error("distribution: missing 'distribution' field");

    auto get = [&](const char* k, double def) {
        return spec[k].asNumber(def);
    };
    auto require = [&](const char* k) {
        if (!spec[k].isNumber())
            throw std::runtime_error(std::string("distribution ") + kind +
                                     ": missing required field '" + k + "'");
        return spec[k].asNumber();
    };

    if (kind == "normal") {
        double mean   = require("mean");
        double stddev = require("stddev");
        double low    = get("low",  NEG_INF);
        double high   = get("high", POS_INF);
        return std::unique_ptr<Distribution>(new NormalDist(mean, stddev, low, high));
    }
    if (kind == "lognormal") {
        // Accept either "mean" / "stddev" (interpreted as log-space) or
        // explicit "mu_log" / "sigma_log".  ModelCenter uses the former.
        double mu_log    = spec["mu_log"].isNumber() ?
                           spec["mu_log"].asNumber() : require("mean");
        double sigma_log = spec["sigma_log"].isNumber() ?
                           spec["sigma_log"].asNumber() : require("stddev");
        double low    = get("low",  NEG_INF);
        double high   = get("high", POS_INF);
        return std::unique_ptr<Distribution>(new LognormalDist(mu_log, sigma_log, low, high));
    }
    if (kind == "uniform") {
        double low  = require("low");
        double high = require("high");
        return std::unique_ptr<Distribution>(new UniformDist(low, high));
    }
    if (kind == "weibull") {
        double scale = require("scale");
        double shape = require("shape");
        double high  = get("high", POS_INF);
        return std::unique_ptr<Distribution>(new WeibullDist(scale, shape, high));
    }
    if (kind == "triangular") {
        double low  = require("low");
        double mode = require("mode");
        double high = require("high");
        return std::unique_ptr<Distribution>(new TriangularDist(low, mode, high));
    }
    if (kind == "exponential") {
        double rate = require("rate");
        return std::unique_ptr<Distribution>(new ExponentialDist(rate));
    }
    if (kind == "constant") {
        double value = require("value");
        return std::unique_ptr<Distribution>(new ConstantDist(value));
    }
    throw std::runtime_error("distribution: unknown type '" + kind + "'");
}

} // namespace rocket6dof
