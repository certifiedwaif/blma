// priors.cpp

#include <Rcpp.h>
#include <cmath>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// [[Rcpp::depends(RcppGSL)]]
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>

#include "GaussLegendre.h"
#include "priors.h"

using namespace std;
using namespace Rcpp;

using Eigen::VectorXd;
using std::string;

// Which one was this one again? Maruyama?
double maruyama(const int n, const int p, const double R2, int p_gamma)
{
  const auto sigma2 = 1. - R2;
  const auto a = 1.;
  const auto b = p;
  #ifdef DEBUG
  Rcpp::Rcout << "sigma2 " << sigma2 << std::endl;
  Rcpp::Rcout << "n " << n << std::endl;
  #endif
  const auto log_sigma2 = std::log(sigma2);
  const auto log_n = std::log(n);
  #ifdef DEBUG
  Rcpp::Rcout << "log_sigma2 " << log_sigma2 << std::endl;
  Rcpp::Rcout << "p_gamma " << p_gamma << std::endl;
  Rcpp::Rcout << "log_n " << log_n << std::endl;
  Rcpp::Rcout << "Rf_lbeta(" << a + p_gamma << ", " << b + p - p_gamma << ") " << Rf_lbeta(a + p_gamma, b + p - p_gamma) << std::endl;
  #endif
  double log_p;
  log_p = -n / 2. * log_sigma2 - p_gamma / 2. * log_n + Rf_lbeta(a + p_gamma, b + p - p_gamma);
  if (sigma2 == 0. || std::isnan(sigma2)) {
    log_p = -INFINITY;
    // throw std::range_error("Invalid sigma2");
  }
  #ifdef DEBUG
  Rcpp::Rcout << "log_p " << log_p << std::endl;
  #endif
  // if (isnan(log_p)) {
  // // If a floating point number is not equal to itself, then it must be NaN
  // // if (log_p != log_p) {
  //  log_p = -INFINITY;
  // }
  #ifdef DEBUG
  Rcpp::Rcout << "log_p " << log_p << std::endl;
  #endif

  return log_p;
}


double BIC(const int n, const int p, double R2, int vp_gamma)
{
  const auto sigma2 = 1. - R2;
  #ifdef DEBUG
  Rcpp::Rcout << "n " << n << " p " << p << " R2 " << R2 << " vp_gamma " << vp_gamma << std::endl;
  #endif
  auto BIC = n * log(sigma2) + vp_gamma * log(n);
  if (sigma2 == 0. || std::isnan(sigma2)) {
    BIC = -INFINITY;
    // throw std::range_error("Invalid sigma2");
  }
  #ifdef DEBUG
  Rcpp::Rcout << "n * log(1 - R2) " << n * log(1 - R2) << std::endl;
  Rcpp::Rcout << "vp_gamma * log(n) " << vp_gamma * log(n) << std::endl;
  #endif
  return -0.5 * BIC;
}


double ZE(const int n, const int p, double R2, int p_gamma)
{
  auto a = -0.75;
  auto b = 0.5 * (n - p_gamma - 5) - a;
  auto c = 0.5 * (n - 1);
  auto d = 0.5 * p_gamma + a;

  auto log_p = -(b+1)*log(1 - R2) + Rf_lbeta(d+1,b+1) - Rf_lbeta(a+1,b+1);
  auto ZE = -2*log_p;
  return log_p;
}


double log_hyperg_2F1(double b, double c, double x)
{
  if (x == 0.)
    return 0.;
  auto val = 0.;
  val += log(c-1);
  val += (1-c)*log(x);
  val += (c-b-1)*log(1-x);
  val += Rf_lbeta(c-1,b-c+1);
  val += Rf_pbeta(x, (c-1), (b-c+1), true, true);
  return val;
}


double log_hyperg_2F1_naive(double b, double c, double x)
{
  auto val = log(gsl_sf_hyperg_2F1( b, 1, c, x));
  return val;
}


// Liang's hyper g-prior
double liang_g1(const int n, const int p, double R2, int p_gamma)
{
  auto a = 3.;
  double log_p_g;
  log_p_g = log(a - 2) - log(p_gamma + a - 2) + log(gsl_sf_hyperg_2F1(0.5*(n-1), 1, 0.5*(p_gamma + a), R2));
  return log_p_g;
}


// Liang's g prior
double liang_g2(const int n, const int p, double R2, int p_gamma)
{
  auto a = 3.;
  auto log_vp_g2 = log(a - 2) - log(p_gamma + a - 2) + log_hyperg_2F1( 0.5*(n-1), 0.5*(p_gamma+a), R2);
  return log_vp_g2;
}


// Liang's g/n prior Appell
double liang_g_n_appell(const int n, const int p, double R2, int p_gamma)
{
  auto a = 3.;

  Rcpp::Environment appell("package:appell");
  Rcpp::Function appellf1_r = appell["appellf1"];
  Rcpp::ComplexVector val(1);
  try {
    Rcpp::List res = appellf1_r(Rcpp::_["a"] = 1.,
                                Rcpp::_["b1"] = a / 2.,
                                Rcpp::_["b2"] = (n - 1.)/2.,
                                Rcpp::_["c"] = (p_gamma + a) / 2.,
                                Rcpp::_["x"] = 1. - 1. / n,
                                Rcpp::_["y"] = R2,
                                Rcpp::_["userflag"] = 1,
                                Rcpp::_["hyp2f1"] = "michel.stoitsov");
    val(0) = res["val"];
  } catch (...) {
    val = Rcpp::ComplexVector(1);
    val(0).r = NA_REAL;
    val(0).i = NA_REAL;
  }

  auto result = log(a - 2.) - log(n) - log(p_gamma + a - 2.) + log(val(0).r);
  #ifdef DEBUG
  Rcpp::Rcout << "liang_g_n_appell(" << n << ", " << p << ", " << R2 << ", " << p_gamma << ") = " << result << std::endl;
  #endif
  return result;
}


// Trapezoidal integration over a potentially irregular grid
double trapint(const VectorXd& xgrid, const VectorXd& fgrid)
{
  auto sum = 0.0;

  #pragma omp simd reduction(+:sum)
  for (auto i = 0; i < xgrid.size() - 1; i++) {
    sum += 0.5 * (xgrid(i + 1) - xgrid(i)) * (fgrid(i) + fgrid(i + 1));
  }
  // Rcpp::Rcout << "sum " << sum << std::endl;

  return sum;
}


// Liang's g/n prior quadrature
double liang_g_n_quad(const int n, const int p, double R2, int p_gamma)
{
  auto a = 3.;
  const int NUM_POINTS = 10000;
  VectorXd xgrid(NUM_POINTS);
  VectorXd fgrid(NUM_POINTS);
  for (int i = 0; i < NUM_POINTS; i++) {
    double u = static_cast<double>(i) / static_cast<double>(NUM_POINTS);
    xgrid(i) = u;
    fgrid(i) = exp((p_gamma / 2. + a / 2. - 2.) * log(1 - u) + -a/2. * log(1. - u * (1. - 1. / n)) + (-(n-1.)/2.) * log(1 - u*R2));
  }
  auto result = log(a - 2.) - log(2. * n) + log(trapint(xgrid, fgrid));
  #ifdef DEBUG
  Rcpp::Rcout << "liang_g_n_quad(" << n << ", " << p << ", " << R2 << ", " << p_gamma << ") = " << result << std::endl;
  #endif
  return result;
}


// Liang's g/n prior approximation
double liang_g_n_approx(const int n, const int p, double R2, int p_gamma)
{
  // #ifdef DEBUG
  // Rcpp::Rcout << "n " << n << " p " << p << " R2 " << R2 << " p_gamma " << p_gamma << std::endl;
  // #endif
  if (p_gamma == 0)
    return 0.;
  if (p_gamma == 1 || p_gamma == 2)
    return liang_g_n_quad(n, p, R2, p_gamma);

  auto a = 3.;

  auto shape1 = 0.5*(p_gamma - 1.);
  auto shape2 = 0.5*(n - p_gamma + 1.);

  auto log_y = log(a - 2) - log(n) - log(R2) - log(1 - R2);
  log_y = log_y + ::Rf_pbeta(R2,shape1,shape2,true,true);
  log_y = log_y - ::Rf_dbeta(R2,shape1,shape2,true);

  #ifdef DEBUG
  Rcpp::Rcout << "liang_g_n_approx(" << n << ", " << p << ", " << R2 << ", " << p_gamma << ") = " << log_y << std::endl;
  #endif
  return log_y;
}


double robust_bayarri1(const int n, const int p, double R2, int p_gamma)
{
  // Rcpp::Rcout << "n " << n << " R2 " << R2 << " p_gamma " << p_gamma << std::endl;
  double r = (1. + n) / (1. + p_gamma);
  double L = r - 1.;

  const int NUM_POINTS = 10000;
  VectorXd x(NUM_POINTS);
  x.setLinSpaced(NUM_POINTS, L, 10000);

  double sigma2 = 1 - R2;
  double beta   = (1 + sigma2*L)/sigma2;

  VectorXd log_f(NUM_POINTS);
  for (int i = 0; i < NUM_POINTS; i++) {
    log_f(i) = -log(2.) + 0.5 * log(r) + 0.5 * (n - p_gamma - 4.) * log(1. + x(i)) - 0.5 * (n - 1.) * log(1. + x(i) * (1. - R2));
    // log_f(i) = -log(2.) + 0.5 * log(r) - 0.5 * (n - 1.)*log(sigma2) + 0.5 * (n - p_gamma - 4.) * log(r + x[i]) - 0.5 * (n - 1.) * log(beta + x[i]);
    // #ifdef DEBUG
    // Rcpp::Rcout << "-log(2) " << -log(2) << std::endl;
    // Rcpp::Rcout << "0.5 * log(r) " << 0.5 * log(r) << std::endl;
    // Rcpp::Rcout << "-0.5 * (n - 1) * log(sigma2) " << -0.5 * (n - 1) * log(sigma2) << std::endl;
    // Rcpp::Rcout << "0.5 * (n - p_gamma - 4) * log(r + x[i]) " << 0.5 * (n - p_gamma - 4) * log(r + x[i]) << std::endl;
    // Rcpp::Rcout << "0.5 * (n - 1.) * log(beta + x[i])  " << 0.5 * (n - 1.) * log(beta + x[i])  << std::endl;
    // // if (i < 10)
    // //   Rcpp::Rcout << "x(" << i << ") " << x(i) << " log_f(i) " << log_f(i) << std::endl;
    // #endif
  }
  double result = log(trapint(x, log_f.array().exp()));
  // Rcpp::Rcout << "result " << result << std::endl;
  return result;
}


double robust_bayarri2(const int n, const int p, double R2, int p_gamma)
{
  #ifdef DEBUG
  Rcpp::Rcout << "n " << n;
  Rcpp::Rcout << " p " << p;
  Rcpp::Rcout << " R2 " << R2;
  Rcpp::Rcout << " p_gamma " << p_gamma;
  #endif
  auto sigma2 = 1. - R2;
  auto L = (1. + n)/(1. + p_gamma) - 1.;
  auto z = R2/(1. + L*sigma2);

  if (p_gamma == 0)
    return 0.;
  double log_vp_gprior7 = 0.5*(n - p_gamma - 1)*log( n + 1 );
  #ifdef DEBUG
  Rcpp::Rcout << " log_vp_gprior7 1 " << log_vp_gprior7 << std::endl;
  #endif
  log_vp_gprior7 -= 0.5*(n - p_gamma - 1)*log( p_gamma + 1);
  #ifdef DEBUG
  Rcpp::Rcout << " log_vp_gprior7 1 " << log_vp_gprior7 << std::endl;
  #endif
  log_vp_gprior7 -= 0.5*(n - 1)*log(1 + L*sigma2);
  #ifdef DEBUG
  Rcpp::Rcout << " log_vp_gprior7 1 " << log_vp_gprior7 << std::endl;
  #endif
  double R2_gamma_tilde = R2 / (1. + L * sigma2);
  log_vp_gprior7 -= log(2 * R2_gamma_tilde * (1. - R2_gamma_tilde));
  #ifdef DEBUG
  Rcpp::Rcout << " log_vp_gprior7 1 " << log_vp_gprior7 << std::endl;
  #endif
  log_vp_gprior7 += ::Rf_pbeta(R2_gamma_tilde, 0.5 * (p_gamma + 1.), 0.5 * (n - p_gamma - 2.), true, true) - ::Rf_dbeta(R2_gamma_tilde, 0.5 * (p_gamma + 1.), 0.5 * (n - p_gamma - 2.), true);
  #ifdef DEBUG
  Rcpp::Rcout << " log_vp_gprior7 1 " << log_vp_gprior7 << std::endl;
  #endif
  return log_vp_gprior7;
}

double log_BF_g_on_n_integrand (double vu, int n, int p, double R2, double a)
{
  double vals = 0.;
  vals += log (a - 2);
  vals -= log (2 * n);
  vals += 0.5 * (p + a - 4) * log (1 - vu);
  vals -= 0.5 * a * log (1 - vu * (1 - 1 / n));
  vals -= 0.5 * (n - 1) * log (1 - vu * R2);

  return (vals);
}

double log_BF_g_on_n_quad (const int n, const int p, const double R2, const int a)
{
  auto f=[=](double x) {
    return log_BF_g_on_n_integrand (x, n, p, R2, a);
  };
  Rosetta::GaussLegendreQuadrature < 1000 > gauss_legendre;
  return log(gauss_legendre.integrate (0., 1., f));
}

void set_log_prob(const string prior, log_prob_fn& log_prob)
{
  if (prior == "maruyama") {
    log_prob = maruyama;
  } else if (prior == "BIC") {
    log_prob = BIC;
  } else if (prior == "ZE") {
    log_prob = ZE;
  } else if (prior == "liang_g1") {
    log_prob = liang_g1;
  } else if (prior == "liang_g2") {
    log_prob = liang_g2;
  } else if (prior == "liang_g_n_appell") {
    log_prob = liang_g_n_appell;
  } else if (prior == "liang_g_n_approx") {
    log_prob = liang_g_n_approx;
  } else if (prior == "liang_g_n_quad") {
    log_prob = liang_g_n_quad;
  } else if (prior == "robust_bayarri1") {
    log_prob = robust_bayarri1;
  } else if (prior == "robust_bayarri2") {
    log_prob = robust_bayarri2;
  } else if (prior == "hyper_g_n_gauss_legendre") {
    log_prob = log_BF_g_on_n_quad;
  } else {
    stringstream ss;
    ss << "Prior " << prior << " unknown";
    Rcpp::stop(ss.str());
  }
}
