// priors.cpp

#include <Rcpp.h>
#include <cmath>
#include <ccomplex>
#include <inttypes.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#define EIGEN_USE_BLAS

// [[Rcpp::depends(RcppGSL)]]
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>

#include "GaussLegendre.h"
#include "laguerre_rule.h"
#include "priors.h"

using namespace Rcpp;

using Eigen::VectorXd;
using std::string;


//' BIC prior
//'
//' @param n The sample size, an integer
//' @param p_gamma The number of covariates in the model gamma
//' @param R2 The correlation co-efficient, a number between -1 and 1
//' @return The log of the Bayes Factor
//' @export
// [[Rcpp::export]]
double BIC(const int n, const int p_gamma, double R2)
{
  	const auto sigma2 = 1. - R2;
#ifdef DEBUG
  	Rcpp::Rcout << "n " << n << " p " << p << " R2 " << R2 << " p_gamma " << p_gamma << std::endl;
#endif
  	auto BIC = n * log(sigma2) + p_gamma * log(n);
  	if (sigma2 == 0. || std::isnan(sigma2)) {
    	BIC = -INFINITY;
    	// throw std::range_error("Invalid sigma2");
  	}
#ifdef DEBUG
  	Rcpp::Rcout << "n * log(1 - R2) " << n * log(1 - R2) << std::endl;
  	Rcpp::Rcout << "p_gamma * log(n) " << p_gamma * log(n) << std::endl;
#endif
  	return -0.5 * BIC;
}


//' ZE prior
//'
//' @param n The sample size, an integer
//' @param p_gamma The number of covariates in the model gamma
//' @param R2 The correlation co-efficient, a number between -1 and 1
//' @return The log of the Bayes Factor
//' @export
// [[Rcpp::export]]
double ZE(const int n, const int p_gamma, const double R2)
{
  	auto a = -0.75;
  	auto b = 0.5 * (n - p_gamma - 5) - a;
  	auto d = 0.5 * p_gamma + a;

  	auto log_p = -(b+1)*log(1 - R2) + Rf_lbeta(d+1,b+1) - Rf_lbeta(a+1,b+1);
  	return log_p;
}


//' log_hyperg_2F1 prior
//'
//' @param b
//' @param c
//' @param x
//' @return
//' @export
// [[Rcpp::export]]
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


//' log_hyperg_2F1_naive
//'
//' @param b
//' @param c
//' @param x
//' @return The log of the Bayes Factor
//' @export
// [[Rcpp::export]]
double log_hyperg_2F1_naive(double b, double c, double x)
{
  	auto val = log(gsl_sf_hyperg_2F1( b, 1, c, x));
  	return val;
}


//' Liang's hyper g-prior
//'
//' @param n The sample size, an integer
//' @param p_gamma The number of covariates in the model gamma
//' @param R2 The correlation co-efficient, a number between -1 and 1
//' @return The log of the Bayes Factor
//' @export
// [[Rcpp::export]]
double liang_g1(const int n, const int p_gamma, const double R2)
{
  	auto a = 3.;
  	double log_p_g;
  	log_p_g = log(a - 2) - log(p_gamma + a - 2) + log(gsl_sf_hyperg_2F1(0.5*(n-1), 1, 0.5*(p_gamma + a), R2));
  	return log_p_g;
}


//' Liang's g prior
//'
//' @param n The sample size, an integer
//' @param p_gamma The number of covariates in the model gamma
//' @param R2 The correlation co-efficient, a number between -1 and 1
//' @return The log of the Bayes Factor
//' @export
// [[Rcpp::export]]
double liang_g2(const int n, const int p_gamma, const double R2)
{
  	auto a = 3.;
  	auto log_vp_g2 = log(a - 2) - log(p_gamma + a - 2) + log_hyperg_2F1( 0.5*(n-1), 0.5*(p_gamma+a), R2);
  	return log_vp_g2;
}


extern "C" void f1(std::complex<double>* a, std::complex<double>* b1,
		std::complex<double>* b2, std::complex<double>* c, double* x, double*
		y, int32_t* algoflag, int32_t* userflag, bool* debug,
		std::complex<double>* val, int* hyp2f1);


//' Liang's g/n prior Appell
//'
//' @param n The sample size, an integer
//' @param p_gamma The number of covariates in the model gamma
//' @param R2 The correlation co-efficient, a number between -1 and 1
//' @return The log of the Bayes Factor
//' @export
// [[Rcpp::export]]
double liang_g_n_appell(const int n, const int p_gamma, const double R2)
{
	double result = 1;
	if (p_gamma == 0)
		return result;
#ifdef DEBUG
	Rcpp::Rcout << "liang_g_n_appell(" << n << ", " << p << ", " << R2 << ", " << p_gamma << ") = ";
#endif


	double a_prime_dub = 3.;

	std::complex<double> a_prime = 3.;

	std::complex<double> val;
	std::complex<double> a = 1.;
	std::complex<double> b1 = a_prime.real() / 2.;
	std::complex<double> b2 = (n - 1.)/2.;
	std::complex<double> c = (p_gamma + a_prime.real()) / 2.;

	double x = 1. - 1. / n;
	double y = R2;
	int32_t algoflag = 0;
	int32_t userflag = 1;
	int hyp2f1 = 2; // "michel.stoitsov";
	bool debug = false;
#pragma omp master
	{
	f1(&a, &b1, &b2, &c, &x, &y, &algoflag, &userflag, &debug, &val, &hyp2f1);
	}
	result = log(a_prime_dub - 2.) - log(n) - log(p_gamma + a_prime_dub -
			2.) + log(val.real());
#ifdef DEBUG
	Rcpp::Rcout << result << std::endl;
#endif
	return result;
}


// Trapezoidal integration over a potentially irregular grid
double trapint(const VectorXd& xgrid, const VectorXd& fgrid)
{
  	auto sum = 0.;

#pragma omp parallel for simd reduction(+:sum)
  	for (auto i = 0; i < xgrid.size() - 1; i++) {
    	sum += 0.5 * (xgrid(i + 1) - xgrid(i)) * (fgrid(i) + fgrid(i + 1));
  	}
  	// Rcpp::Rcout << "sum " << sum << std::endl;

  	return sum;
}


double liang_g_n_quad_integrand(const int n, const int p_gamma, const double R2,
								const double u)
{
  	const auto a = 3.;
	return exp((p_gamma / 2. + a / 2. - 2.) * log(1 - u) + -a/2. * log(1. - u * (1. - 1. / n)) + (-(n-1.)/2.) * log(1 - u*R2));
}

//' Liang's g/n prior quadrature
//'
//' @param n The sample size, an integer
//' @param p_gamma The number of covariates in the model gamma
//' @param R2 The correlation co-efficient, a number between -1 and 1
//' @return The log of the Bayes Factor
//' @export
// [[Rcpp::export]]
double liang_g_n_quad(const int n, const int p_gamma, const double R2)
{
  	const auto a = 3.;
  	const int NUM_POINTS = 10000;
  	auto sum = 0.;
#pragma omp parallel for simd\
	reduction(+:sum)\
	default(none)
  	for (int i = 1; i < NUM_POINTS; i++) {
    	double u_prev = static_cast<double>(i - 1) / static_cast<double>(NUM_POINTS);
    	double u = static_cast<double>(i) / static_cast<double>(NUM_POINTS);
    	auto fgrid = liang_g_n_quad_integrand(n, p_gamma, R2, u_prev);
    	auto fgrid_prev = liang_g_n_quad_integrand(n, p_gamma, R2, u);
    	sum += 0.5 * (fgrid_prev + fgrid) / static_cast<double>(NUM_POINTS);
  	}

  	auto result = log(a - 2.) - log(2. * n) + log(sum);
#ifdef DEBUG
  	Rcpp::Rcout << "liang_g_n_quad(" << n << ", " << p << ", " << R2 << ", " << p_gamma << ") = " << result << std::endl;
#endif
  	return result;
}


//' Liang's g/n prior approximation
//' Liang's g/n prior quadrature
//'
//' @param n The sample size, an integer
//' @param p_gamma The number of covariates in the model gamma
//' @param R2 The correlation co-efficient, a number between -1 and 1
//' @return The log of the Bayes Factor
//' @export
// [[Rcpp::export]]
double liang_g_n_approx(const int n, const int p_gamma, const double R2)
{
  	// #ifdef DEBUG
  	// Rcpp::Rcout << "n " << n << " p " << p << " R2 " << R2 << " p_gamma " << p_gamma << std::endl;
  	// #endif
  	if (p_gamma == 0)
    	return 0.;
  	if (p_gamma == 1 || p_gamma == 2)
    	return liang_g_n_quad(n, p_gamma, R2);

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


//' Robust Bayarri 1
//'
//' @param n The sample size, an integer
//' @param p_gamma The number of covariates in the model gamma
//' @param R2 The correlation co-efficient, a number between -1 and 1
//' @return The log of the Bayes Factor
//' @export
// [[Rcpp::export]]
double robust_bayarri1(const int n, const int p_gamma, const double R2)
{
  	// Rcpp::Rcout << "n " << n << " R2 " << R2 << " p_gamma " << p_gamma << std::endl;
  	double r = (1. + n) / (1. + p_gamma);
  	double L = r - 1.;

  	const int NUM_POINTS = 10000;
  	VectorXd x(NUM_POINTS);
  	x.setLinSpaced(NUM_POINTS, L, 10000);

  	VectorXd log_f(NUM_POINTS);
#pragma omp parallel for\
	shared(log_f, x, r)\
	default(none)
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


//' Robust Bayarri 2
//'
//' @param n The sample size, an integer
//' @param p_gamma The number of covariates in the model gamma
//' @param R2 The correlation co-efficient, a number between -1 and 1
//' @return The log of the Bayes Factor
//' @export
// [[Rcpp::export]]
double robust_bayarri2(const int n, const int p_gamma, const double R2)
{
#ifdef DEBUG
  	Rcpp::Rcout << "n " << n;
  	Rcpp::Rcout << " p " << p;
  	Rcpp::Rcout << " R2 " << R2;
  	Rcpp::Rcout << " p_gamma " << p_gamma;
#endif
  	auto sigma2 = 1. - R2;
  	auto L = (1. + n)/(1. + p_gamma) - 1.;

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


//' log_BF_g_on_n_integrand
//'
//' @param vu The argument, vu
//' @param n The sample size, an integer
//' @param p_gamma The number of covariates in the model gamma
//' @param R2 The correlation co-efficient, a number between -1 and 1
//' @return The log of the Bayes Factor
//' @export
// [[Rcpp::export]]
double log_BF_g_on_n_integrand (const double vu, const int n, const int p_gamma, const double R2)
{
  	auto a = 3.;
  	auto result = 0.;
  	result += log (a - 2.);
  	result -= log (2. * n);
  	result += 0.5 * (p_gamma + a - 4.) * log (1. - vu);
  	result -= 0.5 * a * log (1. - vu * (1. - 1. / n));
  	result -= 0.5 * (n - 1.) * log (1. - vu * R2);

  	return result;
}

//' hyper-g/n Gauss-Legendre quadrature
//'
//' @param n The sample size, an integer
//' @param p_gamma The number of covariates in the model gamma
//' @param R2 The correlation co-efficient, a number between -1 and 1
//' @return The log of the Bayes Factor
//' @export
// [[Rcpp::export]]
double log_BF_g_on_n_quad (const int n, const int p_gamma, const double R2)
{
  	auto f = [=](double x) {
    	return log_BF_g_on_n_integrand(x, n, p_gamma, R2);
  	};
  	static Rosetta::GaussLegendreQuadrature < 1000 > gauss_legendre;
  	return gauss_legendre.integrate(0., 1., f);
}


//' log_BF_Zellner_Siow_integrand
//'
//' @param x The argument, x
//' @param n The sample size, an integer
//' @param p_gamma The number of covariates in the model gamma
//' @param R2 The correlation co-efficient, a number between -1 and 1
//' @return The log of the Bayes Factor
//' @export
// [[Rcpp::export]]
double log_BF_Zellner_Siow_integrand(double x, const int n, const int p_gamma, const double R2)
{
  	auto sigma2 = 1. - R2;
  	auto vz = 0.5 * n * sigma2;
  	auto result = 0.5 * (p_gamma - 1.) * log(2. * x / n);
  	result += 0.5 * (n - p_gamma - 1.) * log(1. + 2. * x / n);
  	result -= 0.5 * (n - 1.) * log(1. + (x / vz));
  	result -= 0.5 * (n - 1) * log(sigma2);

  	return result;
}

const int order = 1000;
bool calculated = false;
double w[order];
double x[order];

//' Zellner-Siow Gauss-Laguerre quadrature
//'
//' @param n The sample size, an integer
//' @param p_gamma The number of covariates in the model gamma
//' @param R2 The correlation co-efficient, a number between -1 and 1
//' @return The log of the Bayes Factor
//' @export
// [[Rcpp::export]]
double log_BF_Zellner_Siow_quad(const int n, const int p_gamma, const double R2)
{
  	const double alpha = 0.;
  	const double beta = 0.;
  	const double a = 0.;
  	const double b = 1.;
  	const int kind = 5;
  	
#pragma omp master
	{
  		if (!calculated) {
			cgqf(order, kind, alpha, beta, a, b, x, w);
			calculated = true;
  		}
  	}
  	
  	//for (auto i = 0; i < order; i++) {
  	//	Rcpp::Rcout << "x[" << i << "]=" << x[i] << "w[" << i << "]=" << w[i] << std::endl;
  	//}
  	
  	auto *logf = new double[order];
#pragma omp parallel for
  	for (auto i = 0; i < order; i++) {
  		logf[i] = log_BF_Zellner_Siow_integrand(x[i], n, p_gamma, R2);
  	}

  	double maxval = logf[0];
#pragma omp parallel for \
	reduction(max:maxval)
  	for (auto i = 0; i < order; i++) {
  		if (logf[i] > maxval) {
			maxval = logf[i];
  		}
  	}
  	
  	double sumval = 0;
#pragma omp parallel for \
	reduction(+:sumval)
  	for (auto i = 0; i < order; i++) {
  		sumval += w[i] * exp(logf[i] - maxval);
  	}
  	
  	double pi = 3.14159265358979323846264338327950;
  	double finalval =  maxval + log(sumval) +  0.5*log(0.5*n/pi) - log(0.5*n);
  	
  	delete [] logf;
  	return finalval;
}


void set_log_prob(const string prior, log_prob_fn& log_prob)
{
  	if (prior == "BIC") {
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
  	} else if (prior == "zellner_siow_gauss_laguerre") {
    	log_prob = log_BF_Zellner_Siow_quad;
  	} else {
		std::stringstream ss;
    	ss << "Prior " << prior << " unknown";
    	Rcpp::stop(ss.str());
  	}
}
