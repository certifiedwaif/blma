#include <Rcpp.h>
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include <sstream>
#include <unordered_map>
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS
#include <boost/functional/hash.hpp>

// [[Rcpp::depends(RcppGSL)]]
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>


#include "graycode.h"
#include "correlation.h"

// #define DEBUG

using namespace std;
using namespace Rcpp;


namespace boost
{
	template <typename B, typename A>
	std::size_t hash_value(const boost::dynamic_bitset<B, A>& bs) {
		return boost::hash_value(bs.m_bits);
	}
}


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
	// 	log_p = -INFINITY;
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
		// Rcpp::Rcout << "0.5 * (n - 1.) * log(beta + x[i])	" << 0.5 * (n - 1.) * log(beta + x[i])	<< std::endl;
		// // if (i < 10)
		// // 	Rcpp::Rcout << "x(" << i << ") " << x(i) << " log_f(i) " << log_f(i) << std::endl;
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


double calculate_log_prob(const uint n, const uint p, const double R2, const uint p_gamma,
													const dbitset gamma,
													const std::function<double (const int n, const int p, double vR2, int vp_gamma)> log_prob,
													const std::string modelprior, const VectorXd modelpriorvec)
{
	double result = log_prob(n, p, R2, p_gamma);

	if (modelprior == "beta-binomial") {
		double alpha = modelpriorvec(0);
		double beta = modelpriorvec(1);
		result += ::Rf_lbeta(alpha + p_gamma, beta + p - p_gamma) - ::Rf_lbeta(alpha, beta);
	}

	if (modelprior == "bernoulli") {
		for (auto j = 0; j < p; j++) {
			auto igamma = gamma[j] ? 1. : 0.;
			if (modelpriorvec(j) == 0. || modelpriorvec(j) == 1.)
				continue;
			result += igamma * log(modelpriorvec(j)) + (1 - igamma) * log(1. - modelpriorvec(j));
		}
	}

	return result;
}


void calculate_log_probabilities(const vector< dbitset >& gamma, const VectorXd& sigma2, const int n,
																	VectorXd& log_probs,
																	const std::function<double (const int n, const int p, double vR2, int vp_gamma)> log_prob,
																	const std::string modelprior, const VectorXd modelpriorvec)
{
	auto K = gamma.size();
	auto p = gamma[0].size();
	auto a = 1.;

	for (auto k = 0; k < K; k++) {
		auto p_gamma = gamma[k].count();
		auto b = p;
		log_probs(k) = log_prob(n, p, 1. - sigma2[k], p_gamma);
		#ifdef DEBUG
		// Rcpp::Rcout << "log_probs[" << k << "] " << log_probs[k] << std::endl;
		#endif
		if (modelprior == "beta-binomial") {
			double alpha = modelpriorvec(0);
			double beta = modelpriorvec(1);
			log_probs(k) += ::Rf_lbeta(alpha + p_gamma, beta + p - p_gamma) - ::Rf_lbeta(alpha, beta);
		}
		if (modelprior == "bernoulli") {
			for (auto j = 0; j < p; j++) {
				auto igamma = gamma[k][j] ? 1. : 0.;
				if (modelpriorvec(j) == 0. || modelpriorvec(j) == 1.)
					continue;
				log_probs(k) += igamma * log(modelpriorvec(j)) + (1 - igamma) * log(1. - modelpriorvec(j));
			}
		}
	}
	// Re-normalise
	log_probs = log_probs.array() - log_probs.maxCoeff();
}


void calculate_weights(const VectorXd& sigma2, const VectorXd& log_probs, VectorXd& w)
{
	const auto K = log_probs.size();
	for (auto k = 0; k < K; k++) {
		// Need to use log-sum-exp trick here
		// Have to skip 0 probabilities somehow.
		if (sigma2[k] == 0.) {
			w[k] = 0.;
		} else {
			w[k] = exp(log_probs[k]) / log_probs.array().exp().sum();
		}
		#ifdef DEBUG
		Rcpp::Rcout << "w[" << k << "] " << w[k] << std::endl;
		#endif
	}
}


double calculate_entropy(const VectorXd& w)
{
	const auto K = w.size();
	auto H = 0.;
	for (auto k = 0; k < K; k++) {
		// Don't add 0 weights
		if (w[k] == 0.) {
			continue;
		}	else {
			H += -w[k] * log(w[k]);
		}
	}
	return H;
}


double calculate_w_dot_prob(const VectorXd& w, const VectorXd& log_probs)
{
	const auto K = w.size();
	auto w_dot_prob = 0.;
	for (auto k = 0; k < K; k++) {
		if (w[k] == 0.) {
			continue;
		} else {
			w_dot_prob += w[k] * exp(log_probs[k]);
		}
		#ifdef DEBUG
		Rcpp::Rcout << "w[k] " << w[k] << " log_probs[k] " << log_probs[k] << " w_dot_prob " << k << " " << w_dot_prob << " " << std::endl;
		#endif
	}
	return w_dot_prob;
}


void gamma_to_NumericMatrix(const vector< dbitset >& gamma, NumericMatrix& nm)
{
	auto K = gamma.size();
	auto p = gamma[0].size();
	for (auto k = 0; k < K; k++) {
		for (auto j = 0; j < p; j++) {
			nm(k, j) = gamma[k][j] ? 1. : 0.;
		}
	}
}


//' Run a Collapsed Variational Approximation to find the K best linear models
//'
//' @param vy_in Vector of responses
//' @param mX_in The matrix of covariates which may or may not be included in each model
//' @param mGamma_in Matrix of initial models, a K by p logical matrix
//' @param prior -- the choice of mixture $g$-prior used to perform Bayesian model averaging. The choices
//' available include:
//' 	\itemize{
//' 		\item{"BIC"}{-- the Bayesian information criterion obtained by using the cake prior
//' 		of Ormerod et al. (2017).}
//'
//' 		\item{"ZE"}{-- special case of the prior structure described by Maruyama and George (2011).}
//'
//' 		\item{"liang_g1"}{-- the mixture \eqn{g}-prior of Liang et al. (2008) with prior hyperparameter
//'     \eqn{a=3} evaluated directly using Equation (10) of Greenaway and Ormerod (2018) where the Gaussian
//'			hypergeometric function is evaluated using the {gsl} library. Note: this option can lead to numerical problems and is only
//''    meant to be used for comparative purposes.}
//'
//' 		\item{"liang_g2"}{-- the mixture \eqn{g}-prior of Liang et al. (2008) with prior hyperparameter
//' 		 \eqn{a=3} evaluated directly using Equation (11)  of Greenaway and Ormerod (2018).}
//'
//' 		\item{"liang_g_n_appell"}{-- the mixture \eqn{g/n}-prior of Liang et al. (2008) with prior
//'			 hyperparameter \eqn{a=3} evaluated using the {appell R} package.}
//'
//' 		\item{"liang_g_approx"}{-- the mixture \eqn{g/n}-prior of Liang et al. (2008) with prior hyperparameter
//'      \eqn{a=3} using the approximation Equation (15)  of Greenaway and Ormerod (2018) for model with more
//' 		  than two covariates and numerical quadrature (see below) for models with one or two covariates.}
//'
//' 		\item{"liang_g_n_quad"}{-- the mixture \eqn{g/n}-prior of Liang et al. (2008) with prior hyperparameter
//'			 \eqn{a=3} evaluated using a composite trapezoid rule.}
//'
//' 		\item{"robust_bayarri1"}{-- the robust prior of Bayarri et al. (2012) using default prior hyper
//'			parameter choices evaluated directly using Equation (18)  of Greenaway and Ormerod (2018) with the
//'     {gsl} library.}
//'
//' 		\item{"robust_bayarri2"}{-- the robust prior of Bayarri et al. (2012) using default prior hyper
//'			parameter choices evaluated directly using Equation (19) of Greenaway and Ormerod (2018).}
//' }
//' @param modelprior The model prior to use. The choices of model prior are "uniform", "beta-binomial" or
//' "bernoulli". The choice of model prior dictates the meaning of the parameter modelpriorvec.
//' @param modelpriorvec If modelprior is "uniform", then the modelpriorvec is ignored and can be null.
//'
//' If
//' the modelprior is "beta-binomial" then modelpriorvec should be length 2 with the first element containing
//' the alpha hyperparameter for the beta prior and the second element containing the beta hyperparameter for
//' the beta prior.
//'
//' If modelprior is "bernoulli", then modelpriorvec must be of the same length as the number
//' of columns in mX. Each element i of modelpriorvec contains the prior probability of the the ith covariate
//' being included in the model.
//' @param bUnique Whether to ensure uniqueness in the population of particles or not. Defaults to true.
//' @param lambda The weighting factor for the entropy in f_lambda. Defaults to 1.
//' @return A list containing the named element models, which is a K by p matrix of the models
//'					selected by the algorithm, and the named element trajectory, which includes a list
//'					of the populations of models for each iteration of the algorithm until it converged
//' @examples
//' mD <- MASS::UScrime
//' notlog <- c(2,ncol(MASS::UScrime))
//' mD[,-notlog] <- log(mD[,-notlog])
//'
//' for (j in 1:ncol(mD)) {
//'   mD[,j] <- (mD[,j] - mean(mD[,j]))/sd(mD[,j])
//' }
//'
//' varnames <- c(
//'   "log(AGE)",
//'   "S",
//'   "log(ED)",
//'   "log(Ex0)",
//'   "log(Ex1)",
//'   "log(LF)",
//'   "log(M)",
//'   "log(N)",
//'   "log(NW)",
//'   "log(U1)",
//'   "log(U2)",
//'   "log(W)",
//'   "log(X)",
//'   "log(prison)",
//'   "log(time)")
//'
//' y.t <- mD$y
//' X.f <- data.matrix(cbind(mD[1:15]))
//' colnames(X.f) <- varnames
//' K <- 100
//' p <- ncol(X.f)
//' initial_gamma <- matrix(rbinom(K * p, 1, .5), K, p)
//' cva_result <- cva(y.t, X.f, initial_gamma, prior = "BIC", modelprior = "uniform",
//'                   modelpriorvec=c(0.))
//' @references
//' Bayarri, M. J., Berger, J. O., Forte, A., Garcia-Donato, G., 2012. Criteria for Bayesian
//' model choice with application to variable selection. Annals of Statistics 40 (3), 1550-
//' 1577.
//'
//' Greenaway, M. J., J. T. Ormerod (2018) Numerical aspects of Bayesian linear models averaging using mixture
//' g-priors.
//'
//' Liang, F., Paulo, R., Molina, G., Clyde, M. a., Berger, J. O., 2008. Mixtures of g priors for
//' Bayesian variable selection. Journal of the American Statistical Association 103 (481),
//' 410-423.
//'
//' Ormerod, J. T., Stewart, M., Yu, W., Romanes, S. E., 2017. Bayesian hypothesis tests
//' with diffuse priors: Can we have our cake and eat it too?
//' @export
// [[Rcpp::export]]
List cva(const NumericVector vy_in, const NumericMatrix mX_in,
				 const NumericMatrix mGamma_in, 
				 const std::string prior,
				 const std::string modelprior, const NumericVector modelpriorvec_in,
				 const bool bUnique = true,
				 const double lambda = 1.)
{
	VectorXd vy(vy_in.length());   // = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(vy_in);
	for (auto i = 0; i < vy_in.length(); i++) vy[i] = vy_in[i];
																 // = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(mX_in);
	MatrixXd mX(mX_in.nrow(), mX_in.ncol());
	for (auto i = 0; i < mX_in.nrow(); i++)
		for (auto j = 0; j < mX_in.ncol(); j++)
			mX(i, j) = mX_in(i, j);
	VectorXd modelpriorvec(modelpriorvec_in.length());
	for (auto i = 0; i < modelpriorvec_in.length(); i++)
		modelpriorvec(i) = modelpriorvec_in(i);
	MatrixXd mGamma(mGamma_in.nrow(), mGamma_in.ncol());
	for (auto i = 0; i < mGamma_in.nrow(); i++)
		for (auto j = 0; j < mGamma_in.ncol(); j++)
			mGamma(i, j) = mGamma_in(i, j);
	const auto n = mX.rows();
	const auto p = mX.cols();
	const MatrixXd mXTX = mX.transpose() * mX;
	const MatrixXd mXTy = mX.transpose() * vy;
	const uint K = mGamma.rows();
	vector< dbitset > gamma(K);
	VectorXd log_probs(K);
	VectorXd w(K);
	vector< vector< dbitset > > trajectory;
	vector< VectorXd > trajectory_probs;
	vector< MatrixXd > mXTX_inv(K);
	VectorXd sigma2(K);
	std::unordered_map< std::size_t, bool > hash;
	// const gsl_rng_type *T;
	// gsl_rng *r;
	const auto RHO = 0.1;
	auto f_lambda_prev = 0.;
	const auto EPSILON = 1e-8;

	std::function<double (const int n, const int p, double vR2, int vp_gamma)> log_prob;
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
	} else {
		stringstream ss;
		ss << "Prior " << prior << " unknown";
		Rcpp::stop(ss.str());
	}
	// Initialise population of K particles randomly
	// Rcpp::Rcout << "Generated" << std::endl;
	// for (auto k = 0; k < K; k++) {
	// 	gamma[k].resize(p);
	// 	// Empty models are pointless, so we keep regenerating gammas until we get a non-zero one
	// 	while (gamma[k].count() == 0) {
	// 		for (auto j = 0; j < p; j++) {
	// 			if (gsl_rng_uniform(r) <= RHO)
	// 				gamma[k][j] = true;
	// 			else
	// 				gamma[k][j] = false;
	// 		}
	// 	}
	// 	Rcpp::Rcout << "gamma[" << k << "] " << gamma[k] << std::endl;
	// }

	if (mGamma.rows() != K || mGamma.cols() != p) {
		stringstream ss;
		ss << "initial_gamma is " << to_string(mGamma.rows()) << " by " << to_string(mGamma.cols()) << ", expected " << K << " by " << p;
		stop(ss.str());
	}

	#ifdef DEBUG
	Rcpp::Rcout << "initial_gamma" << std::endl;
	#endif
	for (auto k = 0; k < K; k++) {
		gamma[k].resize(p);
		for (auto j = 0; j < p; j++) {
			if (mGamma(k, j) == 1.) {
				gamma[k][j] = true;
			}
			else {
				gamma[k][j] = false;
			}
		}
		#ifdef DEBUG
		Rcpp::Rcout << "gamma[" << k << "] " << gamma[k] << std::endl;
		#endif
		if (bUnique) {
			hash.insert({boost::hash_value(gamma[k]), true});
		}
	}

	// Initialise mXTX_inv
	// Initialise sigma2
	for (auto k = 0; k < K; k++) {
		auto p_gamma = gamma[k].count();
		if (p_gamma == 0) {
			stringstream ss;
			ss << "gamma[" << k + 1 << "] has no bits set" << endl;
			throw domain_error(ss.str());
		}
		MatrixXd mX_gamma(n, p_gamma);
		get_cols(mX, gamma[k], mX_gamma);
		MatrixXd mX_gamma_Ty(p_gamma, 1);
		get_rows(mXTy, gamma[k], mX_gamma_Ty);
		mXTX_inv[k] = (mX_gamma.transpose() * mX_gamma).inverse();
		sigma2[k] = 1. - (mX_gamma_Ty.transpose() * mXTX_inv[k] * mX_gamma_Ty).value() / n;
		#ifdef DEBUG
		Rcpp::Rcout << "sigma2[" << k << "] " << sigma2[k] << std::endl;
		#endif
	}
	calculate_log_probabilities(gamma, sigma2, n, log_probs, log_prob, modelprior, modelpriorvec);
	trajectory.push_back(gamma);
	trajectory_probs.push_back(log_probs.array().exp());

	// Loop until convergence
	auto converged = false;
	auto iteration = 1;
	while (!converged) {
		#ifdef DEBUG
		Rcpp::Rcout << "Iteration " << iteration << std::endl;
		#endif

		// #pragma omp parallel for
		for (auto k = 0; k < K; k++) {
			#ifdef DEBUG
			Rcpp::Rcout << "gamma[" << k << "] " << gamma[k] << std::endl;
			#endif
			for (auto j = 0; j < p; j++) {
				dbitset gamma_prime = gamma[k];
				gamma_prime[j] = !gamma_prime[j];
				auto p_gamma = gamma[k].count();
				auto p_gamma_prime = gamma_prime.count();
				if ((p_gamma_prime == 0) || (p_gamma_prime >= n - 1))
					continue;
				bool bUpdate = !gamma[k][j];

				// If we've seen this bitstring before, don't do the update
				if (bUnique) {
					auto h = boost::hash_value(gamma_prime);
					auto search = hash.find(h);
					if (search != hash.end()) {
						continue;
					}	else {
						hash.insert({h, true});
					}
				}
				#ifdef DEBUG
				if (bUpdate)
					Rcpp::Rcout << "Updating " << j << std::endl;
				else
					Rcpp::Rcout << "Downdating " << j << std::endl;
				#endif

				// Update or downdate mXTX_inv
				MatrixXd mXTX_inv_prime;
				mXTX_inv_prime.resize(p_gamma_prime, p_gamma_prime);
				bool bLow;             // gamma_1 or gamma[k]?
				uint min_idx = std::min(gamma[k].find_first(), gamma_prime.find_first());
				#ifdef DEBUG
				Rcpp::Rcout << "min_idx " << min_idx << std::endl;
				#endif

				// This is totally evil
				// Explanation: mXTX is addressed absolutely, but mXTX_inv is addressed relatively. To account for
				// this, we abuse fixed to adjust for the gaps in gamma_prime
				int fixed = 0;
				for (auto idx = min_idx; idx < j; idx++) {
					if (!gamma_prime[idx])
						fixed--;
				}

				if (bUpdate) {
					rank_one_update(gamma[k], j, min_idx,
						fixed,
						mXTX,
						mXTX_inv[k],
						mXTX_inv_prime,
						bLow);
						#ifdef DEBUG
						Rcpp::Rcout << "bLow " << bLow << std::endl;
						#endif
				} else {
					rank_one_downdate(j, min_idx, fixed,
						mXTX_inv[k], mXTX_inv_prime);
				}

				// Calculate sigma2_prime
				double sigma2_prime;

				MatrixXd mX_gamma_prime(n, p_gamma_prime);
				get_cols(mX, gamma_prime, mX_gamma_prime);
				MatrixXd mX_gamma_prime_Ty = mX_gamma_prime.transpose() * vy;
				double nR2 = (mX_gamma_prime_Ty.transpose() * mXTX_inv_prime * mX_gamma_prime_Ty).value();
				sigma2_prime = 1. - nR2 / n;
				// #ifdef DEBUG
				// It's mathematically impossible that nR2 can be that high. Therefore, our inverse must be bad.
				// Recalculate it from scratch and try again.

				// Rcpp::Rcout << "gamma[" << k << "]    " << gamma[k] << std::endl;
				// Rcpp::Rcout << "gamma_prime " << gamma_prime << std::endl;
				#ifdef DEBUG
				MatrixXd mXTX_inv_prime_check = (mX_gamma_prime.transpose() * mX_gamma_prime).inverse();
				if (!mXTX_inv_prime.isApprox(mXTX_inv_prime_check, 1e-8)) {
					Rcpp::Rcout << "gamma[k]    " << gamma[k] << std::endl;
					Rcpp::Rcout << "gamma_prime " << gamma_prime << std::endl;
					Rcpp::Rcout << "mXTX_inv_prime " << std::endl << mXTX_inv_prime << std::endl;
					Rcpp::Rcout << "mXTX_inv_prime_check " << std::endl << mXTX_inv_prime_check << std::endl;
					Rcpp::Rcout << "Difference " << std::endl << mXTX_inv_prime - mXTX_inv_prime_check << std::endl;
					// throw std::range_error("Rank one update failed. I wonder why?");
				}

				double nR2_check = (mX_gamma_prime_Ty.transpose() * mXTX_inv_prime * mX_gamma_prime_Ty).value();
				double sigma2_prime_check = 1. - nR2_check / n;
				if (abs(sigma2_prime - sigma2_prime_check) > EPSILON) {
					Rcpp::Rcout << "sigma2_prime " << sigma2_prime << " sigma2_prime_check " << sigma2_prime_check << std::endl;
				}
				mXTX_inv_prime = mXTX_inv_prime_check;
				#endif
				// Rcpp::Rcout << "sigma2_1 " << sigma2_1 << std::endl;
				// #endif

				#ifdef DEBUG
				Rcpp::Rcout << "sigma2[k] " << sigma2[k] << std::endl;
				Rcpp::Rcout << "sigma2_prime " << sigma2_prime << std::endl;
				#endif
				double log_p_0;
				double log_p_1;
				if (bUpdate) {
					// log_p_0 = log_prob(n, p, 1. - sigma2[k], p_gamma);
					// log_p_1 = log_prob(n, p, 1. - sigma2_prime, p_gamma_prime);
					log_p_0 = calculate_log_prob(n, p, 1. - sigma2[k], p_gamma, gamma[k], log_prob, modelprior, modelpriorvec);
					log_p_1 = calculate_log_prob(n, p, 1. - sigma2_prime, p_gamma_prime, gamma_prime, log_prob, modelprior, modelpriorvec);
				} else {
					// log_p_0 = log_prob(n, p, 1. - sigma2_prime, p_gamma_prime);
					// log_p_1 = log_prob(n, p, 1. - sigma2[k], p_gamma);
					log_p_0 = calculate_log_prob(n, p, 1. - sigma2_prime, p_gamma_prime, gamma_prime, log_prob, modelprior, modelpriorvec);
					log_p_1 = calculate_log_prob(n, p, 1. - sigma2[k], p_gamma, gamma[k], log_prob, modelprior, modelpriorvec);
				}
				#ifdef DEBUG
				Rcpp::Rcout << "log_p_0 " << log_p_0;
				Rcpp::Rcout << " log_p_1 " << log_p_1 << std::endl;
				#endif
				if ((log_p_0 > log_p_1 && !bUpdate) || (log_p_1 > log_p_0 && bUpdate)) {
					if (bUnique) {
						hash.erase(boost::hash_value(gamma[k]));
					}
					gamma[k][j] = bUpdate;
					if (bUnique) {
						hash.insert({boost::hash_value(gamma[k]), true});
					}
					#ifdef DEBUG
					if (bUpdate)
						Rcpp::Rcout << "Keep update" << std::endl;
					else
						Rcpp::Rcout << "Keep downdate" << std::endl;
					#endif
					sigma2[k] = sigma2_prime;
					mXTX_inv[k] = mXTX_inv_prime;
				} else {
					#ifdef DEBUG
					if (bUpdate)
						Rcpp::Rcout << "Don't keep update" << std::endl;
					else
						Rcpp::Rcout << "Don't keep downdate" << std::endl;
					#endif
				}
			}
		}

		calculate_log_probabilities(gamma, sigma2, n, log_probs, log_prob, modelprior, modelpriorvec);

		// Calculate weights
		calculate_weights(sigma2, log_probs, w);

		// Calculate entropy
		auto H = calculate_entropy(w);

		// Rcpp::Rcout << "w.dot(probs) " << w.dot(probs) << std::endl;
		auto w_dot_prob = calculate_w_dot_prob(w, log_probs);
		#ifdef  DEBUG
		Rcpp::Rcout << "w_dot_prob " << w_dot_prob << std::endl;
		Rcpp::Rcout << "H " << H << std::endl;
		#endif

		// Calculate f_lambda
		double f_lambda = w_dot_prob + lambda * H;
		#ifdef DEBUG
		Rcpp::Rcout << "f_lambda_prev " << f_lambda_prev << " f_lambda " << f_lambda << std::endl;
		#endif

		// Check for convergence - is f_lambda changed from the last iteration?
		if ((f_lambda - f_lambda_prev) > EPSILON) {
			f_lambda_prev = f_lambda;
		}
		else {
			converged = true;
		}
		for (auto k = 0; k < K; k++) {
			#ifdef DEBUG
			Rcpp::Rcout << "gamma[" << k + 1 << "] " << gamma[k] << std::endl;
			#endif
		}
		iteration++;
		trajectory.push_back(gamma);
		trajectory_probs.push_back(log_probs.array().exp());
	}
	#ifdef DEBUG
	Rcpp::Rcout << "Converged" << std::endl;
	#endif
	NumericMatrix bitstrings(K, p);
	gamma_to_NumericMatrix(gamma, bitstrings);

	List trajectory_bitstrings;
	NumericMatrix trajectory_probabilities(K, trajectory.size());
	for (auto i = 0; i < trajectory.size(); i++) {
		NumericMatrix bitstrings2(K, p);
		gamma_to_NumericMatrix(trajectory[i], bitstrings2);
		trajectory_bitstrings.push_back(bitstrings2);
		for (auto k = 0; k < K; k++) {
			trajectory_probabilities(k, i) = trajectory_probs[i](k);
		}
	}

	List result = List::create(Named("models") = bitstrings,
															Named("trajectory") = trajectory_bitstrings,
															Named("trajectory_probs") = trajectory_probabilities);
	return result;
}
