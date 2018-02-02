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

#ifdef _OPENMP
	#include <omp.h>
#endif

#include "graycode.h"
#include "correlation.h"
#include "priors.h"

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


double calculate_log_prob(const uint n, const uint p, const double R2, const uint p_gamma,
													const dbitset& gamma,
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
																	const std::string& modelprior, const VectorXd& modelpriorvec)
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


void gamma_to_MatrixXd(const vector< dbitset >& gamma, MatrixXd& m)
{
	auto K = gamma.size();
	auto p = gamma[0].size();
	for (auto k = 0; k < K; k++) {
		for (auto j = 0; j < p; j++) {
			m(k, j) = gamma[k][j] ? 1. : 0.;
		}
	}
}


void calculate_mXTX_inv_prime(const dbitset& gamma, const dbitset& gamma_prime, int j,
															const MatrixXd& mXTX, const MatrixXd& mXTX_inv, MatrixXd& mXTX_inv_prime,
															bool bUpdate)
{
	bool bLow;
	uint min_idx = std::min(gamma.find_first(), gamma_prime.find_first());
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
		rank_one_update(gamma, j, min_idx,
			fixed,
			mXTX,
			mXTX_inv,
			mXTX_inv_prime,
			bLow);
			#ifdef DEBUG
			Rcpp::Rcout << "bLow " << bLow << std::endl;
			#endif
	} else {
		rank_one_downdate(j, min_idx, fixed,
			mXTX_inv, mXTX_inv_prime);
	}
}


double calculate_sigma2_prime(const uint n, const uint p_gamma_prime,
															const MatrixXd& mX, const dbitset& gamma_prime,
															const VectorXd& vy, MatrixXd& mXTX_inv_prime)
{
	MatrixXd mX_gamma_prime(n, p_gamma_prime);
	get_cols(mX, gamma_prime, mX_gamma_prime);
	MatrixXd mX_gamma_prime_Ty = mX_gamma_prime.transpose() * vy;
	double nR2 = (mX_gamma_prime_Ty.transpose() * mXTX_inv_prime * mX_gamma_prime_Ty).value();
	double sigma2_prime = 1. - nR2 / n;
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
	sigma2_prime = sigma2_prime_check;
	#endif
	// Rcpp::Rcout << "sigma2_1 " << sigma2_1 << std::endl;
	// #endif

	return sigma2_prime;
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
//' @param modelpriorvec_in If modelprior is "uniform", then the modelpriorvec is ignored and can be null.
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
//' @param cores The number of cores to use. Defaults to 1.
//' @return The object returned is a list containing:
//' \itemize{
//'		\item{"mGamma"}{-- A K by p binary matrix containing the final population of models}
//'
//'		\item{"vgamma.hat"}{-- The most probable model found by cva}
//'
//'		\item{"vBF"}{-- The null-based Bayes factor for each model in the population}
//'
//'		\item{"posterior_model_probabilities"}{-- The estimated posterior model parameters for each model in
//'		the population.}
//'
//'		\item{"posterior_inclusion_probabilities"}{-- The estimated variable inclusion probabilities for each
//'		model in the population.}
//'
//'		\item{"vR2"}{-- The fitted R-squared values for each model in the population.}
//'
//'		\item{"vp"}{-- The model size for each model in the population.}
//'	}
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
//'                   modelpriorvec_in=NULL)
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
				 const std::string modelprior, const Nullable<NumericVector> modelpriorvec_in = R_NilValue,
				 const bool bUnique = true,
				 const double lambda = 1.,
				 const int cores = 1)
{
	VectorXd vy(vy_in.length());   // = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(vy_in);
	for (auto i = 0; i < vy_in.length(); i++) vy[i] = vy_in[i];
																 // = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(mX_in);
	MatrixXd mX(mX_in.nrow(), mX_in.ncol());
	for (auto i = 0; i < mX_in.nrow(); i++)
		for (auto j = 0; j < mX_in.ncol(); j++)
			mX(i, j) = mX_in(i, j);
	NumericVector modelpriorvec_r(0);
	if (!modelpriorvec_in.isNull()) {
		modelpriorvec_r = modelpriorvec_in.get();
	}
	VectorXd modelpriorvec(modelpriorvec_r.length());
	for (auto i = 0; i < modelpriorvec_r.length(); i++)
		modelpriorvec(i) = modelpriorvec_r(i);
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
	// vector< vector< dbitset > > trajectory;
	// vector< VectorXd > trajectory_probs;
	vector< MatrixXd > mXTX_inv(K);
	VectorXd sigma2(K);
	// std::unordered_map< std::size_t, bool > hash;
	vector<dbitset> vm(K);
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

	#ifdef _OPENMP
		omp_set_num_threads(cores);
	#endif

	#pragma omp parallel for ordered
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
			// #pragma omp ordered
			// hash.insert({boost::hash_value(gamma[k]), true});
			#pragma omp ordered
			vm[k] = gamma[k];
		}
	}

	// Initialise mXTX_inv
	// Initialise sigma2
	#pragma omp parallel for
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
	// trajectory.push_back(gamma);
	// trajectory_probs.push_back(log_probs.array().exp());

	// Loop until convergence
	auto converged = false;
	auto iteration = 1;
	while (!converged) {
		#ifdef DEBUG
		Rcpp::Rcout << "Iteration " << iteration << std::endl;
		#endif

		#pragma omp parallel for ordered
		for (auto k = 0; k < K; k++) {
			#ifdef DEBUG
			Rcpp::Rcout << "gamma[" << k << "] " << gamma[k] << std::endl;
			#endif
			// Try to update the 0s
			for (auto j = 0; j < p; j++) {
				if (gamma[k][j])
					continue;
				dbitset gamma_prime = gamma[k];
				gamma_prime[j] = !gamma_prime[j];
				// #pragma omp ordered {
					auto found = find(vm.begin(), vm.end(), gamma_prime);
					if (found != vm.end()) continue;
				// }
				auto p_gamma = gamma[k].count();
				auto p_gamma_prime = gamma_prime.count();
				if ((p_gamma_prime == 0) || (p_gamma_prime >= n - 1))
					continue;
				bool bUpdate = !gamma[k][j];

				// If we've seen this bitstring before, don't do the update
				// if (bUnique) {
				// 	auto h = boost::hash_value(gamma_prime);
				// 	// #pragma omp ordered
				// 	auto search = hash.find(h);
				// 	if (search != hash.end()) {
				// 		continue;
				// 	}	else {
				// 		#pragma omp ordered
				// 		hash.insert({h, true});
				// 	}
				// }
				#ifdef DEBUG
				Rcpp::Rcout << "Updating " << j << std::endl;
				#endif

				// Update or downdate mXTX_inv
				MatrixXd mXTX_inv_prime(p_gamma_prime, p_gamma_prime);
				calculate_mXTX_inv_prime(gamma[k], gamma_prime, j, mXTX, mXTX_inv[k], mXTX_inv_prime, bUpdate);

				// Calculate sigma2_prime
				double sigma2_prime = calculate_sigma2_prime(n, p_gamma_prime, mX, gamma_prime, vy, mXTX_inv_prime);

				#ifdef DEBUG
				Rcpp::Rcout << "sigma2[k] " << sigma2[k] << std::endl;
				Rcpp::Rcout << "sigma2_prime " << sigma2_prime << std::endl;
				#endif
				double log_p_0;
				double log_p_1;
				log_p_0 = calculate_log_prob(n, p, 1. - sigma2[k], p_gamma, gamma[k], log_prob, modelprior, modelpriorvec);
				log_p_1 = calculate_log_prob(n, p, 1. - sigma2_prime, p_gamma_prime, gamma_prime, log_prob, modelprior, modelpriorvec);
				#ifdef DEBUG
				Rcpp::Rcout << "log_p_0 " << log_p_0;
				Rcpp::Rcout << " log_p_1 " << log_p_1 << std::endl;
				#endif
				if (log_p_1 > log_p_0 && bUpdate) {
					if (bUnique) {
						// #pragma omp ordered
						// hash.erase(boost::hash_value(gamma[k]));
					}
					gamma[k][j] = bUpdate;
					if (bUnique) {
						// #pragma omp ordered
						// hash.insert({boost::hash_value(gamma[k]), true});
						#pragma omp ordered
						vm[k] = gamma[k];
					}
					#ifdef DEBUG
					Rcpp::Rcout << "Keep update" << std::endl;
					#endif
					sigma2[k] = sigma2_prime;
					mXTX_inv[k] = mXTX_inv_prime;
				} else {
					#ifdef DEBUG
					Rcpp::Rcout << "Don't keep update" << std::endl;
					#endif
				}
			}

			// Try to downdate the 1s
			for (auto j = 0; j < p; j++) {
				if (!gamma[k][j])
					continue;
				dbitset gamma_prime = gamma[k];
				gamma_prime[j] = !gamma_prime[j];
				// #pragma omp ordered {
					auto found = find(vm.begin(), vm.end(), gamma_prime);
					if (found != vm.end()) continue;
				// }
				auto p_gamma = gamma[k].count();
				auto p_gamma_prime = gamma_prime.count();
				if ((p_gamma_prime == 0) || (p_gamma_prime >= n - 1))
					continue;
				bool bUpdate = !gamma[k][j];

				// If we've seen this bitstring before, don't do the update
				// if (bUnique) {
				// 	auto h = boost::hash_value(gamma_prime);
				// 	// #pragma omp ordered
				// 	auto search = hash.find(h);
				// 	if (search != hash.end()) {
				// 		continue;
				// 	}	else {
				// 		#pragma omp ordered
				// 		hash.insert({h, true});
				// 	}
				// }
				#ifdef DEBUG
				Rcpp::Rcout << "Downdating " << j << std::endl;
				#endif

				// Downdate mXTX_inv
				MatrixXd mXTX_inv_prime(p_gamma_prime, p_gamma_prime);
				calculate_mXTX_inv_prime(gamma[k], gamma_prime, j, mXTX, mXTX_inv[k], mXTX_inv_prime, bUpdate);

				// Calculate sigma2_prime
				double sigma2_prime = calculate_sigma2_prime(n, p_gamma_prime, mX, gamma_prime, vy, mXTX_inv_prime);

				#ifdef DEBUG
				Rcpp::Rcout << "sigma2[k] " << sigma2[k] << std::endl;
				Rcpp::Rcout << "sigma2_prime " << sigma2_prime << std::endl;
				#endif
				double log_p_0;
				double log_p_1;
				log_p_0 = calculate_log_prob(n, p, 1. - sigma2_prime, p_gamma_prime, gamma_prime, log_prob, modelprior, modelpriorvec);
				log_p_1 = calculate_log_prob(n, p, 1. - sigma2[k], p_gamma, gamma[k], log_prob, modelprior, modelpriorvec);
				#ifdef DEBUG
				Rcpp::Rcout << "log_p_0 " << log_p_0;
				Rcpp::Rcout << " log_p_1 " << log_p_1 << std::endl;
				#endif
				if (log_p_0 > log_p_1 && !bUpdate) {
					if (bUnique) {
						// #pragma omp ordered
						// hash.erase(boost::hash_value(gamma[k]));
					}
					gamma[k][j] = bUpdate;
					if (bUnique) {
						// #pragma omp ordered
						// hash.insert({boost::hash_value(gamma[k]), true});
						#pragma omp ordered
						vm[k] = gamma[k];
					}
					#ifdef DEBUG
					Rcpp::Rcout << "Keep downdate" << std::endl;
					#endif
					sigma2[k] = sigma2_prime;
					mXTX_inv[k] = mXTX_inv_prime;
				} else {
					#ifdef DEBUG
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
		// trajectory.push_back(gamma);
		// trajectory_probs.push_back(log_probs.array().exp());
	}

	#ifdef DEBUG
	Rcpp::Rcout << "Converged" << std::endl;
	#endif
	MatrixXd mGamma_prime(K, p);
	gamma_to_MatrixXd(gamma, mGamma_prime);

	// List trajectory_bitstrings;
	// NumericMatrix trajectory_probabilities(K, trajectory.size());
	// for (auto i = 0; i < trajectory.size(); i++) {
	// 	NumericMatrix bitstrings2(K, p);
	// 	gamma_to_NumericMatrix(trajectory[i], bitstrings2);
	// 	trajectory_bitstrings.push_back(bitstrings2);
	// 	for (auto k = 0; k < K; k++) {
	// 		trajectory_probabilities(k, i) = trajectory_probs[i](k);
	// 	}
	// }

	// TODO: Add vBF, vR2, posterior inclusion probabilities, vp

	VectorXd vp_gamma(K);
	for (auto k = 0; k < K; k++) {
		vp_gamma(k) = gamma[k].count();
	}

	auto M = log_probs.array().maxCoeff(); // Maximum log-likelihood
	// Rcpp::Rcout << "M " << M << std::endl;
	VectorXd vmodel_prob = (log_probs.array() - M).array().exp() / (log_probs.array() - M).array().exp().sum();
	VectorXd vinclusion_prob = mGamma_prime.transpose() * vmodel_prob;

	uint max_idx;
	auto max_prob = vmodel_prob.maxCoeff(&max_idx);
	VectorXd vgamma_hat = mGamma.row(max_idx);

	List result = List::create(Named("mGamma") = mGamma_prime,
														 Named("vgamma.hat") = vgamma_hat,
														 Named("vBF") = log_probs,
														 Named("posterior_model_probabilities") = vmodel_prob,
														 Named("posterior_inclusion_probabilities") = vinclusion_prob,
														 Named("vR2") = 1. - sigma2.array(),
														 Named("vp") = vp_gamma);
														 // Named("trajectory") = trajectory_bitstrings,
														 // Named("trajectory_probs") = trajectory_probabilities);
	return result;
}
