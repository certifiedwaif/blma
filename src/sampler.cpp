#include <Rcpp.h>
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppGSL)]]
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
#include "pva.h"

// #define DEBUG

using namespace std;
using namespace Rcpp;


VectorXd gamma_to_row(const dbitset& gamma)
{
  	auto p = gamma.size();
  	VectorXd v(p);
  	for (int j = 0; j < p; j++) {
    	v(j) = gamma[j] ? 1. : 0.;
  	}
  	return v;
}


//' sampler
//'
//' @usage sampler(iterations, vy, mX, prior = "BIC", modelprior = "uniform",
//' modelpriorvec_in=NULL, cores=1)
//' @param iterations The number of iterations to run the MCMC sampler for
//' @param vy_in Vector of responses
//' @param mX_in The matrix of covariates which may or may not be included in each model
//' @param prior -- the choice of mixture $g$-prior used to perform Bayesian model
//' averaging. The choices available include:
//'   \itemize{
//'     \item{"BIC"}{-- the Bayesian information criterion obtained by using the cake
//'     prior of Ormerod et al. (2017).}
//'
//'     \item{"ZE"}{-- special case of the prior structure described by Maruyama and
//'	George (2011).}
//'
//'     \item{"liang_g1"}{-- the mixture \eqn{g}-prior of Liang et al. (2008) with prior 
//'	hyperparameter \eqn{a=3} evaluated directly using Equation (10) of Greenaway and 
//'	Ormerod (2018) where the Gaussian hypergeometric function is evaluated using the 
//'	{gsl} library. Note: this option can lead to numerical problems and is only
//'     meant to be used for comparative purposes.}
//'
//'     \item{"liang_g2"}{-- the mixture \eqn{g}-prior of Liang et al. (2008) with prior 
//'	hyperparameter \eqn{a=3} evaluated directly using Equation (11) of Greenaway and 
//'	Ormerod (2018).}
//'
//'     \item{"liang_g_n_appell"}{-- the mixture \eqn{g/n}-prior of Liang et al. (2008)
//'	with prior hyperparameter \eqn{a=3} evaluated using the {appell R} package.}
//'
//'     \item{"liang_g_approx"}{-- the mixture \eqn{g/n}-prior of Liang et al. (2008)
//'	with prior hyperparameter eqn{a=3} using the approximation Equation (15)  of
//'	Greenaway and Ormerod (2018) for model with more than two covariates and
//'	numerical quadrature (see below) for models with one or two covariates.}
//'
//'     \item{"liang_g_n_quad"}{-- the mixture \eqn{g/n}-prior of Liang et al. (2008)
//'	with prior hyperparameter eqn{a=3} evaluated using a composite trapezoid rule.}
//'
//'     \item{"robust_bayarri1"}{-- the robust prior of Bayarri et al. (2012) using
//'	default prior hyper choices evaluated directly using Equation (18)  of Greenaway 
//'	and Ormerod (2018) with the {gsl} library.}
//'
//'     \item{"robust_bayarri2"}{-- the robust prior of Bayarri et al. (2012) using 
//'	default prior hyper choices evaluated directly using Equation (19) of Greenaway
//'	and Ormerod (2018).}
//' }
//' @param modelprior The model prior to use. The choices of model prior are "uniform",
//' "beta-binomial" or "bernoulli". The choice of model prior dictates the meaning of the
//' parameter modelpriorvec.
//' @param modelpriorvec If modelprior is "uniform", then the modelpriorvec is ignored
//' and can be null.
//'
//' If modelprior is "beta-binomial" then modelpriorvec should be length 2 with the first
//' element containing alpha hyperparameter for the beta prior and the second element
//' containing the beta hyperparameter for beta prior.
//'
//' If modelprior is "bernoulli", then modelpriorvec must be of the same length as the
//' number columns in mX. Each element i of modelpriorvec contains the prior probability 
//' of the the ith covariate being included in the model.
//' @param cores The number of cores to use. Defaults to 1.
//' @return The object returned is a list containing:
//' \itemize{
//'   \item{"mGamma"}{-- An iterations by p binary matrix containing the sampled models.}
//'   \item{"vinclusion_prob"}{-- The vector of inclusion probabilities.}
//'   \item{"vlogBF"}{-- The vector of logs of the Bayes Factors of the models in mGamma.}
//' }
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
//' vy <- mD$y
//' mX <- data.matrix(cbind(mD[1:15]))
//' colnames(mX) <- varnames
//' K <- 100
//' p <- ncol(mX)
//' sampler_result <- sampler(10000, vy, mX, prior = "BIC",
//'					  modelprior = "uniform",
//'                                   modelpriorvec_in=NULL)
//'
//' @references
//' Bayarri, M. J., Berger, J. O., Forte, A., Garcia-Donato, G., 2012. Criteria for
//' Bayesian model choice with application to variable selection. Annals of Statistics
//' 40 (3), 1550-1577.
//'
//' Greenaway, M. J., J. T. Ormerod (2018) Numerical aspects of Bayesian linear models
//' averaging using mixture g-priors.
//'
//' Liang, F., Paulo, R., Molina, G., Clyde, M. a., Berger, J. O., 2008. Mixtures of g
//' priors for Bayesian variable selection. Journal of the American Statistical
//' Association 103 (481), 410-423.
//'
//' Ormerod, J. T., Stewart, M., Yu, W., Romanes, S. E., 2017. Bayesian hypothesis tests
//' with diffuse priors: Can we have our cake and eat it too?
//' @export
// [[Rcpp::export]]
List sampler(const int iterations,
        const NumericVector vy_in,
        const NumericMatrix mX_in,
        const std::string prior,
        const std::string modelprior,
        const Nullable<NumericVector> modelpriorvec_in = R_NilValue,
        const int cores = 1L)
{
  	// Try using the parallelisation in Eigen. This is an inherently serial algorithm,
  	// and I don't think OpenMP is going to help us here. More cores might
  	// help us with the linear algebra and numerical integration though.
#ifdef _OPENMP
	Eigen::initParallel();
    omp_set_num_threads(cores);
	Eigen::setNbThreads(cores);
#endif

  	VectorXd vy(vy_in.length());   // = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(vy_in);
  	for (auto i = 0; i < vy_in.length(); i++) vy[i] = vy_in[i];
    // = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(mX_in);
  	MatrixXd mX(mX_in.nrow(), mX_in.ncol());
  	for (auto i = 0; i < mX_in.nrow(); i++)
    	for (auto j = 0; j < mX_in.ncol(); j++)
      		mX(i, j) = mX_in(i, j);
  	const auto n = mX.rows();
  	const auto p = mX.cols();
  	// Normalise vy and mX
  	Normed normed = normalise(vy, mX);
  	vy = normed.vy;
  	mX = normed.mX;

  	NumericVector modelpriorvec_r(0);
  	if (!modelpriorvec_in.isNull()) {
    	modelpriorvec_r = modelpriorvec_in.get();
  	}
  	VectorXd modelpriorvec(modelpriorvec_r.length());
  	for (auto i = 0; i < modelpriorvec_r.length(); i++)
    	modelpriorvec(i) = modelpriorvec_r(i);

  	MatrixXd mXTX = mX.transpose() * mX;
  	MatrixXd mXTy = mX.transpose() * vy;
  	int q;
  	double R2;
  	MatrixXd mGamma(iterations, p);
    VectorXd vlogBF(iterations);
  	dbitset gamma(p);
  	mGamma.row(0) = gamma_to_row(gamma);
  	double log_BF_curr;
  	MatrixXd mA;
  	MatrixXd mA_inv;
  	MatrixXd mA_inv_prop;
  	MatrixXd vb;
  	log_prob_fn log_prob;
  	set_log_prob(prior, log_prob);

  	q = gamma.count();
  	if (q == 0) {
    	R2 = 0.;
  	} else {
    	mA.resize(q, q);
    	mA = get_rows_and_cols(mXTX, gamma, gamma, mA);
    	vb.resize(q, 1);
    	vb = get_rows(mXTy, gamma, vb);
    	mA_inv = mA.inverse();
    	R2 = (vb.transpose() * mA_inv * vb / n).value();
  	}
#ifdef DEBUG
  	Rcpp::Rcout << "R2 " << R2 << std::endl;
#endif
  	log_BF_curr = calculate_log_prob(n, p, R2, q, gamma, log_prob, modelprior, modelpriorvec);
#ifdef DEBUG
  	Rcpp::Rcout << "log_BF_curr " << log_BF_curr << std::endl;
#endif
  	for (auto i = 0; i < iterations - 1; i++) {
    	for (auto j = 0; j < p; j++) {
      		dbitset gamma_prop = gamma;
      		gamma_prop[j] = !gamma[j];
#ifdef DEBUG
      		Rcpp::Rcout << "gamma " << gamma << std::end;
      		Rcpp::Rcout << "gamma_prop " << gamma_prop << std::endl;
#endif

      		q = gamma_prop.count();
      		if (q == 0) {
        		R2 = 0.;
      		} else {
        		mA_inv_prop.resize(q, q);
				bool bUpdate = gamma_prop.count() > gamma.count();
#ifdef DEBUG
      			Rcpp::Rcout << "gamma " << gamma << std::endl;
      			Rcpp::Rcout << " gamma_prop " << gamma_prop << std::endl;
      			Rcpp::Rcout << "j " << j << std::endl;
      			Rcpp::Rcout << "mA_inv " << mA_inv << std::endl;
      			Rcpp::Rcout << "mA_inv_prop " << mA_inv_prop << std::endl;
#endif
        		calculate_mXTX_inv_prime(gamma, gamma_prop, j, mXTX, mA_inv, mA_inv_prop, bUpdate);
        		vb.resize(q, 1);
        		vb = get_rows(mXTy, gamma_prop, vb);
        		R2 = (vb.transpose() * mA_inv_prop * vb / n).value();
      		}
#ifdef DEBUG
      		Rcpp::Rcout << "R2 " << R2 << std::endl;
#endif
      		auto log_BF_prop = calculate_log_prob(n, p, R2, q, gamma_prop, log_prob,
                    modelprior, modelpriorvec);
#ifdef DEBUG
      		Rcpp::Rcout << "log_BF_prop " << log_BF_prop << std::endl;
#endif

      		auto r = 1. / (1. + exp(log_BF_prop - log_BF_curr));

      		if (r < R::runif(0., 1.)) {
        		gamma = gamma_prop;
        		log_BF_curr = log_BF_prop;
				mA_inv = mA_inv_prop;
      		}
    	}

    	mGamma.row(i) = gamma_to_row(gamma);
      vlogBF(i) = log_BF_curr;
  	}
  	VectorXd vinclusion_prob(p);
#pragma omp parallel for \
	shared(mGamma, vinclusion_prob) \
	default(none)
  	for (auto i = 0; i < p; i++) {
  		vinclusion_prob(i) = mGamma.col(i).mean();
  	}

  	return List::create(Named("mGamma") = mGamma,
						Named("vinclusion_prob") = vinclusion_prob,
            Named("vlogBF") = vlogBF);
}
