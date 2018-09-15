// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <RcppGSL.h>
#include <Rcpp.h>

using namespace Rcpp;

// maruyama
double maruyama(const int n, const int p, const double R2, const int p_gamma);
RcppExport SEXP _blma_maruyama(SEXP nSEXP, SEXP pSEXP, SEXP R2SEXP, SEXP p_gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< const int >::type p_gamma(p_gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(maruyama(n, p, R2, p_gamma));
    return rcpp_result_gen;
END_RCPP
}
// BIC
double BIC(const int n, const int p, double R2, const int p_gamma);
RcppExport SEXP _blma_BIC(SEXP nSEXP, SEXP pSEXP, SEXP R2SEXP, SEXP p_gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< const int >::type p_gamma(p_gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(BIC(n, p, R2, p_gamma));
    return rcpp_result_gen;
END_RCPP
}
// ZE
double ZE(const int n, const int p, const double R2, const int p_gamma);
RcppExport SEXP _blma_ZE(SEXP nSEXP, SEXP pSEXP, SEXP R2SEXP, SEXP p_gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< const int >::type p_gamma(p_gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(ZE(n, p, R2, p_gamma));
    return rcpp_result_gen;
END_RCPP
}
// log_hyperg_2F1
double log_hyperg_2F1(double b, double c, double x);
RcppExport SEXP _blma_log_hyperg_2F1(SEXP bSEXP, SEXP cSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(log_hyperg_2F1(b, c, x));
    return rcpp_result_gen;
END_RCPP
}
// log_hyperg_2F1_naive
double log_hyperg_2F1_naive(double b, double c, double x);
RcppExport SEXP _blma_log_hyperg_2F1_naive(SEXP bSEXP, SEXP cSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(log_hyperg_2F1_naive(b, c, x));
    return rcpp_result_gen;
END_RCPP
}
// liang_g1
double liang_g1(const int n, const int p, const double R2, const int p_gamma);
RcppExport SEXP _blma_liang_g1(SEXP nSEXP, SEXP pSEXP, SEXP R2SEXP, SEXP p_gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< const int >::type p_gamma(p_gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(liang_g1(n, p, R2, p_gamma));
    return rcpp_result_gen;
END_RCPP
}
// liang_g2
double liang_g2(const int n, const int p, const double R2, const int p_gamma);
RcppExport SEXP _blma_liang_g2(SEXP nSEXP, SEXP pSEXP, SEXP R2SEXP, SEXP p_gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< const int >::type p_gamma(p_gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(liang_g2(n, p, R2, p_gamma));
    return rcpp_result_gen;
END_RCPP
}
// liang_g_n_appell
double liang_g_n_appell(const int n, const int p, const double R2, const int p_gamma);
RcppExport SEXP _blma_liang_g_n_appell(SEXP nSEXP, SEXP pSEXP, SEXP R2SEXP, SEXP p_gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< const int >::type p_gamma(p_gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(liang_g_n_appell(n, p, R2, p_gamma));
    return rcpp_result_gen;
END_RCPP
}
// liang_g_n_quad
double liang_g_n_quad(const int n, const int p, const double R2, const int p_gamma);
RcppExport SEXP _blma_liang_g_n_quad(SEXP nSEXP, SEXP pSEXP, SEXP R2SEXP, SEXP p_gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< const int >::type p_gamma(p_gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(liang_g_n_quad(n, p, R2, p_gamma));
    return rcpp_result_gen;
END_RCPP
}
// liang_g_n_approx
double liang_g_n_approx(const int n, const int p, const double R2, const int p_gamma);
RcppExport SEXP _blma_liang_g_n_approx(SEXP nSEXP, SEXP pSEXP, SEXP R2SEXP, SEXP p_gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< const int >::type p_gamma(p_gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(liang_g_n_approx(n, p, R2, p_gamma));
    return rcpp_result_gen;
END_RCPP
}
// robust_bayarri1
double robust_bayarri1(const int n, const int p, const double R2, const int p_gamma);
RcppExport SEXP _blma_robust_bayarri1(SEXP nSEXP, SEXP pSEXP, SEXP R2SEXP, SEXP p_gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< const int >::type p_gamma(p_gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(robust_bayarri1(n, p, R2, p_gamma));
    return rcpp_result_gen;
END_RCPP
}
// robust_bayarri2
double robust_bayarri2(const int n, const int p, const double R2, const int p_gamma);
RcppExport SEXP _blma_robust_bayarri2(SEXP nSEXP, SEXP pSEXP, SEXP R2SEXP, SEXP p_gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< const int >::type p_gamma(p_gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(robust_bayarri2(n, p, R2, p_gamma));
    return rcpp_result_gen;
END_RCPP
}
// log_BF_g_on_n_integrand
double log_BF_g_on_n_integrand(const double vu, const int n, const int p, const double R2, const double a);
RcppExport SEXP _blma_log_BF_g_on_n_integrand(SEXP vuSEXP, SEXP nSEXP, SEXP pSEXP, SEXP R2SEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type vu(vuSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(log_BF_g_on_n_integrand(vu, n, p, R2, a));
    return rcpp_result_gen;
END_RCPP
}
// log_BF_g_on_n_quad
double log_BF_g_on_n_quad(const int n, const int p, const double R2, const int a);
RcppExport SEXP _blma_log_BF_g_on_n_quad(SEXP nSEXP, SEXP pSEXP, SEXP R2SEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< const int >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(log_BF_g_on_n_quad(n, p, R2, a));
    return rcpp_result_gen;
END_RCPP
}
// pva
List pva(const NumericVector vy_in, const NumericMatrix mX_in, const NumericMatrix mGamma_in, const std::string prior, const std::string modelprior, const Nullable<NumericVector> modelpriorvec_in, const bool bUnique, const double lambda, const int cores);
RcppExport SEXP _blma_pva(SEXP vy_inSEXP, SEXP mX_inSEXP, SEXP mGamma_inSEXP, SEXP priorSEXP, SEXP modelpriorSEXP, SEXP modelpriorvec_inSEXP, SEXP bUniqueSEXP, SEXP lambdaSEXP, SEXP coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type vy_in(vy_inSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type mX_in(mX_inSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type mGamma_in(mGamma_inSEXP);
    Rcpp::traits::input_parameter< const std::string >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const std::string >::type modelprior(modelpriorSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector> >::type modelpriorvec_in(modelpriorvec_inSEXP);
    Rcpp::traits::input_parameter< const bool >::type bUnique(bUniqueSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const int >::type cores(coresSEXP);
    rcpp_result_gen = Rcpp::wrap(pva(vy_in, mX_in, mGamma_in, prior, modelprior, modelpriorvec_in, bUnique, lambda, cores));
    return rcpp_result_gen;
END_RCPP
}
// blma
List blma(NumericVector vy, NumericMatrix mX, std::string prior, std::string modelprior, Nullable<NumericVector> modelpriorvec, const unsigned int cores);
RcppExport SEXP _blma_blma(SEXP vySEXP, SEXP mXSEXP, SEXP priorSEXP, SEXP modelpriorSEXP, SEXP modelpriorvecSEXP, SEXP coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vy(vySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mX(mXSEXP);
    Rcpp::traits::input_parameter< std::string >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< std::string >::type modelprior(modelpriorSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type modelpriorvec(modelpriorvecSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type cores(coresSEXP);
    rcpp_result_gen = Rcpp::wrap(blma(vy, mX, prior, modelprior, modelpriorvec, cores));
    return rcpp_result_gen;
END_RCPP
}
// blma_fixed
List blma_fixed(NumericVector vy, NumericMatrix mX, NumericMatrix mZ, std::string prior, std::string modelprior, Nullable<NumericVector> modelpriorvec, const unsigned int cores);
RcppExport SEXP _blma_blma_fixed(SEXP vySEXP, SEXP mXSEXP, SEXP mZSEXP, SEXP priorSEXP, SEXP modelpriorSEXP, SEXP modelpriorvecSEXP, SEXP coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vy(vySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mX(mXSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mZ(mZSEXP);
    Rcpp::traits::input_parameter< std::string >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< std::string >::type modelprior(modelpriorSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type modelpriorvec(modelpriorvecSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type cores(coresSEXP);
    rcpp_result_gen = Rcpp::wrap(blma_fixed(vy, mX, mZ, prior, modelprior, modelpriorvec, cores));
    return rcpp_result_gen;
END_RCPP
}
// graycode
IntegerMatrix graycode(unsigned int varying, unsigned int fixed);
RcppExport SEXP _blma_graycode(SEXP varyingSEXP, SEXP fixedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type varying(varyingSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type fixed(fixedSEXP);
    rcpp_result_gen = Rcpp::wrap(graycode(varying, fixed));
    return rcpp_result_gen;
END_RCPP
}
// sampler
List sampler(const int iterations, const NumericVector vy_in, const NumericMatrix mX_in, const std::string prior, const std::string modelprior, const Nullable<NumericVector> modelpriorvec_in, const int cores);
RcppExport SEXP _blma_sampler(SEXP iterationsSEXP, SEXP vy_inSEXP, SEXP mX_inSEXP, SEXP priorSEXP, SEXP modelpriorSEXP, SEXP modelpriorvec_inSEXP, SEXP coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type vy_in(vy_inSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type mX_in(mX_inSEXP);
    Rcpp::traits::input_parameter< const std::string >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const std::string >::type modelprior(modelpriorSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector> >::type modelpriorvec_in(modelpriorvec_inSEXP);
    Rcpp::traits::input_parameter< const int >::type cores(coresSEXP);
    rcpp_result_gen = Rcpp::wrap(sampler(iterations, vy_in, mX_in, prior, modelprior, modelpriorvec_in, cores));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_blma_maruyama", (DL_FUNC) &_blma_maruyama, 4},
    {"_blma_BIC", (DL_FUNC) &_blma_BIC, 4},
    {"_blma_ZE", (DL_FUNC) &_blma_ZE, 4},
    {"_blma_log_hyperg_2F1", (DL_FUNC) &_blma_log_hyperg_2F1, 3},
    {"_blma_log_hyperg_2F1_naive", (DL_FUNC) &_blma_log_hyperg_2F1_naive, 3},
    {"_blma_liang_g1", (DL_FUNC) &_blma_liang_g1, 4},
    {"_blma_liang_g2", (DL_FUNC) &_blma_liang_g2, 4},
    {"_blma_liang_g_n_appell", (DL_FUNC) &_blma_liang_g_n_appell, 4},
    {"_blma_liang_g_n_quad", (DL_FUNC) &_blma_liang_g_n_quad, 4},
    {"_blma_liang_g_n_approx", (DL_FUNC) &_blma_liang_g_n_approx, 4},
    {"_blma_robust_bayarri1", (DL_FUNC) &_blma_robust_bayarri1, 4},
    {"_blma_robust_bayarri2", (DL_FUNC) &_blma_robust_bayarri2, 4},
    {"_blma_log_BF_g_on_n_integrand", (DL_FUNC) &_blma_log_BF_g_on_n_integrand, 5},
    {"_blma_log_BF_g_on_n_quad", (DL_FUNC) &_blma_log_BF_g_on_n_quad, 4},
    {"_blma_pva", (DL_FUNC) &_blma_pva, 9},
    {"_blma_blma", (DL_FUNC) &_blma_blma, 6},
    {"_blma_blma_fixed", (DL_FUNC) &_blma_blma_fixed, 7},
    {"_blma_graycode", (DL_FUNC) &_blma_graycode, 2},
    {"_blma_sampler", (DL_FUNC) &_blma_sampler, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_blma(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
