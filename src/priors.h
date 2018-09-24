// priors.h

#pragma once

#include <string>

using std::string;

double BIC(const int n, const int p_gamma, double R2);
double ZE(const int n, const int p_gamma, double R2);
double log_hyperg_2F1(double b, double c, double x);
double log_hyperg_2F1_naive(double b, double c, double x);
double liang_g1(const int n, const int p_gamma, double R2);
double liang_g2(const int n, const int p_gamma, double R2);
double liang_g_n_appell(const int n, const int p_gamma, double R2);
double liang_g_n_quad(const int n, const int p_gamma, double R2);
double liang_g_n_approx(const int n, const int p_gamma, double R2);
double robust_bayarri1(const int n, const int p_gamma, double R2);
double robust_bayarri2(const int n, const int p_gamma, double R2);

typedef std::function<double (const int N, const int p_gamma, double R2)>  log_prob_fn;
void set_log_prob(const string prior, log_prob_fn& log_prob);
