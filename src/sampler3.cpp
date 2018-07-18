VectorXd vy;
MatrixXd mX;

int n = vy.size();
int p = mX.cols();
MatrixXd mXTX = mX.transpose() * mX;
Matrix mXTy = mX.transpose() * vy;
int q;
double R2;
MatrixXd mGamma(iterations, p);
dbitset gamma(p);
double log_BF_curr;

for (auto i = 0; i < iterations; i++) {
  q = gamma.count();
  if (q == 0) {
    R2 = 0.;
  } else {
    mA = get_rows_and_cols(mXTX, gamma, gamma, mA);
    vb = get_rows(mXTy, gamma, vb);
    R2 = (vb.transpose() * mA.inverse() * vb / n).value();
  }
  log_BF_curr = calculate_log_prob(n, p, R2, q, gamma, log_prob, modelprior, modelpriorvec);

  for (auto j = 0; j < p; j++) {
    dbitset gamma_prop = gamma;
    gamma_prop[j] = !gamma[j];

    q = gamma_prop.count();
    if (q == 0) {
      R2 = 0.;
    } else {
      mA = get_rows_and_cols(mXTX, gamma_prop, gamma_prop, mA);
      vb = get_rows(mXTy, gamma_prop, vb);
      R2 = (vb.transpose() * mA.inverse() * vb / n).value();
    }
    auto log_BF_prop = calculate_log_prob(n, p, R2, q, gamma_prop, log_prob, modelprior, modelpriorvec);

    auto r = 1. / (1. + exp(log_BF_prop - log_BF_curr));

    if (r < R::runif(0., 1.)) {
      gamma = gamma_prop;
      log_BF_curr = log_BF_prop;
    }
  }

  mGamma.row(i) = gamma_to_row(gamma);
}