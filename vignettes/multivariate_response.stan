/*
 * Seemingly unrelated regression with index notation
 * Author: Christof Seiler
 */
data {
  int<lower=1> d;
  int<lower=1> p;
  int<lower=0> n;
  vector[d] y[n];
  int<lower=1,upper=p> x[n];
}
parameters {
  matrix[d,p] beta;
  vector<lower=0>[d] sigma;
  cholesky_factor_corr[d] L_Omega;
}
model {
  vector[d] mu[n];
  matrix[d,d] L_Sigma;
  for (i in 1:n) {
    for (j in 1:d)
      mu[i,j] = beta[j,x[i]];
  }
  L_Sigma = diag_pre_multiply(sigma, L_Omega);
  y ~ multi_normal_cholesky(mu, L_Sigma);
}
generated quantities {
  matrix[d,d] Cor;
  Cor = L_Omega * L_Omega';
}
