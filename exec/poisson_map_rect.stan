/*
 * Multivariate Poisson-log Normal model
 * Author: Christof Seiler
 */
data {
  int<lower=1> n; // num of cells
  int<lower=1> d; // num of markers
  int<lower=1> p; // num of explanatory variables (including intercept)
  int<lower=0> Y[n,d]; // observed cell counts
  row_vector[p] X[n]; // design matrix
  int<lower=1> k; // number of donors
  int<lower=1,upper=k> donor[n]; // donor indicator
  real<lower=0> eta; // parameter of lkj prior
}
parameters {
  vector[p] beta[d]; // fixed coefficients
  vector<lower=0>[d] sigma; // random effects std
  vector<lower=0>[d] sigma_donor; // random effects std
  cholesky_factor_corr[d] L; // cholesky factor protein effects
  cholesky_factor_corr[d] L_donor; // cholesky factor protein effects
  vector[d] z[n]; // random effects
  vector[d] z_donor[k]; // random effects
}
transformed parameters {
  vector[d] b[n]; // random effects
  vector[d] b_donor[k]; // random effects
  {
    matrix[d,d] Sigma; // random effects cov matrix
    Sigma = diag_pre_multiply(sigma, L);
    for (i in 1:n)
      b[i] = Sigma * z[i];
  }
  {
    matrix[d,d] Sigma; // random effects cov matrix
    Sigma = diag_pre_multiply(sigma_donor, L_donor);
    for (i in 1:k)
      b_donor[i] = Sigma * z_donor[i];
  }
}
model {
  // priors
  //for (j in 1:d)
  //  beta[j] ~ normal(0, 5);
  sigma ~ cauchy(0, 2.5);
  sigma_donor ~ cauchy(0, 2.5);
  L ~ lkj_corr_cholesky(eta);
  L_donor ~ lkj_corr_cholesky(eta);
  for (i in 1:n)
    z[i] ~ normal(0, 1);
  for (i in 1:k)
    z_donor[i] ~ normal(0, 1);
  // likelihood
  for (i in 1:n) {
    for (j in 1:d) {
      target += poisson_log_lpmf(Y[i,j] | X[i] * beta[j] + b[i,j] + b_donor[donor[i],j]);
    }
  }
}
generated quantities {
  // correlation matrix
  matrix[d,d] Cor;
  matrix[d,d] Cor_donor;
  Cor = L * L';
  Cor_donor = L_donor * L_donor';
}
