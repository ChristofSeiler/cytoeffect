/*
 * Multivariate Poisson-log Normal model
 * Author: Christof Seiler
 */
data {
  int<lower=1> n; // num of cells
  int<lower=1> d; // num of markers
  int<lower=1> p; // num of explanatory variables (including intercept)
  int<lower=0> Y[n,d]; // observed cell counts
  matrix[n,p] X; // design matrix
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
  for (j in 1:d)
    beta[j] ~ normal(0, 7);
  sigma ~ cauchy(0, 2.5);
  sigma_donor ~ cauchy(0, 2.5);
  L ~ lkj_corr_cholesky(eta);
  L_donor ~ lkj_corr_cholesky(eta);
  for (i in 1:n)
    z[i] ~ normal(0, 1);
  for (i in 1:k)
    z_donor[i] ~ normal(0, 1);
  // likelihood
  for (j in 1:d) {
    // Y[i,j] ~ poisson_log(
    //   X[i] * beta[j] +
    //   gamma_beta_cov_x[i,j] +
    //   b[i,j] +
    //   b_donor[donor[i],j]
    // );
    Y[,j] ~ poisson_log(
      X * beta[j] +
      to_vector(b[,j]) +
      to_vector(b_donor[donor,j])
    );
  }
}
generated quantities {
  // in-sample prediction
  int<lower=0> Y_hat[n,d];
  // correlation matrix
  matrix[d,d] Cor;
  matrix[d,d] Cor_donor;
  //matrix[d,d] Cor_term;
  Cor = L * L';
  Cor_donor = L_donor * L_donor';
  for (j in 1:d) {
    Y_hat[,j] = poisson_log_rng(
      X * beta[j] +
      to_vector(b[,j]) +
      to_vector(b_donor[donor,j])
    );
  }
}
