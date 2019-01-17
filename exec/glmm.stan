/*
 * Generalized Linear Mixed Model with donor and cell type random effects
 * Author: Christof Seiler
 */
data {
  int<lower=0> N; // number of cells
  int<lower=1> P; // number of fixed effects
  int<lower=1> J; // number of donors
  int<lower=1,upper=J> donor[N]; // donor indicator
  row_vector[P] X[N]; // design matrix
  int<lower=0,upper=1> treatment[N]; // infection status
  real<lower=0> eta; // parameter of lkj prior
}
parameters {
  vector[P] beta; // fixed coefficients
  cholesky_factor_corr[P] L_donor; // cholesky factor of donor random effects corr matrix
  vector<lower=0>[P] sigma_donor; // donor random effects std
  vector[P] z_donor[J]; // donor random effects
}
transformed parameters {
  vector[P] b_donor[J]; // random effects
  {
    matrix[P,P] Sigma; // random effects cov matrix
    Sigma = diag_pre_multiply(sigma_donor, L_donor);
    for (i in 1:J)
      b_donor[i] = Sigma * z_donor[i];
  }
}
model {
  // priors
  beta ~ normal(0, 10);
  L_donor ~ lkj_corr_cholesky(eta);
  sigma_donor ~ cauchy(0, 2.5);
  for (j in 1:J)
    z_donor[j] ~ normal(0,1);
  // likelihood
  {
    vector[N] x_beta_b_donor;
    for (i in 1:N)
      x_beta_b_donor[i] = X[i] * beta + X[i] * b_donor[donor[i]];
    treatment ~ bernoulli_logit(x_beta_b_donor);
  }
}
generated quantities {
  // donor correlation matrix
  matrix[P,P] Cor_donor;
  Cor_donor = L_donor * L_donor';
}
