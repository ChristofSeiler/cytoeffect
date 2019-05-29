/*
 * Multivariate Poisson-log Normal model
 * Author: Christof Seiler
 */
data {
  int<lower=1> n; // num of cells
  int<lower=1> d; // num of markers
  int<lower=1> p; // num of explanatory variables (including intercept)
  int<lower=0> Y[n,d]; // observed cell counts
  //matrix[n,p] X; // design matrix
  int<lower=1> k; // number of donors
  int<lower=1,upper=k> donor[n]; // donor indicator
  int<lower=1,upper=p> term[n]; // donor indicator
}
transformed data {
  real eta = 1.0; // parameter of lkj prior
  row_vector[d] zeros = rep_row_vector(0, d);
}
parameters {
  vector[p] beta[d]; // fixed coefficients
  vector<lower=0>[d] sigma; // random effects std
  vector<lower=0>[d] sigma_term; // random effects std
  vector<lower=0>[d] sigma_donor; // random effects std
  cholesky_factor_corr[d] L; // cholesky factor protein effects
  cholesky_factor_corr[d] L_term; // cholesky factor protein effects
  // cholesky_factor_corr[d] L_donor; // cholesky factor protein effects
  vector[d] z[n]; // random effects
  vector[d] z_term[n]; // random effects
  // vector[d] z_donor[k]; // random effects
  vector[d] b_donor[k]; // random effects
}
transformed parameters {
  vector[d] b[n]; // random effects
  //vector[d] b_term[n]; // random effects
  // vector[d] b_donor[k]; // random effects
  {
    matrix[d,d] Sigma; // random effects cov matrix
    matrix[d,d] Sigma_term; // random effects cov matrix
    Sigma = diag_pre_multiply(sigma, L);
    Sigma_term = diag_pre_multiply(sigma_term, L_term);
    for (i in 1:n) {
      if (term[i] == 1)
        b[i] = Sigma * z[i];
      else
        b[i] = Sigma_term * z_term[i];
    }
  }
  // {
  //   matrix[d,d] Sigma; // random effects cov matrix
  //   Sigma = diag_pre_multiply(sigma_term, L_term);
  //   for (i in 1:n)
  //     b_term[i] = X[i,2] * Sigma * z_term[i];
  // }
  // {
  //   matrix[d,d] Sigma; // random effects cov matrix
  //   Sigma = diag_pre_multiply(sigma_donor, L_donor);
  //   for (i in 1:k)
  //     b_donor[i] = Sigma * z_donor[i];
  // }
}
model {
  // priors
  for (j in 1:d)
    beta[j] ~ normal(0, 7);
  sigma ~ cauchy(0, 2.5);
  sigma_term ~ cauchy(0, 2.5);
  sigma_donor ~ cauchy(0, 2.5);
  L ~ lkj_corr_cholesky(eta);
  L_term ~ lkj_corr_cholesky(eta);
  // L_donor ~ lkj_corr_cholesky(eta);
  for (i in 1:n) {
    z[i] ~ std_normal();
    z_term[i] ~ std_normal();
  }
  for (i in 1:k)
    b_donor[i] ~ normal(zeros, sigma_donor);
  // likelihood
  for (j in 1:d) {
    // Y[i,j] ~ poisson_log(
    //   X[i] * beta[j] +
    //   gamma_beta_cov_x[i,j] +
    //   b[i,j] +
    //   b_donor[donor[i],j]
    // );
    Y[,j] ~ poisson_log(
      //X * beta[j] +
      beta[j,term] +
      to_vector(b[,j]) +
      //to_vector(b_term[,j]) +
      to_vector(b_donor[donor,j])
    );
  }
}
generated quantities {
  // in-sample prediction
  // int<lower=0> Y_hat[n,d];
  // correlation matrix
  matrix[d,d] Cor;
  matrix[d,d] Cor_term;
  matrix[d,d] Cor_donor;
  Cor = L * L';
  Cor_term = L_term * L_term';
  // Cor_donor = L_donor * L_donor';
  // for (j in 1:d) {
  //   Y_hat[,j] = poisson_log_rng(
  //     //X * beta[j] +
  //     beta[j,term] +
  //     to_vector(b[,j]) +
  //     //to_vector(b_term[,j]) +
  //     to_vector(b_donor[donor,j])
  //   );
  // }
}
