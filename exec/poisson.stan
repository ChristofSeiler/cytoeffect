/*
 * Multivariate Poisson-log Normal model
 * Author: Christof Seiler
 */
 functions {
   matrix L_cov_exp_quad(matrix Dmat, real alphasq, real rhosq, real delta) {
    int s = dims(Dmat)[1];
    matrix[s,s] exp_Dmat = exp(-rhosq*square(Dmat));
    matrix[s,s] K;
    for (i in 1:(s-1)) {
      K[i,i] = alphasq + delta;
      for (j in (i+1):s) {
        K[i,j] = alphasq * exp_Dmat[i,j];
        K[j,i] = K[i,j];
      }
    }
    K[s,s] = alphasq + delta;
    return cholesky_decompose(K);
  }
}
data {
  int<lower=1> n; // num of cells
  int<lower=1> d; // num of markers
  int<lower=1> p; // num of explanatory variables (including intercept)
  int<lower=0> Y[n,d]; // observed cell counts
  matrix[n,p] X; // design matrix
  int<lower=1> k; // number of donors
  int<lower=1,upper=k> donor[n]; // donor indicator
  real<lower=0> eta; // parameter of lkj prior
  int<lower=1> s; // spatial locations
  int<lower=1,upper=s> subtype[n]; // cell subtype indicator
  matrix[s,s] Dmat;
}
transformed data {
  real delta = 1e-9;
}
parameters {
  vector[p] beta[d]; // fixed coefficients
  vector<lower=0>[d] sigma; // random effects std
  vector<lower=0>[d] sigma_term; // random effects std
  vector<lower=0>[d] sigma_donor; // random effects std
  cholesky_factor_corr[d] L; // cholesky factor protein effects
  cholesky_factor_corr[d] L_term; // cholesky factor protein effects
  cholesky_factor_corr[d] L_donor; // cholesky factor protein effects
  vector[d] z[n]; // random effects
  vector[d] z_term[n]; // random effects
  vector[d] z_donor[k]; // random effects
  vector[s] z_spatial; // random effects
  real<lower=0> alphasq;
  real<lower=0> rhosq;
}
transformed parameters {
  vector[d] b[n]; // random effects
  vector[d] b_term[n]; // random effects
  vector[d] b_donor[k]; // random effects
  vector[s] b_spatial; // random effects
  {
    matrix[d,d] Sigma; // random effects cov matrix
    Sigma = diag_pre_multiply(sigma, L);
    for (i in 1:n)
      b[i] = Sigma * z[i];
  }
  {
    matrix[d,d] Sigma; // random effects cov matrix
    Sigma = diag_pre_multiply(sigma_term, L_term);
    for (i in 1:n)
      b_term[i] = X[i,2] * Sigma * z_term[i];
  }
  {
    matrix[d,d] Sigma; // random effects cov matrix
    Sigma = diag_pre_multiply(sigma_donor, L_donor);
    for (i in 1:k)
      b_donor[i] = Sigma * z_donor[i];
  }
  {
    matrix[s,s] Sigma; // random effects cov matrix
    Sigma = L_cov_exp_quad(Dmat, alphasq, rhosq, delta);
    b_spatial = Sigma * z_spatial;
  }
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
  L_donor ~ lkj_corr_cholesky(eta);
  for (i in 1:n)
    z[i] ~ std_normal();
  for (i in 1:n)
    z_term[i] ~ std_normal();
  for (i in 1:k)
    z_donor[i] ~ std_normal();
  z_spatial ~ std_normal();
  alphasq ~ cauchy(0, 2.5);
  rhosq ~ cauchy(0, 2.5);
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
      to_vector(b_term[,j]) +
      to_vector(b_donor[donor,j]) +
      b_spatial[subtype[s]]
    );
  }
}
generated quantities {
  // in-sample prediction
  int<lower=0> Y_hat[n,d];
  // correlation matrix
  matrix[d,d] Cor;
  matrix[d,d] Cor_term;
  matrix[d,d] Cor_donor;
  //matrix[d,d] Cor_term;
  Cor = L * L';
  Cor_term = L_term * L_term';
  Cor_donor = L_donor * L_donor';
  for (j in 1:d) {
    Y_hat[,j] = poisson_log_rng(
      X * beta[j] +
      to_vector(b[,j]) +
      to_vector(b_term[,j]) +
      to_vector(b_donor[donor,j]) +
      b_spatial[subtype[s]]
    );
  }
}
