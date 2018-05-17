/*
 * Generalized Linear Mixed Model with donor and cell type random effects
 * Author: Christof Seiler
 */
data {
  int<lower=0> N; // number of cells
  int<lower=1> P; // number of fixed effects
  int<lower=1> J; // number of donors
  int<lower=1> K; // number of cell types
  int<lower=1,upper=J> donor[N]; // donor indicator
  int<lower=1,upper=K> celltype[N]; // cell type indicator
  row_vector[P] X[N]; // design matrix
  int<lower=0,upper=1> treatment[N]; // infection status
  real<lower=0> eta; // parameter of lkj prior
}
parameters {
  vector[P] beta; // fixed coefficients
  cholesky_factor_corr[P] Ld; // cholesky factor of donor random effects corr matrix
  vector<lower=0>[P] sigmad; // donor random effects std
  vector[P] zd[J]; // donor random effects
  cholesky_factor_corr[P] Lc; // cholesky factor of cell type random effects corr matrix
  vector<lower=0>[P] sigmac; // cell type random effects std
  vector[P] zc[K]; // cell type random effects
}
transformed parameters {
  vector[P] ud[J]; // donor random effects
  vector[P] uc[K]; // cell type random effects
  {
    matrix[P,P] Sigmad; // donor random effects cov matrix
    Sigmad = diag_pre_multiply(sigmad, Ld);
    for (j in 1:J)
      ud[j] = Sigmad * zd[j];
  }
  {
    matrix[P,P] Sigmac; // cell type random effects cov matrix
    Sigmac = diag_pre_multiply(sigmac, Lc);
    for (k in 1:K)
      uc[k] = Sigmac * zc[k];
  }
}
model {
  // priors
  Ld ~ lkj_corr_cholesky(eta);
  Lc ~ lkj_corr_cholesky(eta);
  for (j in 1:J)
    zd[j] ~ normal(0,1);
  for (k in 1:K)
    zc[k] ~ normal(0,1);
  // likelihood
  for (i in 1:N)
    treatment[i] ~ bernoulli_logit(X[i] * beta + X[i] * ud[donor[i]] + X[i] * uc[celltype[i]]);
}
generated quantities {
  // donor correlation matrix
  matrix[P,P] Cord;
  matrix[P,P] Corc;
  Cord = Ld * Ld';
  Corc = Lc * Lc';
}
