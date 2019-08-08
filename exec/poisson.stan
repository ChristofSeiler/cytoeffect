/*
 * Multivariate Poisson-log Normal model
 * Author: Christof Seiler
 */
functions {
  matrix polar(matrix X) {
    int d = dims(X)[1];
    int r = dims(X)[2];
    matrix[d,r] Q;
    vector[r] eval;
    vector[r] eval_trans;
    matrix[r,r] evec;

    eval = eigenvalues_sym(X'*X);
    for(l in 1:r)
      eval_trans[l] = 1/sqrt(eval[l]);
    evec = eigenvectors_sym(X'*X);
    Q = X*evec*diag_matrix(eval_trans)*evec';
    return Q;
  }
  matrix cov2cor(matrix Q, vector sigma) {
    int d = dims(Q)[1];
    vector[d] sig;
    matrix[d,d] Cov;
    matrix[d,d] Cor;
    Cov = quad_form(diag_matrix(square(sigma)), Q');
    sig = diagonal(Cov);
    for(i in 1:d)
      sig[i] = 1/sqrt(sig[i]);
    Cor = quad_form_diag(Cov, sig);
    return Cor;
  }
  int num_zeros(int[] y) {
    int counter = 0;
    for (i in 1:size(y))
      counter += (y[i] == 0);
    return counter;
  }
}
data {
  int<lower=1> n; // num of cells
  int<lower=1> d; // num of markers
  int<lower=1> p; // num of explanatory variables (including intercept)
  int<lower=0> Y[d,n]; // observed cell counts
  int<lower=1> k; // number of donors
  int<lower=1,upper=k> donor[n]; // donor indicator
  int<lower=1,upper=p> term[n]; // condition indicator
  int<lower=1> r_donor; // rank of latent matrix
}
transformed data {
  real eta = 1.0; // parameter of lkj prior
  int<lower=0> y[d*n] = to_array_1d(Y);
  int<lower=0> n_zero = num_zeros(y);
  int<lower=0> n_nonzero = size(y) - n_zero;
  int<lower=0,upper=size(y)> indices_zero[n_zero];
  int<lower=0,upper=size(y)> indices_nonzero[n_nonzero];
  int<lower=0,upper=size(y)> n_zero_current = 0;
  int<lower=0,upper=size(y)> n_nonzero_current = 0;
  for (i in 1:size(y)) {
    if (y[i] == 0) {
      n_zero_current += 1;
      indices_zero[n_zero_current] = i;
    }
    else {
      n_nonzero_current += 1;
      indices_nonzero[n_nonzero_current] = i;
    }
  }
}
parameters {
  vector[n*d] z; // random effects
  matrix[d,p] beta; // fixed coefficients
  vector<lower=0>[d] sigma; // random effects std
  vector<lower=0>[r_donor] sigma_donor; // random effects std
  cholesky_factor_corr[d] L; // cholesky factor protein effects
  vector[k*r_donor] z_donor; // random effects
  vector[d*r_donor] x_donor; // distribution on X for polar expansion
  real<lower=0, upper=1> theta[d,k]; // mixing proportions
}
transformed parameters {
  matrix[k,d] b_donor;
  matrix[d,r_donor] Q_donor;
  matrix[r_donor,d] Sigma_donor_t;
  {
    matrix[d,r_donor] X;
    matrix[k,r_donor] Z;
    X = to_matrix(x_donor, d, r_donor);
    Q_donor = polar(X);
    Sigma_donor_t = diag_post_multiply(Q_donor, sigma_donor)';
    Z = to_matrix(z_donor, k, r_donor);
    b_donor = Z * Sigma_donor_t;
  }
}
model {
  // priors
  for (j in 1:d)
    beta[j] ~ normal(0, 7);
  sigma ~ cauchy(0, 2.5);
  sigma_donor ~ cauchy(0, 2.5);
  L ~ lkj_corr_cholesky(eta);
  z ~ std_normal();
  z_donor ~ std_normal();
  x_donor ~ std_normal();
  {
    matrix[d,d] Sigma_t;
    matrix[n,d] b;
    matrix[n,d] Z;
    vector[n*d] lambda;
    real theta_vec[n*d];
    Sigma_t = diag_pre_multiply(sigma, L)';
    Z = to_matrix(z, n, d);
    b = Z * Sigma_t;
    // likelihood
    lambda = to_vector(beta'[term,] + b + b_donor[donor,]);
    theta_vec = to_array_1d(theta[,donor]);
    for (i in 1:n_zero) {
      // mixtures cannot be vectorized
      real theta_current = theta_vec[indices_zero[i]];
      target += log_sum_exp(bernoulli_lpmf(1 | theta_current),
                            bernoulli_lpmf(0 | theta_current) +
                            poisson_log_lpmf(0 | lambda[indices_zero[i]]));
    }
    target += bernoulli_lpmf(0 | theta_vec[indices_nonzero]);
    target += poisson_log_lpmf(y[indices_nonzero] | lambda[indices_nonzero]);
  }
}
generated quantities {
  // correlation matrices
  matrix[d,d] Cor;
  matrix[d,d] Cor_donor;
  Cor = L * L';
  Cor_donor = cov2cor(Q_donor, sigma_donor);
}
