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
    for(i in 1:d)
      sig[i] = 1/sqrt(sig[i]);
    Cor = quad_form_diag(Cov, sig);
    return Cor;
  }
}
data {
  int<lower=1> n; // num of cells
  int<lower=1> d; // num of markers
  int<lower=1> p; // num of explanatory variables (including intercept)
  int<lower=0> Y[n,d]; // observed cell counts
  int<lower=1> k; // number of donors
  int<lower=1,upper=k> donor[n]; // donor indicator
  int<lower=1,upper=p> term[n]; // condition indicator
  int<lower=1> r; // Rank of latent matrix
}
//transformed data {
//  real eta = 1.0; // parameter of lkj prior
//}
parameters {
  vector[p] beta[d]; // fixed coefficients
  vector<lower=0>[r] sigma; // random effects std
  vector<lower=0>[r] sigma_term; // random effects std
  vector<lower=0>[r] sigma_donor; // random effects std
  //cholesky_factor_corr[d] L; // cholesky factor protein effects
  //cholesky_factor_corr[d] L_term; // cholesky factor protein effects
  //cholesky_factor_corr[d] L_donor; // cholesky factor protein effects
  vector[d] z[n]; // random effects
  vector[d] z_term[n]; // random effects
  vector[r] z_donor[k]; // random effects
  vector[d*r] x; // distribution on X for polar expansion
  vector[d*r] x_term; // distribution on X for polar expansion
  vector[d*r] x_donor; // distribution on X for polar expansion
}
transformed parameters {
  vector[d] b[n]; // random effects
  vector[d] b_donor[k]; // random effects
  matrix[d,r] Q;
  matrix[d,r] Q_term;
  matrix[d,r] Q_donor;
  {
    matrix[d,r] Sigma;
    matrix[d,r] Sigma_term;
    matrix[d,r] Sigma_donor;
    matrix[d,r] X;
    X = to_matrix(x, d, r);
    Q = polar(X);
    Sigma = diag_post_multiply(Q, sigma);
    X = to_matrix(x_term, d, r);
    Q = polar(X);
    Sigma_term = diag_post_multiply(Q_term, sigma_term);
    for (i in 1:n) {
      if (term[i] == 1)
        b[i] = Sigma * z[i];
      else
        b[i] = Sigma_term * z_term[i];
    }
    X = to_matrix(x_donor, d, r);
    Q = polar(X);
    Sigma_donor = diag_post_multiply(Q_donor, sigma_donor);
    for (i in 1:k)
      b_donor[i] = Sigma_donor * z_donor[i];
  }
}
model {
  // priors
  for (j in 1:d)
    beta[j] ~ normal(0, 7);
  sigma ~ cauchy(0, 2.5);
  sigma_term ~ cauchy(0, 2.5);
  sigma_donor ~ cauchy(0, 2.5);
  //L ~ lkj_corr_cholesky(eta);
  //L_term ~ lkj_corr_cholesky(eta);
  //L_donor ~ lkj_corr_cholesky(eta);
  for (i in 1:n) {
    z[i] ~ std_normal();
    z_term[i] ~ std_normal();
  }
  for (i in 1:k)
    z_donor[i] ~ std_normal();
  x ~ std_normal();
  x_term ~ std_normal();
  x_donor ~ std_normal();
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
  //Cor = L * L';
  //Cor_term = L_term * L_term';
  //Cor_donor = L_donor * L_donor';
  Cor = cov2cor(Q, sigma);
  Cor_term = cov2cor(Q_term, sigma_term);
  Cor_donor = cov2cor(Q_donor, sigma_donor);
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
