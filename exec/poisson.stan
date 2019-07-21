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
}
data {
  int<lower=1> n; // num of cells
  int<lower=1> d; // num of markers
  int<lower=1> p; // num of explanatory variables (including intercept)
  int<lower=0> Y[d,n]; // observed cell counts
  int<lower=1> k; // number of donors
  int<lower=1,upper=k> donor[n]; // donor indicator
  int<lower=1,upper=p> term[n]; // condition indicator
  int<lower=1> r_cell; // rank of latent matrix
  int<lower=1> r_donor; // rank of latent matrix
}
transformed data {
  int<lower=1,upper=n> index_sorted[n] = sort_indices_asc(term);
  int<lower=1,upper=p> term_sorted[n] = term[index_sorted];
  int<lower=1,upper=n> n_ref_cond = rank(term_sorted, n);
  int<lower=1,upper=n> n_tar_cond = n-n_ref_cond;
  int<lower=1,upper=k> donor_sorted[n] = donor[index_sorted];
  int<lower=0> Y_sorted[d,n] = Y[,index_sorted];
  int<lower=0> y_sorted[d*n] = to_array_1d(Y_sorted);
}
parameters {
  matrix[d,p] beta; // fixed coefficients
  vector<lower=0>[r_cell] sigma; // random effects std
  vector<lower=0>[r_cell] sigma_term; // random effects std
  vector<lower=0>[r_donor] sigma_donor; // random effects std
  vector[n_ref_cond*r_cell] z; // random effects
  vector[n_tar_cond*r_cell] z_term; // random effects
  vector[k*r_donor] z_donor; // random effects
  vector[d*r_cell] x; // distribution on X for polar expansion
  vector[d*r_cell] x_term; // distribution on X for polar expansion
  vector[d*r_donor] x_donor; // distribution on X for polar expansion
}
transformed parameters {
  matrix[k,d] b_donor;
  matrix[d,r_cell] Q;
  matrix[r_cell,d] Sigma_t;
  matrix[d,r_cell] Q_term;
  matrix[r_cell,d] Sigma_term_t;
  matrix[d,r_donor] Q_donor;
  matrix[r_donor,d] Sigma_donor_t;
  {
    matrix[d,r_cell] X_ref;
    matrix[d,r_cell] X_tar;
    // reference level
    X_ref = to_matrix(x, d, r_cell);
    Q = polar(X_ref);
    Sigma_t = diag_post_multiply(Q, sigma)';
    // target level
    X_tar = to_matrix(x_term, d, r_cell);
    Q_term = polar(X_tar);
    Sigma_term_t = diag_post_multiply(Q_term, sigma_term)';
  }
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
  sigma_term ~ cauchy(0, 2.5);
  sigma_donor ~ cauchy(0, 2.5);
  z ~ std_normal();
  z_term ~ std_normal();
  z_donor ~ std_normal();
  x ~ std_normal();
  x_term ~ std_normal();
  x_donor ~ std_normal();
  {
    matrix[n,d] b;
    matrix[n_ref_cond,r_cell] Z_ref;
    matrix[n_ref_cond,d] b_ref;
    matrix[n_tar_cond,r_cell] Z_tar;
    matrix[n_tar_cond,d] b_tar;
    // reference level
    Z_ref = to_matrix(z, n_ref_cond, r_cell);
    b_ref = Z_ref * Sigma_t;
    // target level
    Z_tar = to_matrix(z_term, n_tar_cond, r_cell);
    b_tar = Z_tar * Sigma_term_t;
    // combine levels
    b = append_row(b_ref, b_tar);
    // likelihood
    y_sorted ~ poisson_log(
      to_vector(beta'[term_sorted,] + b + b_donor[donor_sorted,])
      );
  }
}
generated quantities {
  // correlation matrix
  matrix[d,d] Cor;
  matrix[d,d] Cor_term;
  matrix[d,d] Cor_donor;
  Cor = cov2cor(Q, sigma);
  Cor_term = cov2cor(Q_term, sigma_term);
  Cor_donor = cov2cor(Q_donor, sigma_donor);
}
