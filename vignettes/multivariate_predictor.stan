/*
 * Logistic regression with index notation
 * Author: Christof Seiler
 */
data {
  int<lower=1> p;
  int<lower=0> n;
  int<lower=0,upper=1> y[n];
  row_vector[p] x[n];
}
parameters {
  vector[p] beta;
}
model {
  for (i in 1:n)
    y[i] ~ bernoulli_logit(x[i] * beta);
}
