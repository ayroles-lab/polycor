data {
  int<lower=0>    N;                     // number of observations
  int<lower=1> N_ids;                    // number of individuals
  array[N] int<lower=1, upper=N_ids> id; // individual IDs
  array[N] real Y;                       // y_2
  array[N] real X;                       // y_1
  cov_matrix[N_ids] A;                   // GRM matrix
}
transformed data{
  matrix[N_ids, N_ids] LA;
  LA = cholesky_decompose(A);
}
parameters {
  real beta_0;                   // slope mean
  real mu_0;                     // y_2 trait mean
  real<lower=0> sigma_y;         // y_2 std
  real<lower=0, upper = 1> h2_y; // variance partition of additive effects
  vector[N_ids]  a_tilde;            // breeding values
}
transformed parameters {
    vector[N_ids] a;
    a = sigma_y * sqrt(h2_y) * (LA * a_tilde);
}
model {
    vector[N] mu;
    
    for(n in 1:N)
      mu[n] = mu_0 + beta_0 * X[n] + a[id[n]];

    Y ~ normal(mu, sigma_y*sqrt(1 - h2_y));

    beta_0 ~ normal(0, 1);
    mu_0 ~ normal(0, 1);
    sigma_y ~ normal(0, 1);
    h2_y ~ beta(2, 4);
    a_tilde ~ normal(0, 1);
}
