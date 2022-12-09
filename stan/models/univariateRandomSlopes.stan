data {
  int<lower=0>    N;                     // number of observations
  int<lower=1> N_ids;                    // number of individuals
  array[N] int<lower=1, upper=N_ids> id; // individual IDs
  array[N] real Y;                       // y_2
  array[N] real X;                       // y_1
  cov_matrix[N_ids] A;                   // GRM matrix
  vector[2] bets_prior_trait;
  vector[2] bets_prior_slope;
}
transformed data{
  matrix[N_ids, N_ids] LA;
  LA = cholesky_decompose(A);
}
parameters {
  real beta_0;                       // slope mean
  real mu_0;                         // y_2 trait mean
  real<lower=0> sigma_y;             // y_2 std
  real<lower=0> sigma_beta;          // beta std
  real<lower=0, upper = 1> h2_y;     // variance partition of additive effects
  real<lower=0, upper = 1> h2_beta;  // variance partition of additive effects
  vector[N_ids]  a_tilde;            // breeding values
  vector[N_ids]  beta_tilde;         // random slopes
  vector[N]      beta;               
}
transformed parameters {
    vector[N_ids] a;
    vector[N_ids] beta_add;
    a = sigma_y * sqrt(h2_y) * (LA * a_tilde);
    beta_add = sigma_beta * sqrt(h2_beta) * (LA * beta_tilde);
}
model {
    vector[N] mu;
    vector[N] mu_beta;
    
    for(n in 1:N){
      mu_beta[n] = beta_0 + beta_add[id[n]];
      mu[n] = mu_0 + beta[n] * X[n] + a[id[n]];
    }
    beta ~ normal(mu_beta, sigma_beta*sqrt(1 - h2_beta));
    Y ~ normal(mu, sigma_y*sqrt(1 - h2_y));

    beta_0 ~ normal(0, 1);
    mu_0 ~ normal(0, 1);
    sigma_y ~ normal(0, 1);
    h2_y ~ beta(beta_prior_trait[1], beta_prior_trait[2]);
    sigma_beta ~ normal(0, 1);
    h2_beta ~ beta(bets_prior_slope[1], bets_prior_slope[2]);
    a_tilde ~ normal(0, 1);
    beta_tilde ~ normal(0, 1);
}
