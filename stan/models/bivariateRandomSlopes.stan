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
  int K = 2;                             // number of traits
  LA = cholesky_decompose(A);
}
parameters {
  real beta_0;                       // slope mean
  real mu_x_0;                       // y_1 trait mean
  real mu_y_0;                       // y_2 trait mean

  real<lower=0> sigma_beta;          // beta sd
  real<lower=0, upper = 1> h2_beta;  // variance partition of additive effects
  vector[N_ids]  beta_tilde;         // random slopes precursor
  vector[N]      beta;     

 
  vector<lower=0>[K] L_sigma;        // total trait sd
  vector<lower=0, upper = 1>[K] h2;  // variance partition of both traits
  matrix[N_ids, K] a_tilde;          // breeding values precursor
  cholesky_factor_corr[K] L_Omega_G; // G matrix         
}
transformed parameters {
  vector[N_ids] beta_add;        // random slopes
  matrix[N_ids, K] a;            // breeding values
  vector<lower=0>[K] L_sigma_G;  // trait genetic sd
  vector<lower=0>[K] L_sigma_R;  // trait residual sd

  for(k in 1:K){
    L_sigma_G[k] = L_sigma[k] * sqrt(h2[k]);
    L_sigma_R[k] = L_sigma[k] * sqrt(1-h2[k]);
  }
  a = (LA * a_tilde) * diag_pre_multiply(L_sigma_G, L_Omega_G)'; // a ~ N(0, A x G)
  beta_add = sigma_beta * sqrt(h2_beta) * (LA * beta_tilde);
}
model {
    vector[N] mu_x;
    vector[N] mu_y;
    vector[N] mu_beta;
    
    // per observation mean parameters
    for(n in 1:N){
      mu_beta[n] = beta_0 + beta_add[id[n]];
      mu_x[n] = mu_x_0 + a[id[n], 1];
      mu_y[n] = mu_y_0 + a[id[n], 2] + beta[n] * X[n];
    }

    // residual variation likelihoods and priors
    beta ~ normal(mu_beta, sigma_beta*sqrt(1 - h2_beta));
    X ~ normal(mu_x, L_sigma_R[1]);
    Y ~ normal(mu_y, L_sigma_R[2]);

    // Intercept priors
    beta_0 ~ normal(0, 1);
    mu_x_0 ~ normal(0, 1);
    mu_y_0 ~ normal(0, 1);

    // sd priors
    sigma_beta ~ normal(0, 1);
    L_sigma ~ normal(0, 1);

    // variance partition priors
    h2_beta ~ beta(2, 4);
    h2 ~ beta(2, 2);

    // random effect precursors
    beta_tilde ~ normal(0, 1);
    to_vector(a_tilde) ~ normal(0, 1);
    
    // trait genetic correlation prior
    L_Omega_G ~ lkj_corr_cholesky(4);
}
generated quantities {
    corr_matrix[K] corrG;
    real rho;

    corrG = multiply_lower_tri_self_transpose(L_Omega_G);
    rho = corrG[1,2];
}
