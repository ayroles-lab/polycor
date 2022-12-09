# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
cmdstanr::install_cmdstan()

pak::pkg_install(c("rmcelreath/rethinking", "bayesplot", "posterior", "ggplot2", "cowplot", "patchwork"))

library(rethinking)
library(cmdstanr)
library(ggplot2)
library(cowplot)
library(posterior)
library(bayesplot)
library(patchwork)
color_scheme_set("brightblue")

source(here::here("simFunctions.R"))

#### Random slopes ####

N_id  = 250
N_rep = 5
N = N_id * N_rep

h2_y = 0.4
sigma_y = 1

G = simulateKinship(N_id)

h2_beta = 0.8
sigma_beta = 1

beta_0 = 1
beta.add = sigma_beta * sqrt(h2_beta) * t(chol(G)) %*% rnorm(N_id)
beta = beta_0 + rep(beta.add, each = N_rep) + rnorm(N, sd = sigma_beta * sqrt(1-h2_beta))

x.add = rnorm(N_id)
x = rep(x.add, each = N_rep) + rnorm(N)

y.add = sigma_y * sqrt(h2_y) * t(chol(G)) %*% rnorm(N_id)
y = rep(y.add, each = N_rep) + beta * x + rnorm(N, sd = sigma_y * sqrt(1-h2_y))

mod_uni <- cmdstan_model(here::here("stan/models/univariateRandomSlopes.stan"))
# data {
#   int<lower=0>    N;                     // number of observations
#   int<lower=1> N_ids;                    // number of individuals
#   array[N] int<lower=1, upper=N_ids> id; // individual IDs
#   array[N] real Y;                       // y_2
#   array[N] real X;                       // y_1
#   cov_matrix[N_ids] A;                   // GRM matrix
# }
data_list_uni = list(N = N, 
                     N_ids = N_id,
                     id = rep(1:N_id, each = N_rep),
                     Y = y - mean(y),
                     X = x,
                     A = as.matrix(G))
fit_uni <- mod_uni$sample(
  data = data_list_uni,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 500 # print update every 500 iters
)
fit_uni$summary()
png("test.png", height = 900, width = 1000)
(mcmc_recover_hist(fit_uni$draws("mu_0"), 0) +
mcmc_recover_hist(fit_uni$draws("beta_0"), 1)) / 
(mcmc_recover_hist(fit_uni$draws("h2_y"), h2_y) +
mcmc_recover_hist(fit_uni$draws("h2_beta"), h2_beta)) +
  plot_layout(guides = 'collect')
dev.off()

png("test.png", height = 450, width = 1000)
data.frame(estimate = colMeans(colMeans(fit_uni$draws("a"))), 
           true = as.vector(y.add)) |> ggplot(aes(true, estimate)) + labs(x = "Simulated True Value", y = "Model Estimate") +
           geom_point() + geom_abline(intercept = 0, slope = 1) + theme_minimal() + ggtitle("Breeding value BLUPs (g.add)") +
data.frame(estimate = colMeans(colMeans(fit_uni$draws("beta_add"))), 
           true = as.vector(beta.add)) |> ggplot(aes(true, estimate)) + labs(x = "Simulated True Value", y = "Model Estimate") +
           geom_point() + geom_abline(intercept = 0, slope = 1) + theme_minimal() + ggtitle("Random Slope BLUPs (beta.add)")
dev.off()

#### Random slopes and genetic correlations ####

N_id  = 250
N_rep = 5
N = N_id * N_rep

h2_x = 0.4
sigma_x = 1

h2_y = 0.2
sigma_y = 1

h2_beta = 0.8
sigma_beta = 1

G = simulateKinship(N_id)

rhog = 0.5
varG = c(sigma_x * h2_x, sigma_y * h2_y)
Sigma_xy = sqrt(varG) %*% t(sqrt(varG)) * matrix(c(1, 0.5, 0.5, 1), 2) 
g_tilde = MASS::mvrnorm(N_id, c(0, 0), Sigma = diag(2))
g.add = as.matrix(t(chol(G)) %*% g_tilde %*% chol(Sigma_xy))

beta_0 = 1
beta.add = sigma_beta * sqrt(h2_beta) * t(chol(G)) %*% rnorm(N_id)
beta = beta_0 + rep(beta.add, each = N_rep) + rnorm(N, sd = sigma_beta * sqrt(1-h2_beta))

x.add = g.add[,1]
x = rep(x.add, each = N_rep) + rnorm(N, sd = sigma_x * sqrt(1-h2_x))

y.add = g.add[,2]
y = rep(y.add, each = N_rep) + rnorm(N, sd = sigma_y * sqrt(1-h2_y)) + beta * x

mod_bi <- cmdstan_model(here::here("stan/models/bivariateRandomSlopes.stan"))
# data {
#   int<lower=0>    N;                     // number of observations
#   int<lower=1> N_ids;                    // number of individuals
#   array[N] int<lower=1, upper=N_ids> id; // individual IDs
#   array[N] real Y;                       // y_2
#   array[N] real X;                       // y_1
#   cov_matrix[N_ids] A;                   // GRM matrix
# }
data_list_bi = list(N = N, 
                    N_ids = N_id,
                    id = rep(1:N_id, each = N_rep),
                    Y = y - mean(y),
                    X = x,
                    A = as.matrix(G))
fit_bi <- mod_bi$sample(
  data = data_list_bi,
  seed = 123,
  chains = 8,
  iter_warmup = 2000,
  iter_sampling = 2000,
  parallel_chains = 8,
  adapt_delta = 0.99, max_treedepth = 13,
  refresh = 500 # print update every 500 iters
)
fit_bi$summary()
png("test.png", height = 900, width = 1000)
(mcmc_recover_hist(fit_bi$draws("mu_x_0"), 0) +
 mcmc_recover_hist(fit_bi$draws("mu_y_0"), 0) + 
 mcmc_recover_hist(fit_bi$draws("beta_0"), 1)) / 
(mcmc_recover_hist(fit_bi$draws("h2"), c(h2_x, h2_y)) +
 mcmc_recover_hist(fit_bi$draws("h2_beta"), h2_beta) + 
 mcmc_recover_hist(fit_bi$draws("rho"), rhog)) +
  plot_layout(guides = 'collect')
dev.off()

{png("test.png", height = 450, width = 1000)
data.frame(estimate = colMeans(colMeans(fit_bi$draws("a"))), 
           true = c(as.vector(x.add), as.vector(y.add)), 
           trait = rep(c("x", "y"), each = N_id)) |> 
           ggplot(aes(true, estimate, color = trait)) + 
           labs(x = "Simulated True Value", y = "Model Estimate") +
           geom_point() + geom_abline(intercept = 0, slope = 1) + theme_minimal() + ggtitle("Breeding value BLUPs (g.add)") +
data.frame(estimate = colMeans(colMeans(fit_bi$draws("beta_add"))), 
           true = as.vector(beta.add)) |> ggplot(aes(true, estimate)) + labs(x = "Simulated True Value", y = "Model Estimate") +
           geom_point() + geom_abline(intercept = 0, slope = 1) + theme_minimal() + ggtitle("Random Slope BLUPs (beta.add)")
dev.off()}
