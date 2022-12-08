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

#### Fixed slopes ####

N_id  = 250
N_rep = 5
N = N_id * N_rep

beta_0 = 1
h2_y = 0.4
sigma_y = 1

G = simulateKinship(N_id)

x.add = rnorm(N_id)
x = rep(x.add, each = N_rep) + rnorm(N)

y.add = sigma_y * sqrt(h2_y) * t(chol(G)) %*% rnorm(N_id)
y = rep(y.add, each = N_rep) + beta_0 * x + rnorm(N, sd = sigma_y * sqrt(1-h2_y))

mod <- cmdstan_model(here::here("stan/models/univariate.stan"))
# data {
#   int<lower=0>    N;                     // number of observations
#   int<lower=1> N_ids;                    // number of individuals
#   array[N] int<lower=1, upper=N_ids> id; // individual IDs
#   array[N] real Y;                       // y_2
#   array[N] real X;                       // y_1
#   cov_matrix[N_ids] A;                   // GRM matrix
# }
data_list = list(N = N, 
                 N_ids = N_id,
                 id = rep(1:N_id, each = N_rep),
                 Y = y,
                 X = x,
                 A = as.matrix(G))

fit <- mod$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 500 # print update every 500 iters
)
fit$summary()
png("test.png", height = 900, width = 1000)
mcmc_recover_hist(fit$draws("mu_0"), 0) +
mcmc_recover_hist(fit$draws("beta_0"), 1)
dev.off()
png("test.png", height = 900, width = 2000)
mcmc_recover_intervals(fit$draws("a"), as.vector(y.add))
dev.off()
png("test.png", height = 900, width = 1000)
data.frame(estimate = colMeans(colMeans(fit$draws("a"))), 
           true = as.vector(y.add)) %>% ggplot(aes(true, estimate)) + geom_point() + geom_abline(intercept = 0, slope = 1) + theme_minimal()
dev.off()


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

mod <- cmdstan_model(here::here("stan/models/univariateRandomSlopes.stan"))
# data {
#   int<lower=0>    N;                     // number of observations
#   int<lower=1> N_ids;                    // number of individuals
#   array[N] int<lower=1, upper=N_ids> id; // individual IDs
#   array[N] real Y;                       // y_2
#   array[N] real X;                       // y_1
#   cov_matrix[N_ids] A;                   // GRM matrix
# }
data_list = list(N = N, 
                 N_ids = N_id,
                 id = rep(1:N_id, each = N_rep),
                 Y = y,
                 X = x,
                 A = as.matrix(G))
fit <- mod$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 500 # print update every 500 iters
)
fit$summary()
png("test.png", height = 900, width = 1000)
(mcmc_recover_hist(fit$draws("mu_0"), 0) +
mcmc_recover_hist(fit$draws("beta_0"), 1)) / 
(mcmc_recover_hist(fit$draws("h2_y"), h2_y) +
mcmc_recover_hist(fit$draws("h2_beta"), h2_beta))
dev.off()
png("test.png", height = 900, width = 2000)
mcmc_recover_intervals(fit$draws("a"), as.vector(y.add))
dev.off()
png("test.png", height = 450, width = 1000)
data.frame(estimate = colMeans(colMeans(fit$draws("a"))), 
           true = as.vector(y.add)) %>% ggplot(aes(true, estimate)) + labs(x = "Simulated True Value", y = "Model Estimate") +
           geom_point() + geom_abline(intercept = 0, slope = 1) + theme_minimal() + ggtitle("Breeding value BLUPs (g.add)") +
data.frame(estimate = colMeans(colMeans(fit$draws("beta_add"))), 
           true = as.vector(beta.add)) %>% ggplot(aes(true, estimate)) + labs(x = "Simulated True Value", y = "Model Estimate") +
           geom_point() + geom_abline(intercept = 0, slope = 1) + theme_minimal() + ggtitle("Random Slope BLUPs (beta.add)")
dev.off()

