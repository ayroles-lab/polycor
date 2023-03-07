# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
cmdstanr::install_cmdstan()



pak::pkg_install(c("rmcelreath/rethinking", "bayesplot", "posterior", "ggplot2", "cowplot", "patchwork", "loo", "tictoc", "rio"))

library(rethinking)
library(cmdstanr)
library(ggplot2)
library(cowplot)
library(posterior)
library(bayesplot)
library(patchwork)
library(loo)
library(tictoc)
library(rio)
color_scheme_set("brightblue")

sim <- new.env()
source(here::here("simFunctions.R"), local=sim)
source(here::here("stan/R/loadModels.R"))


#### Random slopes ####
{
N_id  = 250
N_rep = 10
N = N_id * N_rep

h2_y = 0.3
sigma_y = 1

G = sim$simulateKinship(N_id)

h2_beta = 0.6
sigma_beta = 1

beta_0 = 1
beta.add = sigma_beta * sqrt(h2_beta) * t(chol(G)) %*% rnorm(N_id)
beta = beta_0 + rep(beta.add, each = N_rep) + rnorm(N, sd = sigma_beta * sqrt(1-h2_beta))

x.add = rnorm(N_id)
x = rep(x.add, each = N_rep) + rnorm(N)

y.add = sigma_y * sqrt(h2_y) * t(chol(G)) %*% rnorm(N_id)
y = rep(y.add, each = N_rep) + beta * x + rnorm(N, sd = sigma_y * sqrt(1-h2_y))

data_list_uni = list(N = N, 
                     N_ids = N_id,
                     id = rep(1:N_id, each = N_rep),
                     Y = y - mean(y),
                     X = x,
                     A = as.matrix(G),
                     h2_prior_trait = c(3, 4),
                     h2_prior_slope = c(3, 4))
}

fit_uniRandSlopes <- uniRandSlopes$sample(
  data = data_list_uni,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 0
)
fit_uniRandSlopesiid <- uniRandSlopesiid$sample(
  data = data_list_uni,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 0 
)
fit_uniRandSlopes$summary(c("sigma_y", "sigma_beta", "h2_y", "h2_beta"))
png("test.png", height = 700, width = 600)
mcmc_recover_intervals(fit_uniRandSlopesiid$draws(c("sigma_y", "sigma_beta", "h2_y", "h2_beta")), 
                                                  true = c(sigma_y, sigma_beta, h2_y, h2_beta)) + theme(legend.position = "none") + ggtitle("iid") + 
mcmc_recover_intervals(fit_uniRandSlopes$draws(c("sigma_y", "sigma_beta", "h2_y", "h2_beta")), 
                                                  true = c(sigma_y, sigma_beta, h2_y, h2_beta)) + theme(legend.position = "none") + ggtitle("GRM") +
  plot_layout(guides = 'collect')
dev.off()

a_blups = data.frame(
  true = as.vector(y.add), 
  GRM = colMeans(colMeans(fit_uniRandSlopes$draws("a"))), 
  iid = colMeans(colMeans(fit_uniRandSlopesiid$draws("a")))) 
(corrs = cor(a_blups))
png("test.png", height = 500, width = 1000)
a_blups |> ggplot(aes(true, GRM)) + labs(x = "Simulated True Value", y = "Model Estimate") +
           geom_point() + geom_abline(intercept = 0, slope = 1) + theme_minimal() + ggtitle("Breeding value BLUPs (GRM)") +
a_blups |> ggplot(aes(true, iid)) + labs(x = "Simulated True Value", y = "Model Estimate") +
           geom_point() + geom_abline(intercept = 0, slope = 1) + theme_minimal() + ggtitle("Breeding value BLUPs (iid)") 
dev.off()

slope_blups = data.frame(
  true = as.vector(beta.add), 
  GRM = colMeans(colMeans(fit_uniRandSlopes$draws("beta_add"))), 
  iid = colMeans(colMeans(fit_uniRandSlopesiid$draws("beta_add")))) 
(corrs = cor(slope_blups))
png("slopesGRM_iid.png", height = 1000, width = 1000)
(mcmc_recover_hist(fit_uniRandSlopes$draws(c("sigma_y", "sigma_beta", "h2_y", "h2_beta")), 
                                           true = c(sigma_y, sigma_beta, h2_y, h2_beta)) + 
                                           ggtitle("GRM") +
mcmc_recover_hist(fit_uniRandSlopesiid$draws(c("sigma_y", "sigma_beta", "h2_y", "h2_beta")), 
                                             true = c(sigma_y, sigma_beta, h2_y, h2_beta)) + 
                                             ggtitle("iid") + 
  plot_layout(guides = 'collect') )/ 
(slope_blups |> ggplot(aes(true, GRM)) + 
           labs(x = "Simulated True Value", y = "Model Estimate") +
           geom_point() + geom_abline(intercept = 0, slope = 1) + 
           annotate("text", x = -1, y = 2, 
           label = paste0("Corr = ", round(corrs["GRM", "true"], 3)), size = 6) +
           theme_minimal() + ggtitle("Random Slope BLUPs (GRM)") +
slope_blups |> ggplot(aes(true, iid)) + 
           labs(x = "Simulated True Value", y = "Model Estimate") +
           geom_point() + geom_abline(intercept = 0, slope = 2) + 
           annotate("text", x = -1, y = 2, 
           label = paste0("Corr = ", round(corrs["iid", "true"], 3)), size = 6) +
           theme_minimal() + ggtitle("Random Slope BLUPs (iid)"))
# slope_blups |> ggplot(aes(iid, GRM)) + 
#            labs(x = "iid estimate", y = "GRM estimate") +
#            geom_point() + geom_abline(intercept = 0, slope = 1) + theme_minimal() + ggtitle("Estimate comparison")
dev.off()

#### Random slopes and genetic correlations ####
{
N_id  = 250
N_rep = 5
N = N_id * N_rep

h2_x = 0.4
sigma_x = 1

h2_y = 0.4
sigma_y = 1

h2_beta = 0.8
sigma_beta = 1

G = sim$simulateKinship(N_id)

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

# data {
#   int<lower=0>    N;                     // number of observations
#   int<lower=1> N_ids;                    // number of individuals
#   array[N] int<lower=1, upper=N_ids> id; // individual IDs
#   array[N] real Y;                       // y_2
#   array[N] real X;                       // y_1
#   cov_matrix[N_ids] A;                   // GRM matrix
#   vector[2] h2_prior_trait;
#   vector[2] h2typos_prior_slope;
# }
data_list_bi = list(N = N, 
                    N_ids = N_id,
                    id = rep(1:N_id, each = N_rep),
                    Y = y,
                    X = x,
                    A = as.matrix(G), 
                    h2_prior_trait = c(3, 4),
                    h2_prior_slope = c(3, 4))
}
fit_biRandSlopes <- biRandSlopes$sample(
  data = data_list_bi,
  seed = 123,
  chains = 4,
  iter_warmup = 2000,
  iter_sampling = 1000,
  parallel_chains = 4,
  adapt_delta = 0.9, max_treedepth = 12,
  refresh = 500 # print update every 500 iters
)
fit_biRandSlopes$summary(c("L_sigma", "h2", "h2_beta",  "rho")) 
png("test.png", height = 900, width = 1000)
mcmc_recover_intervals(fit_biRandSlopes$draws(c("h2", "h2_beta",  "rho")), 
                                              true = c(h2_x, h2_y, h2_beta, rhog)) + 
                                              theme(legend.position = "none") + 
                                              ggtitle("GRM - Pleiotropic") + 
  plot_layout(guides = 'collect')
dev.off()
  
png("test.png", height = 450, width = 1000)
data.frame(estimate = colMeans(colMeans(fit_biRandSlopes$draws("a"))), 
           true = c(as.vector(x.add), as.vector(y.add)),
           trait = rep(c("x", "y"), each = N_id)) |> 
           ggplot(aes(true, estimate, color = trait)) + 
           labs(x = "Simulated True Value", y = "Model Estimate") +
           geom_point() + geom_abline(intercept = 0, slope = 1) + theme_minimal() + ggtitle("Breeding value BLUPs (g.add)") +
data.frame(estimate = colMeans(colMeans(fit_biRandSlopes$draws("beta_add"))), 
           true = as.vector(beta.add)) |> ggplot(aes(true, estimate)) + labs(x = "Simulated True Value", y = "Model Estimate") +
           geom_point() + geom_abline(intercept = 0, slope = 1) + theme_minimal() + ggtitle("Random Slope BLUPs (beta.add)")
dev.off()

data.frame(estimate = colMeans(colMeans(fit_bi$draws("beta_add"))), 
           true = as.vector(beta.add)) |> cor()


x = list(
function()
fit_bi <<- biRandSlopes$sample(
  data = data_list_uni,
  seed = 123,
  chains = 4,
  iter_warmup = 2000,
  iter_sampling = 1000,
  parallel_chains = 4,
  max_treedepth = 11,
  refresh = 0 
),
function()
fit_bi_true <<- biRandSlopes$sample(
  data = data_list_bi,
  seed = 123,
  chains = 4,
  iter_warmup = 2000,
  iter_sampling = 1000,
  parallel_chains = 4,
  max_treedepth = 11,
  refresh = 0 
))
library(doMC)
library(parallel)
registerDoMC(2)
mclapply(x, function(x) x())
fit_bi$summary("rho")      # no correlation 
fit_bi_true$summary("rho") # correlation of 0.5
