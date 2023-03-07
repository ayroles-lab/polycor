# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
cmdstanr::install_cmdstan()

pak::pkg_install(c("rmcelreath/rethinking", "bayesplot", "posterior", "ggplot2", "cowplot", "patchwork", "loo", "tictoc"))

library(rethinking)
library(cmdstanr)
library(ggplot2)
library(cowplot)
library(posterior)
library(bayesplot)
library(patchwork)
library(loo)
library(tictoc)
color_scheme_set("brightblue")

sim <- new.env()
source(here::here("simFunctions.R"), local=sim)
source(here::here("stan/R/loadModels.R"))


#### Random slopes ####
{
N_id  = 1000
N_rep = 5
N = N_id * N_rep

h2_y = 0.3
sigma_y = 1

G = sim$simulateKinship(N_id)

h2_beta = 0.6
sigma_beta = 1

beta_0 = 1
beta.add = sigma_beta * sqrt(h2_beta) * t(chol(G)) %*% rnorm(N_id)
beta.sparse = rep(0, N_id); beta.sparse[1] = 10; beta.sparse[10] = 5; beta.sparse[20] = 2.5; beta.sparse[30] = 1.25
beta = beta_0 + rep(beta.add, each = N_rep) + rep(beta.sparse, each = N_rep) + rnorm(N, sd = sigma_beta * sqrt(1-h2_beta))

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
fit_uniPenalizedRandSlopes <- uniPenalizedRandSlopes$sample(
  data = data_list_uni,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 3000,
  iter_sampling = 1000,
  max_treedepth = 11,
  adapt_delta = 0.9,
  refresh = 0 
)
fit_uniRandSlopes$summary(c("sigma_y", "sigma_beta", "h2_y", "h2_beta"))
png("test.png", height = 700, width = 600)
mcmc_recover_intervals(fit_uniPenalizedRandSlopes$draws(c("sigma_y", "sigma_beta", "h2_y", "h2_beta")), 
                                                  true = c(sigma_y, sigma_beta, h2_y, h2_beta)) + theme(legend.position = "none") + ggtitle("Penalized") + 
mcmc_recover_intervals(fit_uniRandSlopes$draws(c("sigma_y", "sigma_beta", "h2_y", "h2_beta")), 
                                                  true = c(sigma_y, sigma_beta, h2_y, h2_beta)) + theme(legend.position = "none") + ggtitle("GRM") +
  plot_layout(guides = 'collect')
dev.off()

a_blups = data.frame(
  true = as.vector(y.add), 
  GRM = colMeans(colMeans(fit_uniRandSlopes$draws("a")))) 
(corrs = cor(a_blups))
# png("test.png", height = 500, width = 1000)
# a_blups |> ggplot(aes(true, GRM)) + labs(x = "Simulated True Value", y = "Model Estimate") +
#            geom_point() + geom_abline(intercept = 0, slope = 1) + theme_minimal() + ggtitle("Breeding value BLUPs (GRM)") +
# a_blups |> ggplot(aes(true, iid)) + labs(x = "Simulated True Value", y = "Model Estimate") +
#            geom_point() + geom_abline(intercept = 0, slope = 1) + theme_minimal() + ggtitle("Breeding value BLUPs (iid)") 
# dev.off()

slope_blups = data.frame(
  true = as.vector(beta.add), 
  GRM = colMeans(colMeans(fit_uniRandSlopes$draws("beta_add"))), 
  HAL = colMeans(colMeans(fit_uniPenalizedRandSlopes$draws("beta_add")))) 
(corrs = cor(slope_blups))
png("test.png", height = 500, width = 1000)
# (mcmc_recover_hist(fit_uniRandSlopes$draws(c("sigma_y", "sigma_beta", "h2_y", "h2_beta")), 
#                                            true = c(sigma_y, sigma_beta, h2_y, h2_beta)) + 
#                                            ggtitle("GRM") +
# mcmc_recover_hist(fit_uniPenalizedRandSlopes$draws(c("sigma_y", "sigma_beta", "h2_y", "h2_beta")), 
#                                              true = c(sigma_y, sigma_beta, h2_y, h2_beta)) + 
#                                              ggtitle("iid") + 
#   plot_layout(guides = 'collect') )/ 
(slope_blups |> ggplot(aes(true, GRM)) + 
           labs(x = "Simulated True Value", y = "Model Estimate") +
           geom_point() + geom_abline(intercept = 0, slope = 1) + 
           annotate("text", x = -1, y = 2, 
           label = paste0("Corr = ", round(corrs["GRM", "true"], 3)), size = 6) +
           theme_minimal() + ggtitle("Random Slope BLUPs (GRM)") +
slope_blups |> ggplot(aes(true, HAL)) + 
           labs(x = "Simulated True Value", y = "Model Estimate") +
           geom_point() + geom_abline(intercept = 0, slope = 2) + 
           annotate("text", x = -1, y = 2, 
           label = paste0("Corr = ", round(corrs["HAL", "true"], 3)), size = 6) +
           theme_minimal() + ggtitle("Random Slope BLUPs (HAL)"))
# slope_blups |> ggplot(aes(iid, GRM)) + 
#            labs(x = "iid estimate", y = "GRM estimate") +
#            geom_point() + geom_abline(intercept = 0, slope = 1) + theme_minimal() + ggtitle("Estimate comparison")
dev.off()

slope_penalized = data.frame(
  true = as.vector(beta.sparse), 
  GRM = colMeans(colMeans(fit_uniRandSlopes$draws("beta_add"))), 
  shrinkage = colMeans(colMeans(fit_uniPenalizedRandSlopes$draws("shrinkage"))), 
  HAL = colMeans(colMeans(fit_uniPenalizedRandSlopes$draws("beta_sparse")))) |> 
  dplyr::mutate(sparse = shrinkage < 0.5)

png("test.png", height = 1000, width = 1000)
(ggplot(slope_penalized, aes(1:N_id, GRM)) + geom_point() + theme_minimal()) /
(ggplot(slope_penalized, aes(1:N_id, HAL)) + geom_point() + theme_minimal()) /
ggplot(slope_penalized, aes(1:N_id, shrinkage, col = sparse)) + geom_point(size= 2) + scale_color_manual(values = 2:1) + theme_minimal() + ylim(0, 1) 
dev.off()
