# in simulation of rep. observations:
# - the genetic part is simulated once (common across rep)
# - the residual part is simulated for each rep.

### inc 
library(tidyverse)

library(Matrix)
library(MASS)

library(nlme)
library(lme4)

library(cowplot)

# change to your path here, e.g. based on `nodemame`
path <- switch(Sys.info()[['nodename']],
  "~/git/hemostat/polycor/" )
  
src <- file.path(path, "relmer.R")
source(src)

### par
testing <- FALSE

sim_rand_gen <- TRUE

vals_n <- c(100, 500, 1e3)
vals_h <- c(0.1, 0.5, 0.9)
vals_rep <- c(2, 5)
vals_sim <- 1:5

if(testing) {
 vals_n <- c(10, 50)
}

### local functions
sim_ids <- function(N) paste0("id", 1:N)

sim_kin <- function(h, N, L = 5)
{
  stopifnot(!(N %% L)) # family size = L
  Nl <- N / L
  
  # kinship of a single family
  K0 <- matrix(0.5, L, L)
  diag(K0) <- 1
  K0[1, 2] <- 0
  K0[2, 1] <- 0
  
  K <- kronecker(Diagonal(Nl), K0)
  
  ids <- sim_ids(N) 
  rownames(K) <- ids
  colnames(K) <- ids
  
  K
}

sim_univar <- function(h, N, rep = 1, L = 5, scale = TRUE, sim_rand_gen = FALSE) 
{
  stopifnot(!(N %% L)) # family size = L
  Nl <- N / L
  
  K <- sim_kin(h, N, L)
  
  # V = [K; I]
  V_gen <- h*K 
  V_res <- (1-h)*Diagonal(N)
  
  if(sim_rand_gen) {
    y_gen <- mvrnorm(1, rep(0, N), V_gen)
  }
  
  lapply(seq(1, rep), function(r) {
    y_res <- mvrnorm(1, rep(0, N), V_res)
    
    if(sim_rand_gen) {
      y_gen <- mvrnorm(1, rep(0, N), V_gen)
    }
    
    y <- y_gen + y_res
    if(scale) {
      y <- scale(y) %>% as.numeric
    }
    
    ids <- sim_ids(N)
    reps <- paste0("rep", r)
    data_frame(id = ids, rid = ids, y = y, rep = reps)
  }) %>% 
  bind_rows
}

### run 
grid <- expand.grid(n = as.integer(vals_n), h = vals_h, nrep = vals_rep, sim = vals_sim) %>%
  as_tibble

out <- lapply(seq(nrow(grid)), function(i)
{
  n <- grid$n[i]
  h <- grid$h[i]
  nrep <- grid$nrep[i]
  sim <- grid$sim[i]
  
  cat(" -", i, "/", nrow(grid), "\n")
  print(grid[i, ])
  
  kin <- sim_kin(h, n/nrep)
  dat <- sim_univar(h, n/nrep, nrep, sim_rand_gen = sim_rand_gen)

  mod <- relmer(y ~ 1 - rep + (1| id) + (1|rid), data = dat, relmat = list(id = kin), calc.derivs = FALSE)
  vf <- VarProp(mod)
  h_hat <- with(vf, prop[grp == "id"])
  e_hat <- with(vf, prop[grp == "rid"] + prop[grp == "Residual"])
  cove_hat <- with(vf, prop[grp == "rid"])
  
  data_frame(sim = sim, n = n, h = h, e = 1-h, cove = e,
    h_hat = h_hat, e_hat = e_hat, cove_hat = cove_hat)
})

tab <- bind_rows(out)

### save
saveRDS(tab, "univar_rep_h.rds")

### plot
p1 <- ggplot(tab, aes(h, h_hat, group = as.factor(h))) + geom_boxplot() + facet_wrap(~ n) + geom_abline(linetype = 3) + coord_equal()
p2 <- ggplot(tab, aes(e, e_hat, group = as.factor(e))) + geom_boxplot() + facet_wrap(~ n) + geom_abline(linetype = 3) + coord_equal()
p3 <- ggplot(tab, aes(cove, cove_hat, group = as.factor(cove))) + geom_boxplot() + facet_wrap(~ n) + geom_abline(linetype = 3) + coord_equal()

plot_grid(p1, p2, p3, ncol = 1, align = "h")

