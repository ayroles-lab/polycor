### inc 
library(tidyverse)
library(cowplot)

library(Matrix)
library(MASS)

library(nlme)
library(lme4)

# change to your path here, e.g. based on `nodemame`
path <- switch(Sys.info()[['nodename']],
  "~/git/hemostat/polycor/" )
  
src <- file.path(path, "relmer.R")
source(src)

### par
testing <- FALSE

nrep <- 4
vals_n <- c(100, 500)
vals_rhoe <- c(0.1, 0.5, 0.9)
vals_sim <- 1:5

if(testing) {
  nrep <- 2
  vals_n <- c(10, 50)
}

e1 <- 0.5
e2 <- 0.5

### local functions
sim_ids <- function(N) paste0("id", 1:N)
sim_obs <- function(N) paste0("obs", 1:N)

sim_bivar_env <- function(e1, e2, re, N, rep = 1) 
{
  # V = [I]
  e12 <- re*sqrt(e1)*sqrt(e2)
  Sigma_resid <- matrix(c(e1, e12, e12, e2), 2, 2)
  V_resid <- kronecker(Sigma_resid, Diagonal(N))

  V <- V_resid
  
  lapply(seq(1, rep), function(r) {
    y12 <- mvrnorm(1, rep(0, 2*N), V)
    y1 <- y12[seq(1, N)]
    y2 <- y12[N + seq(1, N)]
    
    ids <- sim_ids(N)
    reps <- paste0("rep", r)
    data_frame(id = ids, rid = ids, y1 = y1, y2 = y2, rep = reps)
  }) %>% 
  bind_rows %>%
  mutate(obs = sim_obs(n()))
}

### run 1: no rep.
grid <- expand.grid(n = as.integer(vals_n), rhoe = vals_rhoe, sim = vals_sim) %>%
  as_tibble

out <- lapply(seq(nrow(grid)), function(i)
{
  n <- grid$n[i]
  rhoe <- grid$rhoe[i]
  sim <- grid$sim[i]
  
  cat(" -", i, "/", nrow(grid), "\n")
  print(grid[i, ])
  
  dat <- sim_bivar_env(e1, e2, rhoe, n/nrep, nrep)
  bdat <- gather(dat, tname, tvalue, y1, y2) %>%
    mutate(tindex = (as.integer(as.factor(tname)) - 1))

  mod <- relmer(tvalue ~ tname - 1 + (0 + tname | id) + (0 + tname | obs), data = bdat, weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE)
  vf <- VarCorr(mod) %>% as.data.frame %>% filter(grp == "obs" & var1 == "tnamey1" & var2 == "tnamey2")
  rhoe_hat <- vf$sdcor
  
  vf <- VarCorr(mod) %>% as.data.frame %>% filter(grp == "id" & var1 == "tnamey1" & var2 == "tnamey2")
  rhoe0_hat <- vf$sdcor

  mod2 <- relmer(tvalue ~ tname - 1 + (0 + dummy(tname, "y1") | id) + (0 + dummy(tname, "y2") | id) + (0 + tname | obs), data = bdat, weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE)

  vf <- VarCorr(mod2) %>% as.data.frame %>% filter(grp == "obs" & var1 == "tnamey1" & var2 == "tnamey2")
  rhoe2_hat <- vf$sdcor
 
  data_frame(sim = sim, n = n/nrep, nrep = nrep, rhoe = rhoe, rhoe_hat = rhoe_hat, rhoe0_hat = rhoe0_hat, rhoe2_hat = rhoe2_hat, e1 = e1, e2 = e2)
})

tab <- bind_rows(out)

### save
#saveRDS(tab, "env_rep_rhoe.rds")

# plot
p1 <- ggplot(tab, aes(rhoe, rhoe_hat, group = as.factor(rhoe))) + geom_boxplot() + facet_wrap(~ n, nrow = 1) + geom_abline(linetype = 3) + coord_equal() + theme_minimal()

p2 <- ggplot(tab, aes(rhoe, rhoe0_hat, group = as.factor(rhoe))) + geom_boxplot() + facet_wrap(~ n, nrow = 1) + geom_abline(linetype = 3) + coord_equal() + theme_minimal()

p3 <- ggplot(tab, aes(rhoe, rhoe2_hat, group = as.factor(rhoe))) + geom_boxplot() + facet_wrap(~ n, nrow = 1) + geom_abline(linetype = 3) + coord_equal() + theme_minimal()


