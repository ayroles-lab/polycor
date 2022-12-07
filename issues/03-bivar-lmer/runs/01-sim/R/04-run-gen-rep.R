### inc 
library(tidyverse)

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
vals_n <- 4e3 # c(100, 500, 1e3)
vals_rhog <- c(0.1, 0.5, 0.9)
vals_sim <- 1:5

if(testing) {
  nrep <- 2
  vals_n <- c(10, 50)
}

h1 <- 0.75
h2 <- 0.75
rhoe <- 0.1

### local functions
sim_ids <- function(N) paste0("id", 1:N)
sim_obs <- function(N) paste0("obs", 1:N)

sim_kin <- function(h1, h2, rg, re, N, L = 5)
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

sim_bivar <- function(h1, h2, rg, re, N, rep = 1, L = 5) 
{
  stopifnot(!(N %% L)) # family size = L
  Nl <- N / L
  
  K <- sim_kin(h1, h2, rg, re, N, L)
  
  # V = [K; I]
  h12 <- rg*sqrt(h1)*sqrt(h2)
  Sigma_gen <- matrix(c(h1, h12, h12, h2), 2, 2)
  V_gen <- kronecker(Sigma_gen, K)

  e12 <- re*sqrt(1-h2)*sqrt(1-h2)
  Sigma_resid <- matrix(c(1-h1, e12, e12, 1-h2), 2, 2)
  V_resid <- kronecker(Sigma_resid, Diagonal(N))

  V <- V_gen + V_resid
  
  lapply(seq(1, rep), function(r) {
    y12 <- mvrnorm(1, rep(0, 2*N), V)
    y1 <- y12[seq(1, N)]
    y2 <- y12[N + seq(1, N)]
    
    ids <- sim_ids(N)
    reps <- paste0("rep", r)
    data_frame(id = ids, rid = ids, y1 = y1, y2 = y2, n/nrep, nrep)
  }) %>% 
  bind_rows %>%
  mutate(obs = sim_obs(n()), robs = obs)
}

### run: no rep.
grid <- expand.grid(n = as.integer(vals_n), rhog = vals_rhog, sim = vals_sim) %>%
  as_tibble

out <- lapply(seq(nrow(grid)), function(i)
{
  n <- grid$n[i]
  rhog <- grid$rhog[i]
  sim <- grid$sim[i]
  
  cat(" -", i, "/", nrow(grid), "\n")
  print(grid[i, ])
  
  dat <- sim_bivar(h1, h2, rhog, rhoe, n/nrep, nrep)
  bdat <- gather(dat, tname, tvalue, y1, y2) %>%
    mutate(tindex = (as.integer(as.factor(tname)) - 1))

  kin <- sim_kin(h1, h2, rhog, rhoe, n/nrep)
  kinobs <- kin[dat$id, dat$id] # kinship for observations with reptetive IDs
  rownames(kinobs) <- colnames(kinobs) <- dat$obs

  mod <- relmer(tvalue ~ tname - 1 + (0 + tname | id) + (0 + tname | obs) + (0 + tname | rid) + (0 + tname | robs), data = bdat, weights = rep(1e10, nrow(bdat)), relmat = list(id = kin, obs = kinobs), calc.derivs = FALSE, verbose = 2)
  vf <- VarCorr(mod) %>% as.data.frame %>% filter(grp == "obs" & var1 == "tnamey1" & var2 == "tnamey2")
  rhog_hat <- vf$sdcor
  
  data_frame(sim = sim, n = n/nrep, nrep = nrep, rhog = rhog, rhog_hat = rhog_hat, h1 = h1, h2 = h2, rhoe = rhoe)
})

tab <- bind_rows(out)

### save
saveRDS(tab, "gen_rep_rhog.rds")

# plot
# ggplot(tab, aes(rhog, rhog_hat, group = as.factor(rhog))) + geom_boxplot() + facet_wrap(~ n, nrow = 1) + geom_abline(linetype = 3) + coord_equal() + theme_minimal()

