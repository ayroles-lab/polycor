# https://stats.stackexchange.com/questions/13166/rs-lmer-cheat-sheet
# http://lme4.r-forge.r-project.org/slides/2011-03-16-Amsterdam/2Longitudinal.pdf
# https://rpsychologist.com/r-guide-longitudinal-lme-lmer

### inc 
library(tidyverse)

library(Matrix)
library(MASS)

library(lme4qtl)

path <- switch(Sys.info()[['nodename']],
  "~/git/hemostat/polycor/" )
src <- file.path(path, "relmer.R")
source(src)

### par
h1 <- 0.75
h2 <- 0.75
rg <- 0.75
re <- 0.1

sim_ids <- function(N) paste0("ID", 1:N)

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
    data_frame(id = ids, rid = ids, y1 = y1, y2 = y2, rep = reps)
  }) %>% 
  bind_rows %>%
  mutate(obs = paste0("obs", seq(n())),
    robs = obs)
}

# Start with a single trait
N <- 1000
kin <- sim_kin(h1, h2, rg, re, N)
dat <- sim_bivar(h1, h2, rg, re, N)

m <- relmer(y1 ~ (1|id), dat, relmat = list(id = kin))

# two traits & no replication
N <- 1000
kin <- sim_kin(h1, h2, rg, re, N)
dat <- sim_bivar(h1, h2, rg, re, N)
bdat <- gather(dat, tname, tvalue, y1, y2) 

m1 <- relmer(tvalue ~ tname - 1 + (0 + tname | id) + (0 + tname | rid), data = bdat, weights = rep(1e10, nrow(bdat)), relmat = list(id = kin), calc.derivs = FALSE)

# two traits + replication
N <- 250
R <- 2
kin <- sim_kin(h1, h2, rg, re, N)
dat <- sim_bivar(h1, h2, rg, re, N, R)
bdat <- gather(dat, tname, tvalue, y1, y2) 

kinobs <- kin[dat$id, dat$id] # kinship for observations with reptetive IDs
rownames(kinobs) <- colnames(kinobs) <- dat$obs

m2 <- relmer(tvalue ~ tname - 1 + (0 + tname | rid) + (0 + tname | robs), data = bdat, weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE)

m3 <- relmer(tvalue ~ tname - 1 + (0 + tname | obs) + (0 + tname | rid) + (0 + tname | robs), data = bdat, weights = rep(1e10, nrow(bdat)), relmat = list(obs = kinobs), calc.derivs = FALSE, verbose = 2)

m4 <- relmer(tvalue ~ tname - 1 + (0 + tname | id) + (0 + tname | rid) + (0 + tname | robs), data = bdat, weights = rep(1e10, nrow(bdat)), relmat = list(id = kin), calc.derivs = FALSE, verbose = 2)

m5 <- relmer(tvalue ~ tname - 1 + (0 + tname | id) + (0 + tname | obs) + (0 + tname | rid) + (0 + tname | robs), data = bdat, weights = rep(1e10, nrow(bdat)), relmat = list(id = kin, obs = kinobs), calc.derivs = FALSE, verbose = 2)


stop()

### testing
source("../../relmer.R")
m <- relmer(y1 ~ (1|id), dat, relmat = list(id = kin))

relmat <- list(id = kin)
names(relmat) <- "id:rep"
m <- relmer(y1 ~ (1|id:rep), dat, relmat = relmat)

stop()

# Example of data sim. & model fit (no rep)
N <- 100
kin <- sim_kin(h1, h2, rg, re, N)
dat <- sim_bivar(h1, h2, rg, re, N)
bdat <- gather(dat, tname, trait, y1, y2)

#m1 <- relmer(trait ~ (0 + tname|id) + (0 + tname|rid), bdat, weights = rep(1e10, nrow(bdat)), relmat = list(id = kin), calc.derivs = FALSE)
m1 <- relmer(trait ~ (0 + tname|id), bdat, relmat = list(id = kin))

# Example of data sim. & model fit (with rep)

N <- 100
R <- 10
kin <- sim_kin(h1, h2, rg, re, N)
dat <- sim_bivar(h1, h2, rg, re, N, R)
bdat <- gather(dat, tname, trait, y1, y2)

#m2 <- relmer(trait ~ (0 + tname|id:rep) + (0 + tname|rid:rep), bdat, weights = rep(1e10, nrow(bdat)), relmat = list(id = kin), calc.derivs = FALSE)
m2 <- relmer(trait ~ (0 + tname|id:rep), bdat, relmat = list(id = kin))
