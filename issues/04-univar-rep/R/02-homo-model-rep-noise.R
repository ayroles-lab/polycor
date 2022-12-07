# in simulation of rep. observations:
# - the genetic part is simulated once (common across rep)
# - the residual part is simulated for each rep.

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
testing <- TRUE

vals_n <- 500
vals_rep <- 4
vals_h <- c(0.1, 0.5, 0.9)
vals_sim <- 1:5

if(testing) {
  vals_n <- 40
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

sim_univar <- function(h, N, rep = 1, L = 5, scale = TRUE) 
{
  stopifnot(!(N %% L)) # family size = L
  Nl <- N / L
  
  K <- sim_kin(h, N, L)
  
  # V = [K; I]
  V_gen <- h*K 
  V_res <- (1-h)*Diagonal(N)
  
  y_gen <- mvrnorm(1, rep(0, N), V_gen)
  lapply(seq(1, rep), function(r) {
    y_res <- mvrnorm(1, rep(0, N), V_res)
    
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

### no rep
N <- 1000; h <- 0.8; w <- rep(1e10, N)
kin <- sim_kin(h, N)
dat <- sim_univar(h, N)

(mod <- relmer(y ~ (1|id) + (1|rid), dat, relmat = list(id = kin), weights = w, calc.derivs = FALSE))

### rep
N <- 1000; R <- 2; h <- 0.8; w <- rep(1e10, N*R)
kin <- sim_kin(h, N)
dat <- sim_univar(h, N, R)

(mod1 <- relmer(y ~ (1|id) + (1|rid), dat, relmat = list(id = kin), weights = w, calc.derivs = FALSE))

(mod2 <- relmer(y ~ (1|id) + (1|rid), dat, relmat = list(id = kin)))

### rep
N <- 250; R <- 4; h <- 0.8; w <- rep(1e10, N*R)
kin <- sim_kin(h, N)
dat <- sim_univar(h, N, R)

(mod3 <- relmer(y ~ (1|id) + (1|rid), dat, relmat = list(id = kin), weights = w, calc.derivs = FALSE))

(mod4 <- relmer(y ~ (1|id) + (1|rid), dat, relmat = list(id = kin)))
