### inc 
library(tidyverse)
library(magrittr)

library(Matrix)
library(MASS)

library(lme4)

# change to your path here, e.g. based on `nodemame`
path <- switch(Sys.info()[['nodename']],
  "~/git/hemostat/polycor/" )
  
src <- file.path(path, "relmer.R")
source(src)

src <- file.path(path, "simFunctions.R")
source(src)

### par
h2 <- 0.75
nrep <- 1
n <- 20
rhog <- 0.75

### simulate data
set.seed(10)
simdat <- pleiotropy.simulatePairKin(n_ind = n, n_perInd = nrep, 
    h2 = h2, overlap = rhog)

dat <- simdat$pheno %>% as_tibble %>%
  mutate(
    trait1 = as.numeric(scale(trait1)), trait2 = as.numeric(scale(trait2)),
    id = as.character(id), id1 = id, id2 = id, rid = id)
G <- simdat$G

bdat <- gather(dat, tname, tvalue, trait1, trait2)

### (2) fit Bivariate Model 1 (see issue #5)
mod1 <- relmer(tvalue ~ -1 + tname + (0 + tname | id) + (0 + tname | rid), 
  data = bdat, relmat = list(id = G), weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE)

# alternative way to fit `mod1`:
#mod1 <- relmer(tvalue ~ -1 + tname + (0 + tname | id) + (1 | rid), data = bdat, relmat = list(id = G), calc.derivs = FALSE)


# 
# G = ZZ' matrix for trait1 (not scaled by h2)
getME(mod1, "Ztlist") %>% names
getME(mod1, "Ztlist") %>% .[[1]] %>% crossprod %>% image

ranef(mod1) %>% str

### (2) fit Bivariate Model 1 (dropping the term with genetic covariance)
mod2 <- relmer(tvalue ~ -1 + tname + 
    (0 + dummy(tname, "trait1") | id1) + 
    (0 + dummy(tname, "trait2") | id2) +
    (0 + tname | rid), 
  data = bdat, relmat = list(id1 = G, id2 = G), 
  weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE)

# check there is no genetic cov. term
VarCorr(mod2)

### (3) fit Bivariate Model 1 (separating three genetic terms)
# - here we observe the variance leakage: gcov parts vs. marginal h1 & h2
bdat <- mutate(bdat, id12 = paste(id, tname, sep = "_"))

sigma12 <- matrix(c(0, 1, 1, 0), 2, 2)
G12 <- kronecker(sigma12, G)
names12 <- lapply(unique(bdat$tname), function(t) paste(rownames(G), t, sep = "_")) %>% unlist
rownames(G12) <- colnames(G12) <- names12

mod3 <- relmer(tvalue ~ -1 + tname + 
    (0 + dummy(tname, "trait1") | id1) + 
    (0 + dummy(tname, "trait2") | id2) +
    (1 | id12) + 
    (0 + tname | rid),
  data = bdat, relmat = list(id1 = G, id2 = G, id12 = G12), 
  weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE,
  method.relfac = "svd")

# check variance components & comare to `mod1`
VarCorr(mod3)

### (4) fit Bivariate Model 1+3 (pleiotropy + gcov)
# - here we observe that the additional covariance term (1|id12) adds little to the (1|id)
#   (that how we simulated data under pleiotropy model)
mod4 <- relmer(tvalue ~ -1 + tname + 
    (0 + tname | id) + (1 | id12) + 
    (0 + tname | rid),
  data = bdat, relmat = list(id = G, id12 = G12), 
  weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE,
  method.relfac = "svd")

# check variance components & comare to `mod1`
VarCorr(mod4)

### (5.1) 
set.seed(10)
simdat <- genCov.simulatePairKin(n_ind = n, n_perInd = nrep, 
  h2.cov = 0.8, h2.add = 0.1)

dat <- simdat$pheno %>% as_tibble %>%
  mutate(
    trait1 = as.numeric(scale(trait1)), trait2 = as.numeric(scale(trait2)),
    id = as.character(id), id1 = id, id2 = id, rid = id)
G <- simdat$G

bdat <- gather(dat, tname, tvalue, trait1, trait2)

m1 <- relmer(tvalue ~ -1 + tname + (0 + tname | id) + (0 + tname | rid), 
  data = bdat, relmat = list(id = G), weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE)

VarCorr(m1)

### (5.2)
bdat <- mutate(bdat, id12 = paste(id, tname, sep = "_"))

sigma12 <- matrix(c(0, 1, 1, 0), 2, 2)
G12 <- kronecker(sigma12, G)
names12 <- lapply(unique(bdat$tname), function(t) paste(rownames(G), t, sep = "_")) %>% unlist
rownames(G12) <- colnames(G12) <- names12

m2 <- relmer(tvalue ~ -1 + tname + 
    (0 + tname | id) + (1 | id12) + 
    (0 + tname | rid),
  data = bdat, relmat = list(id = G, id12 = G12), 
  weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE,
  method.relfac = "svd")

# check variance components & comare to `mod1`
VarCorr(m2)

