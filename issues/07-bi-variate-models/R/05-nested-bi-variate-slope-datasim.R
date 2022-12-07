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
h2 <- 0.8
rhog <- 0.8

h2.cov <- 0.4
h2.add <- 0.4

nrep <- 5
n <- 50

### data simulation (random-slope model)
set.seed(20)
simdat <- genCov.simulatePairKin(n_ind = n, n_perInd = nrep, 
  h2.cov = h2.cov, h2.add = h2.add)

dat <- simdat$pheno %>% as_tibble %>%
  mutate(
    trait1sc = as.numeric(scale(trait1)), trait2sc = as.numeric(scale(trait2_withAdd)),
    id = as.factor(as.character(id)), obs = as.factor(as.character(obs)), 
    id1 = id, id2 = id, 
    obs1 = obs, obs2 = obs,
    rid = id, robs = obs)

G <- simdat$G

Gobs <- G[dat$id, dat$id] # kinship for observations with reptetive IDs
rownames(Gobs) <- colnames(Gobs) <- dat$obs

bdat <- gather(dat, tname, tvalue, trait1sc, trait2sc) %>%
  mutate(
    idcov = paste(id, tname, sep = "_"),
    obscov = paste(obs, tname, sep = "_"))

sigmacov <- matrix(c(0, 1, 1, 0), 2, 2)
Gcov <- kronecker(sigmacov, G)
namescov <- lapply(unique(bdat$tname), function(t) paste(rownames(G), t, sep = "_")) %>% unlist
rownames(Gcov) <- colnames(Gcov) <- namescov

Gcovobs <- Gcov[bdat$idcov, bdat$idcov] # kinship for observations with reptetive IDs
rownames(Gcovobs) <- colnames(Gcovobs) <- bdat$obscov

### (1) the standard bi-variate
m1 <- relmer(tvalue ~ -1 + tname + (0 + tname | id) + (1 | rid), 
  data = bdat, relmat = list(id = G), calc.derivs = FALSE)

VarCorr(m1)

### (2) like model 1, but obs are used instead of id
m2 <- relmer(tvalue ~ -1 + tname + (0 + tname | obs) + (1 | robs), 
  data = bdat, relmat = list(obs = Gobs), calc.derivs = FALSE)

VarCorr(m2)

### (3) nested bi-variate: id + obs
m3 <- relmer(tvalue ~ -1 + tname + (0 + tname | id) + (1 | rid) + (0 + tname | obs) + (1 | robs), 
  data = bdat, relmat = list(id = G, obs = Gobs), calc.derivs = FALSE)

VarCorr(m3)

### (4) nested bi-variate + extended (gen. cov is a separate term)
m4 <- relmer(tvalue ~ -1 + tname + 
    (0 + dummy(tname, "trait1sc") | id1) + (0 + dummy(tname, "trait2sc") | id2) + 
    (1 | idcov) + 
    (1 | rid) + 
    (0 + dummy(tname, "trait1sc") | obs1) + (0 + dummy(tname, "trait2sc") | obs2) + 
    (1 | obscov) + 
    (1 | robs), 
  data = bdat, relmat = list(id1 = G, id2 = G, idcov = Gcov,
    obs1 = Gobs, obs2 = Gobs, obscov = Gcovobs), 
  calc.derivs = FALSE, method.relfac = "svd")

VarProp(m4) %>% as_data_frame
