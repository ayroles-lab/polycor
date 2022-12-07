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
h2.cov <- 0.8

nrep <- 5
n <- 50

### data simulation (pleiotropy)
set.seed(20)
simdat <- pleiotropy.simulatePairKin(n_ind = n, n_perInd = nrep, 
  h2 = h2, overlap = rhog)

dat <- simdat$pheno %>% as_tibble %>%
  mutate(
    trait1sc = as.numeric(scale(trait1)), trait2sc = as.numeric(scale(trait2)),
    id = as.factor(as.character(id)), obs = as.factor(as.character(obs)), 
    id1 = id, id2 = id, rid = id, robs = obs,
    obscov = paste(obs, tname, sep = "_"))

G <- simdat$G
Gobs <- G[dat$id, dat$id] # kinship for observations with reptetive IDs
rownames(Gobs) <- colnames(Gobs) <- dat$obs

bdat <- gather(dat, tname, tvalue, trait1sc, trait2sc)

### (1) the standard bi-variate
m1 <- relmer(tvalue ~ -1 + tname + (0 + tname | id) + (1 | rid), 
  data = bdat, relmat = list(id = G), calc.derivs = FALSE)

VarCorr(m1)

### (2) like model 1, but obs are used instead of id
m2 <- relmer(tvalue ~ -1 + tname + (0 + tname | obs) + (1 | robs), 
  data = bdat, relmat = list(obs = Gobs), calc.derivs = FALSE)

VarCorr(m2)

### (3) neste bi-variate: id + obs
m3 <- relmer(tvalue ~ -1 + tname + (0 + tname | id) + (1 | rid) + (0 + tname | obs) + (1 | robs), 
  data = bdat, relmat = list(id = G, obs = Gobs), calc.derivs = FALSE)

VarCorr(m3)
