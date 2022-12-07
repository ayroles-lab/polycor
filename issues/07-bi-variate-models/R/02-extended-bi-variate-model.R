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
h2.add <- 0.4
h2.cov <- 0.4

nrep <- 1
n <- 200

### data simulation (random-slope model)
set.seed(10)
simdat <- genCov.simulatePairKin(n_ind = n, n_perInd = nrep, 
  h2.cov = h2.cov, h2.add = h2.add)

dat <- simdat$pheno %>% as_tibble %>%
  mutate(
    trait1sc = as.numeric(scale(trait1)), trait2sc = as.numeric(scale(trait2_withAdd)),
    id = as.character(id), id1 = id, id2 = id, rid = id)
G <- simdat$G

bdat <- gather(dat, tname, tvalue, trait1sc, trait2sc)

### (0) random-slope model (data fitting model = data sim. model)
#m0 <- relmer(trait2 ~ (0 + trait1|id) + (1|rid), dat,
#  relmat = list(id = G), calc.derivs = FALSE)
#
#getME(m0, "Ztlist") %>% names
#getME(m0, "Ztlist") %>% .[[1]] %>% crossprod %>% image


### (1) the standard bi-variate
m1 <- relmer(tvalue ~ -1 + tname + (0 + tname | id) + (1 | rid), 
  data = bdat, relmat = list(id = G), calc.derivs = FALSE)

VarCorr(m1)

### (2) extended bi-variate #1
bdat <- mutate(bdat, id12 = paste(id, tname, sep = "_"))

sigma12 <- matrix(c(0, 1, 1, 0), 2, 2)
G12 <- kronecker(sigma12, G)
names12 <- lapply(unique(bdat$tname), function(t) paste(rownames(G), t, sep = "_")) %>% unlist
rownames(G12) <- colnames(G12) <- names12

m2 <- relmer(tvalue ~ -1 + tname + 
    (0 + dummy(tname, "trait1sc") | id1) + (0 + dummy(tname, "trait2sc") | id2) + 
    (1 | id12) + 
    (1 | rid),
  data = bdat, relmat = list(id1 = G, id2 = G, id12 = G12), 
  calc.derivs = FALSE,
  method.relfac = "svd")

# check variance components & comare to `mod1`
VarProp(m2)

### (3) extended bi-variate #2
bdat <- mutate(bdat, id12 = paste(id, tname, sep = "_"), id21 = id12)

sigma12 <- matrix(c(0, 0, 1, 0), 2, 2)
G12 <- kronecker(sigma12, G)
names12 <- lapply(unique(bdat$tname), function(t) paste(rownames(G), t, sep = "_")) %>% unlist
rownames(G12) <- colnames(G12) <- names12

sigma21 <- matrix(c(0, 1, 0, 0), 2, 2)
G21 <- kronecker(sigma21, G)
names21<- lapply(unique(bdat$tname), function(t) paste(rownames(G), t, sep = "_")) %>% unlist
rownames(G21) <- colnames(G21) <- names21

m3 <- relmer(tvalue ~ -1 + tname + 
    (0 + dummy(tname, "trait1sc") | id1) + (0 + dummy(tname, "trait2sc") | id2) +
    (1 | id12) + (1 | id21) + 
    (1 | rid),
  data = bdat, relmat = list(id1 = G, id2 = G, id12 = G12, id21 = G21), 
  calc.derivs = FALSE,
  method.relfac = "svd")

# check variance components & comare to `mod1`
VarProp(m3)

### (4) 
m4 <- relmer(tvalue ~ -1 + tname + 
    (0 + tname | id) +  
    (1 | id12) + (1 | id21) + 
    (1 | rid),
  data = bdat, relmat = list(id = G, id12 = G12, id21 = G21), 
  calc.derivs = FALSE,
  method.relfac = "svd")

# check variance components & comare to `mod1`
VarProp(m4)

