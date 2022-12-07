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
h2.cov <- 0.8

nrep <- 5
n <- 50

### data simulation (random-slope model)
set.seed(20)
simdat <- genCov.simulatePairKin(n_ind = n, n_perInd = nrep, 
  h2.cov = h2.cov, h2.add = 0)

dat <- simdat$pheno %>% as_tibble %>%
  mutate(
    trait1sc = as.numeric(scale(trait1)), trait2sc = as.numeric(scale(trait2)),
    id = as.factor(as.character(id)), obs = as.factor(as.character(obs)), 
    id1 = id, id2 = id, rid = id)
G <- simdat$G

bdat <- gather(dat, tname, tvalue, trait1sc, trait2sc)

### (1) random slope
m1 <- relmer(trait2 ~ (0 + trait1|id), data = dat, relmat = list(id = G))

VarProp(m1)

### (2) the standard bi-variate
m2 <- relmer(tvalue ~ -1 + tname + (0 + tname | id) + (1 | rid), 
  data = bdat, relmat = list(id = G), calc.derivs = FALSE)

VarProp(m2)

### (3) extended bi-variate #1
bdat <- mutate(bdat, idcov = paste(id, tname, sep = "_"))

sigmacov <- matrix(c(0, 1, 1, 0), 2, 2)
Gcov <- kronecker(sigmacov, G)
namescov <- lapply(unique(bdat$tname), function(t) paste(rownames(G), t, sep = "_")) %>% unlist
rownames(Gcov) <- colnames(Gcov) <- namescov

m3 <- relmer(tvalue ~ -1 + tname + 
    (0 + dummy(tname, "trait1sc") | id1) + (0 + dummy(tname, "trait2sc") | id2) + 
    (1 | idcov) + 
    (1 | rid),
  data = bdat, relmat = list(id1 = G, id2 = G, idcov = Gcov), 
  calc.derivs = FALSE,
  method.relfac = "svd")

# check variance components & comare to `mod1`
VarProp(m3)


### (4) extended bi-variate #2
bdat <- mutate(bdat, id12 = paste(id, tname, sep = "_"), id21 = id12)

sigma12 <- matrix(c(0, 0, 1, 0), 2, 2)
G12 <- kronecker(sigma12, G)
names12 <- lapply(unique(bdat$tname), function(t) paste(rownames(G), t, sep = "_")) %>% unlist
rownames(G12) <- colnames(G12) <- names12

sigma21 <- matrix(c(0, 1, 0, 0), 2, 2)
G21 <- kronecker(sigma21, G)
names21<- lapply(unique(bdat$tname), function(t) paste(rownames(G), t, sep = "_")) %>% unlist
rownames(G21) <- colnames(G21) <- names21

m4 <- relmer(tvalue ~ -1 + tname + 
    (0 + dummy(tname, "trait1sc") | id1) + (0 + dummy(tname, "trait2sc") | id2) +
    (1 | id12) + (1 | id21) + 
    (1 | rid),
  data = bdat, relmat = list(id1 = G, id2 = G, id12 = G12, id21 = G21), 
  calc.derivs = FALSE,
  method.relfac = "svd")

# check variance components 
VarProp(m4)

### (5) extended bi-variate #3
bdat <- mutate(bdat, idcov = paste(id, tname, sep = "_"))

sigmacov <- matrix(c(0, 1, 1, 0), 2, 2)
Gcov <- kronecker(sigmacov, G)
namescov <- lapply(unique(bdat$tname), function(t) paste(rownames(G), t, sep = "_")) %>% unlist
rownames(Gcov) <- colnames(Gcov) <- namescov

m5 <- relmer(tvalue ~ -1 + tname + 
    (0 + tname|id) + 
    (1 | idcov) + 
    (1 | rid),
  data = bdat, relmat = list(id = G, idcov = Gcov), 
  calc.derivs = FALSE,
  method.relfac = "svd")

# check variance components 
VarProp(m5)

# > anova(m5, update(m5, ~ . -(1 | idcov)))

### (6) extended bi-variate with only cov. part
bdat <- mutate(bdat, idcov = paste(id, tname, sep = "_"))

sigmacov <- matrix(c(0, 1, 1, 0), 2, 2)
Gcov <- kronecker(sigmacov, G)
namescov <- lapply(unique(bdat$tname), function(t) paste(rownames(G), t, sep = "_")) %>% unlist
rownames(Gcov) <- colnames(Gcov) <- namescov

m6 <- relmer(tvalue ~ -1 + tname + 
    (1 | idcov) + 
    (1 | rid),
  data = bdat, relmat = list(idcov = Gcov), 
  calc.derivs = FALSE,
  method.relfac = "svd")

# check variance components 
VarProp(m6)

# > anova(m6, update(m6, ~ . -(1 | idcov)))
