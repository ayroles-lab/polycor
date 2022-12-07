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

src <- file.path(path, "simFunctions.R")
source(src)

### par
testing <- FALSE

h2 <- 0.75
nrep <- 4
vals_n <- c(100, 500, 1e3)
vals_rhog <- c(0.1, 0.5, 0.9)
vals_sim <- 1:5

if(testing) {
  nrep <- 2
  vals_n <- c(200)
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
  
  simdat <- pleiotropy.simulatePairKin(n_ind = n, n_perInd = nrep, 
    h2 = h2, overlap = rhog)

  dat <- simdat$pheno %>% 
    mutate(id = as.character(id), rid = id)
  G <- simdat$G

  bdat <- gather(dat, tname, tvalue, trait1, trait2)
  mod <- relmer(tvalue ~ -1 + obs + tname + (0 + tname | id) + (0 + tname | rid), 
    data = bdat, relmat = list(id = G), calc.derivs = FALSE, verbose = 2)
    
  vf <- VarCorr(mod) %>% as.data.frame %>% 
    filter(grp == "id" & var1 == "tnametrait1" & var2 == "tnametrait2")
  rhog_hat <- vf$sdcor
  
  data_frame(sim = sim, n = n/nrep, nrep = nrep, rhog = rhog, rhog_hat = rhog_hat)
})

tab <- bind_rows(out)

