---
title: "Data simulations: pleiotropy + individual-level correlations"
author: "Simon Forsberg"
date: "2019-02-25"
output:
  html_document:
    theme: default
    toc: true
    keep_md: true
---



Simulations based on actual genotypes in outbred but stratified populations. Implemented in the function genBivarCov.simulatePair_outbred. The genotype based GRMs are dense, which means that the modeling fitting takes much longer. Because of this, I could only simulate relatively low sample sizes (see section Runtime). Ideally, each scenario below should be accompanied by a corresponding causal graph (TO DO).

#### Brief summary
* n simulated SNPs per individual
* GRM calculated from SNPs
* trait1 ~ g.add.t1 + e.t1
* trait2 ~ g.add.t2 + g.cov * trait1 + e.t2

**Where:**

* g.cov = polygenic effect on cor(trait1, trait2)
* g.add = additive polygenic effect on trait1/trait2.
* rhog = Fraction of overlapping causal SNPs

Polygenic effects per individual are the sum of the individual SNP effects

Note that rhog is defined as the fraction of overlapping causal loci. This should be proportional but not necessarily equal to the genetic correlation. When rhog > 0, h2.add.t2 and h2.cov should also be proportional but not equal to the respective true parameter values.

## General Observations
The slope_add model again seem robust in all the tested scenarios. At high rhog, h2.add.t2 estimates start to look a little weird. But I'm not sure what the true h2.add.t2 is in this scenario (see above). TBD

The bivar is a little bit more problematic:

* There is often a downward bias in the estimates of h2.add.t2
* When rhog = 0, the estimated rhog is all over the place, often 1 or -1. Could we do some significance test of this estimate to rule these out?
* When true rhog > 0, the rhog estimates are still pretty shaky (see scenario 3)

## Include


```r
library(tidyverse)
library(magrittr)

library(Matrix)
library(MASS)

library(nlme)
library(lme4)

library(cowplot)
library(wesanderson)
theme_set(theme_cowplot(font_size = 10))

# change to your path here, e.g. based on `nodemame`
path <- switch(Sys.info()[['nodename']],
  "~/git/hemostat/polycor/", 'gen-sflaptop' = '/media/sf9/ExtraDrive1/Dropbox/Research_170228/covQTL/bin/polycor/')
  
src <- file.path(path, "relmer.R")
source(src)

src <- file.path(path, "simFunctions.R")
source(src)

src <- file.path(path, "plotFunctions.R")
source(src)
```

## Parameters
These parameters are common to all the simulations


```r
testing <- F

n_loci <- 200
nrPops <- 10
Fst <- .3

#Sample sizes and number of simulations
nrep <- 10
vals_n <- c(50, 100)
vals_sim <- 1:5

if(testing) {
  nrep <- 2
  vals_n <- c(10, 20)
  vals_sim <- 1:2
}
```

## Functions

A function to fit the bivar and slope_add models to simulated data and wrap up the output. The following estimates are returned from the two models

| **bivar** | **slope_add** |
|--|--|
| rhog | h2.cov |
| h2.add.t1 | -- |
| h2.add.t2 | h2.add.t2 |


```r
run_models <- function(models, dat, G)
{
  lapply(models, function(model) {
    tab <- switch(model,
      "bivar" = {
        bdat <- gather(dat, tname, tvalue, trait1, trait2)
          mod <- relmer(tvalue ~ -1 + tname + (0 + tname | id) + (0 + tname | rid), 
            data = bdat, relmat = list(id = G), calc.derivs = FALSE)
          
          #Estimated variance components
          vf <- VarProp(mod) %>% filter(grp == "id")
          sigma.t1 <- vf$vcov[vf$var1 == 'tnametrait1' & is.na(vf$var2)]
          sigma.t2 <- vf$vcov[vf$var1 == 'tnametrait2' & is.na(vf$var2)]
          sigma.t1t2 <- vf$vcov[which(vf$var1 == 'tnametrait1' & vf$var2 == 'tnametrait2')]
          
          #Estimated h2s and rhog
          rhog <- sigma.t1t2 / sqrt(sigma.t1 * sigma.t2)
          h2.add.t1 <- vf$prop[vf$var1 == 'tnametrait1' & is.na(vf$var2)]
          h2.add.t2 <- vf$prop[vf$var1 == 'tnametrait2' & is.na(vf$var2)]
          
          data_frame(rhog = rhog, h2.add.t1 = h2.add.t1, h2.add.t2 = h2.add.t2)
      },
      "slope_add" = {
          mod <- relmer(trait2 ~ (1|rid) + (0 + trait1|id), data = dat, 
            relmat = list(id = G, rid = G))
        
          #Estimates of h2
          vf <- VarProp(mod)
          h2.add.t2 <- vf$prop[vf$grp == 'rid']
          h2.cov <- vf$prop[vf$grp == 'id']
          
          data_frame(h2.cov = h2.cov, h2.add.t2 = h2.add.t2)
      },
      stop("error in switch: unknown `model`"))
        
    mutate(tab, model = model)
  }) %>% bind_rows
}
```

# Simulation 1: h2.cov + h2.add.t2


```r
#Parameters
h2_cov <- c(0.3, 0.6)
h2_add_t1 <- 0
h2_add_t2 <- c(1e-4, 0.1, 0.2, 0.3)
rhog <- 0
resFile <- '../results/sim1.RData'

if(!file.exists(resFile)){
  grid <- expand.grid(n = as.integer(vals_n), sim = vals_sim, h2_add_t1 = h2_add_t1, h2_add_t2 = h2_add_t2, h2_cov = h2_cov, rhog = rhog)
  SimRes <- data.frame(grid, runtime = NA,
                       h2.add.t2_slopeMod = NA, 
                       h2.add.t1_bivar = NA, h2.add.t2_bivar = NA,
                       h2.cov_slopeMod = NA, rhog_bivar = NA)
  
  for(i in 1:nrow(grid)){
    n <- grid$n[i]
    sim <- grid$sim[i]
    h2.add.t2 <- grid$h2_add_t2[i]
    h2.cov <- grid$h2_cov[i]
    message(i, "/", nrow(grid), "\n")
    
    #Run simulation
    simdat <- genBivarCov.simulatePair_outbred(n_ind = n, n_perInd = nrep, 
                                               h2.add.t1 = h2_add_t1, h2.add.t2 = h2.add.t2,
                                               h2.cov = h2.cov, rhog = rhog, 
                                               nrPops = nrPops, Fst = Fst, postproc = T)
    
    #Fit models
    models <- c("bivar", "slope_add")
    start <- Sys.time()
    simdat.modFit <- run_models(models, simdat$pheno, simdat$G)
    end <- Sys.time()
    time <- end - start
    
    #Store results
    SimRes$h2.add.t2_slopeMod[i] <- simdat.modFit$h2.add.t2[simdat.modFit$model == 'slope_add']
    SimRes$h2.add.t1_bivar[i] <- simdat.modFit$h2.add.t1[simdat.modFit$model == 'bivar']
    SimRes$h2.add.t2_bivar[i] <- simdat.modFit$h2.add.t2[simdat.modFit$model == 'bivar']
    SimRes$h2.cov_slopeMod[i] <- simdat.modFit$h2.cov[simdat.modFit$model == 'slope_add']
    SimRes$rhog_bivar[i] <- simdat.modFit$rhog[simdat.modFit$model == 'bivar']
    SimRes$runtime[i] <- as.numeric(time, unit = 'mins')
  }
  save(list = 'SimRes', file = resFile)
}else
  load(resFile)
```


## h2.add.t2
* x-axis = h2.add.t2 true
* y-axis = h2.add.t2 estimated

![](../figures/plot_sim1_addT2-1.png)<!-- -->

## h2.cov & rhog
* x-axis = h2.add.t2 true
* y-axis = h2.cov (slope_add model), rhog (bivar model)

Dashed lines indicate simulated values for h2.cov and rhog
![](../figures/plot_sim1_h2cov_rhog-1.png)<!-- -->

# Simulation 2: h2.cov + h2.add.t1 + h2.add.t2


```r
#Parameters
h2_cov <- c(0.3, 0.6)
h2_add_t1 <- .3
h2_add_t2 <- c(1e-4, 0.1, 0.2, 0.3)
rhog <- 0
resFile <- '../results/sim2.RData'

if(!file.exists(resFile)){
  grid <- expand.grid(n = as.integer(vals_n), sim = vals_sim, h2_add_t1 = h2_add_t1, h2_add_t2 = h2_add_t2, h2_cov = h2_cov, rhog = rhog)
  SimRes <- data.frame(grid, runtime = NA,
                       h2.add.t2_slopeMod = NA, 
                       h2.add.t1_bivar = NA, h2.add.t2_bivar = NA,
                       h2.cov_slopeMod = NA, rhog_bivar = NA)
  
  for(i in 1:nrow(grid)){
    n <- grid$n[i]
    sim <- grid$sim[i]
    h2.add.t2 <- grid$h2_add_t2[i]
    h2.cov <- grid$h2_cov[i]
    message(i, "/", nrow(grid), "\n")
    
    #Run simulation
    simdat <- genBivarCov.simulatePair_outbred(n_ind = n, n_perInd = nrep, 
                                               h2.add.t1 = h2_add_t1, h2.add.t2 = h2.add.t2,
                                               h2.cov = h2.cov, rhog = rhog, 
                                               nrPops = nrPops, Fst = Fst, postproc = T)
    
    #Fit models
    models <- c("bivar", "slope_add")
    start <- Sys.time()
    simdat.modFit <- run_models(models, simdat$pheno, simdat$G)
    end <- Sys.time()
    time <- end - start
    
    #Store results
    SimRes$h2.add.t2_slopeMod[i] <- simdat.modFit$h2.add.t2[simdat.modFit$model == 'slope_add']
    SimRes$h2.add.t1_bivar[i] <- simdat.modFit$h2.add.t1[simdat.modFit$model == 'bivar']
    SimRes$h2.add.t2_bivar[i] <- simdat.modFit$h2.add.t2[simdat.modFit$model == 'bivar']
    SimRes$h2.cov_slopeMod[i] <- simdat.modFit$h2.cov[simdat.modFit$model == 'slope_add']
    SimRes$rhog_bivar[i] <- simdat.modFit$rhog[simdat.modFit$model == 'bivar']
    SimRes$runtime[i] <- as.numeric(time, unit = 'mins')
  }
  save(list = 'SimRes', file = resFile)
}else
  load(resFile)
```

## h2.add.t2
* x-axis = h2.add.t2 true
* y-axis = h2.add.t2 estimated

![](../figures/plot_sim2_addT2-1.png)<!-- -->

## h2.cov & rhog
* x-axis = h2.add.t2 true
* y-axis = h2.cov (slope_add model), rhog (bivar model)

Dashed lines indicate simulated values for h2.cov and rhog

![](../figures/plot_sim2_h2cov_rhog-1.png)<!-- -->

# Simulation 3: h2.cov + h2.add.t1 + h2.add.t2 + rhog


```r
#Parameters
h2_cov <- c(0.3, 0.6)
h2_add_t1 <- .3
h2_add_t2 <- c(1e-4, 0.1, 0.2, 0.3)
rhog <- .5
resFile <- '../results/sim3.RData'

if(!file.exists(resFile)){
  grid <- expand.grid(n = as.integer(vals_n), sim = vals_sim, h2_add_t1 = h2_add_t1, h2_add_t2 = h2_add_t2, h2_cov = h2_cov, rhog = rhog)
  SimRes <- data.frame(grid, runtime = NA,
                       h2.add.t2_slopeMod = NA, 
                       h2.add.t1_bivar = NA, h2.add.t2_bivar = NA,
                       h2.cov_slopeMod = NA, rhog_bivar = NA)
  
  for(i in 1:nrow(grid)){
    n <- grid$n[i]
    sim <- grid$sim[i]
    h2.add.t2 <- grid$h2_add_t2[i]
    h2.cov <- grid$h2_cov[i]
    message(i, "/", nrow(grid), "\n")
    
    #Run simulation
    simdat <- genBivarCov.simulatePair_outbred(n_ind = n, n_perInd = nrep, 
                                               h2.add.t1 = h2_add_t1, h2.add.t2 = h2.add.t2,
                                               h2.cov = h2.cov, rhog = rhog, 
                                               nrPops = nrPops, Fst = Fst, postproc = T)
    
    #Fit models
    models <- c("bivar", "slope_add")
    start <- Sys.time()
    simdat.modFit <- run_models(models, simdat$pheno, simdat$G)
    end <- Sys.time()
    time <- end - start
    
    #Store results
    SimRes$h2.add.t2_slopeMod[i] <- simdat.modFit$h2.add.t2[simdat.modFit$model == 'slope_add']
    SimRes$h2.add.t1_bivar[i] <- simdat.modFit$h2.add.t1[simdat.modFit$model == 'bivar']
    SimRes$h2.add.t2_bivar[i] <- simdat.modFit$h2.add.t2[simdat.modFit$model == 'bivar']
    SimRes$h2.cov_slopeMod[i] <- simdat.modFit$h2.cov[simdat.modFit$model == 'slope_add']
    SimRes$rhog_bivar[i] <- simdat.modFit$rhog[simdat.modFit$model == 'bivar']
    SimRes$runtime[i] <- as.numeric(time, unit = 'mins')
  }
  save(list = 'SimRes', file = resFile)
}else
  load(resFile)
```

## h2.add.t2
* x-axis = h2.add.t2 true
* y-axis = h2.add.t2 estimated

![](../figures/plot_sim3_addT2-1.png)<!-- -->

## h2.cov & rhog
* x-axis = h2.add.t2 true
* y-axis = h2.cov (slope_add model), rhog (bivar model)

Dashed lines indicate simulated values for h2.cov and rhog

![](../figures/plot_sim3_h2cov_rhog-1.png)<!-- -->

# Simulation 4: h2.cov + h2.add.t1 + h2.add.t2 + rhog

Same as simulation 3 but with rhog = 0.9


```r
#Parameters
h2_cov <- c(0.3, 0.6)
h2_add_t1 <- .3
h2_add_t2 <- c(1e-4, 0.1, 0.2, 0.3)
rhog <- .9
resFile <- '../results/sim4.RData'

if(!file.exists(resFile)){
  grid <- expand.grid(n = as.integer(vals_n), sim = vals_sim, h2_add_t1 = h2_add_t1, h2_add_t2 = h2_add_t2, h2_cov = h2_cov, rhog = rhog)
  SimRes <- data.frame(grid, runtime = NA,
                       h2.add.t2_slopeMod = NA, 
                       h2.add.t1_bivar = NA, h2.add.t2_bivar = NA,
                       h2.cov_slopeMod = NA, rhog_bivar = NA)
  
  for(i in 1:nrow(grid)){
    n <- grid$n[i]
    sim <- grid$sim[i]
    h2.add.t2 <- grid$h2_add_t2[i]
    h2.cov <- grid$h2_cov[i]
    message(i, "/", nrow(grid), "\n")
    
    #Run simulation
    simdat <- genBivarCov.simulatePair_outbred(n_ind = n, n_perInd = nrep, 
                                               h2.add.t1 = h2_add_t1, h2.add.t2 = h2.add.t2,
                                               h2.cov = h2.cov, rhog = rhog, 
                                               nrPops = nrPops, Fst = Fst, postproc = T)
    
    #Fit models
    models <- c("bivar", "slope_add")
    start <- Sys.time()
    simdat.modFit <- run_models(models, simdat$pheno, simdat$G)
    end <- Sys.time()
    time <- end - start
    
    #Store results
    SimRes$h2.add.t2_slopeMod[i] <- simdat.modFit$h2.add.t2[simdat.modFit$model == 'slope_add']
    SimRes$h2.add.t1_bivar[i] <- simdat.modFit$h2.add.t1[simdat.modFit$model == 'bivar']
    SimRes$h2.add.t2_bivar[i] <- simdat.modFit$h2.add.t2[simdat.modFit$model == 'bivar']
    SimRes$h2.cov_slopeMod[i] <- simdat.modFit$h2.cov[simdat.modFit$model == 'slope_add']
    SimRes$rhog_bivar[i] <- simdat.modFit$rhog[simdat.modFit$model == 'bivar']
    SimRes$runtime[i] <- as.numeric(time, unit = 'mins')
  }
  save(list = 'SimRes', file = resFile)
}else
  load(resFile)
```

## h2.add.t2
* x-axis = h2.add.t2 true
* y-axis = h2.add.t2 estimated

![](../figures/plot_sim4_addT2-1.png)<!-- -->

## h2.cov & rhog
* x-axis = h2.add.t2 true
* y-axis = h2.cov (slope_add model), rhog (bivar model)

Dashed lines indicate simulated values for h2.cov and rhog

![](../figures/plot_sim4_h2cov_rhog-1.png)<!-- -->

# Runtime
Time to fit both models (run_models function above). To finish this in a reasonable time, I kept the sample sizes pretty low. Judging by just these two sample sizes, runtime does not scale linearly with n. We might thus run into problems when trying to fit the models to large datasets. 

![](../figures/plot_runtime-1.png)<!-- -->

Average runtimes (minutes):

| **n = 50** | **n = 100** |
|--|--|
| 0.109 | 0.76 |
