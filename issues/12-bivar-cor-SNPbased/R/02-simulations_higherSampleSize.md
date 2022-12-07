---
title: "Data simulations: pleiotropy + individual-level correlations. Genotype based"
author: "Simon Forsberg"
date: "2019-03-12"
output:
  html_document:
    theme: default
    toc: true
    keep_md: true
---



Simulations based on actual genotypes in outbred but stratified populations. Implemented in the function genBivarCov.simulatePair_outbred. Here, I simulated larger sample sizes (n = 500 x 10 & 1000 x 10) and the models can take really long to fit. See the Runtime section below.

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
Largely the same conclusions as at lower sample sizes. The bivar model takes especially long to fit. In the worst case, one model fit took ~50h for n = 1000 x 10. See runtime section below.

### slope_add
Quite robust in all scenarios. If you squint, there might be a slight correlation between h2.add.t2 and estimated h2.cov when rhog > 0 (simulations 3 & 4). Not sure whats going on there but it dosn't look like it's a big problem. 

### bivar
* Often a downward bias in the estimates of h2.add.t2
* When rhog = 0 AND the heritability of one of the traits = 0, the estimated rhog is all over the place, often 1 or -1. 
    * I think the reason is that you divide by something close to zero when h2.add.t1 = 0 or h2.add.t2 = 0.
    * When both h2.add.t1 > 0 & h2.add.t2 > 0, the estimates of rhog are around zero (simulation 2).

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
vals_n <- c(500, 1000)
vals_sim <- 1:3

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
                    start <- Sys.time()
                    mod <- relmer(tvalue ~ -1 + tname + (0 + tname | id) + (0 + tname | rid), 
                                  data = bdat, relmat = list(id = G), calc.derivs = FALSE)
                    end <- Sys.time()
                    time <- end - start
                    
                    #Estimated variance components
                    vf <- VarProp(mod) %>% filter(grp == "id")
                    sigma.t1 <- vf$vcov[vf$var1 == 'tnametrait1' & is.na(vf$var2)]
                    sigma.t2 <- vf$vcov[vf$var1 == 'tnametrait2' & is.na(vf$var2)]
                    sigma.t1t2 <- vf$vcov[which(vf$var1 == 'tnametrait1' & vf$var2 == 'tnametrait2')]
                    
                    #Estimated h2s and rhog
                    rhog <- sigma.t1t2 / sqrt(sigma.t1 * sigma.t2)
                    h2.add.t1 <- vf$prop[vf$var1 == 'tnametrait1' & is.na(vf$var2)]
                    h2.add.t2 <- vf$prop[vf$var1 == 'tnametrait2' & is.na(vf$var2)]
                    
                    data_frame(rhog = rhog, h2.add.t1 = h2.add.t1, h2.add.t2 = h2.add.t2, runtime = as.numeric(time, unit = 'mins'))
                  },
                  "slope_add" = {
                    start <- Sys.time()
                    mod <- relmer(trait2 ~ (1|rid) + (0 + trait1|id), data = dat, 
                                  relmat = list(id = G, rid = G))
                    end <- Sys.time()
                    time <- end - start
                    
                    #Estimates of h2
                    vf <- VarProp(mod)
                    h2.add.t2 <- vf$prop[vf$grp == 'rid']
                    h2.cov <- vf$prop[vf$grp == 'id']
                    
                    data_frame(h2.cov = h2.cov, h2.add.t2 = h2.add.t2, runtime = as.numeric(time, unit = 'mins'))
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
resFile <- '../results/02-highN/sim1.RData'

if(!file.exists(resFile)){
  grid <- expand.grid(n = as.integer(vals_n), sim = vals_sim, h2_add_t1 = h2_add_t1, h2_add_t2 = h2_add_t2, h2_cov = h2_cov, rhog = rhog)
  SimRes <- data.frame(grid, runtime_slope = NA, runtime_bivar = NA,
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
    simdat.modFit <- run_models(models, simdat$pheno, simdat$G)

    #Store results
    SimRes$h2.add.t2_slopeMod[i] <- simdat.modFit$h2.add.t2[simdat.modFit$model == 'slope_add']
    SimRes$h2.add.t1_bivar[i] <- simdat.modFit$h2.add.t1[simdat.modFit$model == 'bivar']
    SimRes$h2.add.t2_bivar[i] <- simdat.modFit$h2.add.t2[simdat.modFit$model == 'bivar']
    SimRes$h2.cov_slopeMod[i] <- simdat.modFit$h2.cov[simdat.modFit$model == 'slope_add']
    SimRes$rhog_bivar[i] <- simdat.modFit$rhog[simdat.modFit$model == 'bivar']
    SimRes$runtime_slope[i] <- simdat.modFit$runtime[simdat.modFit$model == 'slope_add']
    SimRes$runtime_bivar[i] <- simdat.modFit$runtime[simdat.modFit$model == 'bivar']
  }
  save(list = 'SimRes', file = resFile)
}else
  load(resFile)
```


## h2.add.t2
* x-axis = h2.add.t2 true
* y-axis = h2.add.t2 estimated

![](../figures/02-highN/plot_sim1_addT2-1.png)<!-- -->

## h2.cov & rhog
* x-axis = h2.add.t2 true
* y-axis = h2.cov (slope_add model), rhog (bivar model)

Dashed lines indicate simulated values for h2.cov and rhog
![](../figures/02-highN/plot_sim1_h2cov_rhog-1.png)<!-- -->

# Simulation 2: h2.cov + h2.add.t1 + h2.add.t2


```r
#Parameters
h2_cov <- c(0.3, 0.6)
h2_add_t1 <- .3
h2_add_t2 <- c(1e-4, 0.1, 0.2, 0.3)
rhog <- 0
resFile <- '../results/02-highN/sim2.RData'

if(!file.exists(resFile)){
  grid <- expand.grid(n = as.integer(vals_n), sim = vals_sim, h2_add_t1 = h2_add_t1, h2_add_t2 = h2_add_t2, h2_cov = h2_cov, rhog = rhog)
  SimRes <- data.frame(grid, runtime_slope = NA, runtime_bivar = NA,
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
    simdat.modFit <- run_models(models, simdat$pheno, simdat$G)

    #Store results
    SimRes$h2.add.t2_slopeMod[i] <- simdat.modFit$h2.add.t2[simdat.modFit$model == 'slope_add']
    SimRes$h2.add.t1_bivar[i] <- simdat.modFit$h2.add.t1[simdat.modFit$model == 'bivar']
    SimRes$h2.add.t2_bivar[i] <- simdat.modFit$h2.add.t2[simdat.modFit$model == 'bivar']
    SimRes$h2.cov_slopeMod[i] <- simdat.modFit$h2.cov[simdat.modFit$model == 'slope_add']
    SimRes$rhog_bivar[i] <- simdat.modFit$rhog[simdat.modFit$model == 'bivar']
    SimRes$runtime_slope[i] <- simdat.modFit$runtime[simdat.modFit$model == 'slope_add']
    SimRes$runtime_bivar[i] <- simdat.modFit$runtime[simdat.modFit$model == 'bivar']
  }
  save(list = 'SimRes', file = resFile)
}else
  load(resFile)
```

## h2.add.t2
* x-axis = h2.add.t2 true
* y-axis = h2.add.t2 estimated

![](../figures/02-highN/plot_sim2_addT2-1.png)<!-- -->

## h2.cov & rhog
* x-axis = h2.add.t2 true
* y-axis = h2.cov (slope_add model), rhog (bivar model)

Dashed lines indicate simulated values for h2.cov and rhog

![](../figures/02-highN/plot_sim2_h2cov_rhog-1.png)<!-- -->

# Simulation 3: h2.cov + h2.add.t1 + h2.add.t2 + rhog


```r
#Parameters
h2_cov <- c(0.3, 0.6)
h2_add_t1 <- .3
h2_add_t2 <- c(1e-4, 0.1, 0.2, 0.3)
rhog <- .5
resFile <- '../results/02-highN/sim3.RData'

if(!file.exists(resFile)){
  grid <- expand.grid(n = as.integer(vals_n), sim = vals_sim, h2_add_t1 = h2_add_t1, h2_add_t2 = h2_add_t2, h2_cov = h2_cov, rhog = rhog)
  SimRes <- data.frame(grid, runtime_slope = NA, runtime_bivar = NA,
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
    simdat.modFit <- run_models(models, simdat$pheno, simdat$G)

    #Store results
    SimRes$h2.add.t2_slopeMod[i] <- simdat.modFit$h2.add.t2[simdat.modFit$model == 'slope_add']
    SimRes$h2.add.t1_bivar[i] <- simdat.modFit$h2.add.t1[simdat.modFit$model == 'bivar']
    SimRes$h2.add.t2_bivar[i] <- simdat.modFit$h2.add.t2[simdat.modFit$model == 'bivar']
    SimRes$h2.cov_slopeMod[i] <- simdat.modFit$h2.cov[simdat.modFit$model == 'slope_add']
    SimRes$rhog_bivar[i] <- simdat.modFit$rhog[simdat.modFit$model == 'bivar']
    SimRes$runtime_slope[i] <- simdat.modFit$runtime[simdat.modFit$model == 'slope_add']
    SimRes$runtime_bivar[i] <- simdat.modFit$runtime[simdat.modFit$model == 'bivar']
  }
  save(list = 'SimRes', file = resFile)
}else
  load(resFile)
```

## h2.add.t2
* x-axis = h2.add.t2 true
* y-axis = h2.add.t2 estimated

![](../figures/02-highN/plot_sim3_addT2-1.png)<!-- -->

## h2.cov & rhog
* x-axis = h2.add.t2 true
* y-axis = h2.cov (slope_add model), rhog (bivar model)

Dashed lines indicate simulated values for h2.cov and rhog

![](../figures/02-highN/plot_sim3_h2cov_rhog-1.png)<!-- -->

# Simulation 4: h2.cov + h2.add.t1 + h2.add.t2 + rhog

Same as simulation 3 but with rhog = 0.9


```r
#Parameters
h2_cov <- c(0.3, 0.6)
h2_add_t1 <- .3
h2_add_t2 <- c(1e-4, 0.1, 0.2, 0.3)
rhog <- .9
resFile <- '../results/02-highN/sim4.RData'

if(!file.exists(resFile)){
  grid <- expand.grid(n = as.integer(vals_n), sim = vals_sim, h2_add_t1 = h2_add_t1, h2_add_t2 = h2_add_t2, h2_cov = h2_cov, rhog = rhog)
  SimRes <- data.frame(grid, runtime_slope = NA, runtime_bivar = NA,
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
    simdat.modFit <- run_models(models, simdat$pheno, simdat$G)

    #Store results
    SimRes$h2.add.t2_slopeMod[i] <- simdat.modFit$h2.add.t2[simdat.modFit$model == 'slope_add']
    SimRes$h2.add.t1_bivar[i] <- simdat.modFit$h2.add.t1[simdat.modFit$model == 'bivar']
    SimRes$h2.add.t2_bivar[i] <- simdat.modFit$h2.add.t2[simdat.modFit$model == 'bivar']
    SimRes$h2.cov_slopeMod[i] <- simdat.modFit$h2.cov[simdat.modFit$model == 'slope_add']
    SimRes$rhog_bivar[i] <- simdat.modFit$rhog[simdat.modFit$model == 'bivar']
    SimRes$runtime_slope[i] <- simdat.modFit$runtime[simdat.modFit$model == 'slope_add']
    SimRes$runtime_bivar[i] <- simdat.modFit$runtime[simdat.modFit$model == 'bivar']
  }
  save(list = 'SimRes', file = resFile)
}else
  load(resFile)
```

## h2.add.t2
* x-axis = h2.add.t2 true
* y-axis = h2.add.t2 estimated

![](../figures/02-highN/plot_sim4_addT2-1.png)<!-- -->

## h2.cov & rhog
* x-axis = h2.add.t2 true
* y-axis = h2.cov (slope_add model), rhog (bivar model)

Dashed lines indicate simulated values for h2.cov and rhog

![](../figures/02-highN/plot_sim4_h2cov_rhog-1.png)<!-- -->

# Runtime
![](../figures/02-highN/plot_runtime-1.png)<!-- -->

Average runtimes (minutes):

| slope_add **n = 500**x10 | bivar **n = 500**x10 | slope_add **n = 1000**x10 | bivar **n = 1000**x10 |
|--------------|--------------|--------------|--------------|
| 18.008 | 95.539| 70.295 | 308.303 |
