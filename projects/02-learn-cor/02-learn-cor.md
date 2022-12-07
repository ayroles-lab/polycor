---
title: "Learning basics about correlations"
author: "Andrey Ziyatdinov"
date: "2019-02-01"
output:
  html_document:
    theme: united
    toc: true
    keep_md: true
---







# Include custom code


```r
library(Matrix)
library(MASS)
library(lme4)

library(tidyverse)
library(magrittr)

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

path <- switch(Sys.info()[['nodename']],
  "~/git/hemostat/polycor/" )
  
file.path(path, "relmer.R") %>% source
file.path(path, "simFunctions.R") %>% source
```

# Two pleitropic traits

We quickly run a single round of data simulation and conclusions differ from simulation to simulation.
Multiple runs and larger smaple sizes are needed to evaluate the accurcy of estimates.

## Two pleitropic traits + bi-variate model

Given the two traits are simulated under a bi-variate model with a substantial portion of genetic covariance,
what estimates of heritabilities are expected from uni-variate models?

Simulate two pleitropic traits:
 

```r
set.seed(1)
simdat <- pleiotropy.simulatePairKin(n_ind = 500, n_perInd = 1, 
  h2 = 0.5, overlap = 0.9)

dat <- simdat$pheno %>% as_tibble %>%
  mutate(
    trait1sc = as.numeric(scale(trait1)), trait2sc = as.numeric(scale(trait2)),
    id = as.factor(as.character(id)), rid = id)
bdat <- gather(dat, tname, tvalue, trait1sc, trait2sc)

G <- simdat$G
```

Fit the standard bi-variate model (weights):



```r
mod <- relmer(tvalue ~ -1 + tname + (0 + tname|id) + (0 + tname|rid), 
  bdat, relmat = list(id = G), 
  weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE)
mod
```

```
Linear mixed model fit by REML ['lmerMod']
Formula: tvalue ~ -1 + tname + (0 + tname | id) + (0 + tname | rid)
   Data: bdat
Weights: rep(1e+10, nrow(bdat))
REML criterion at convergence: 2674.111
Random effects:
 Groups   Name          Std.Dev. Corr 
 id       tnametrait1sc 0.7964        
          tnametrait2sc 0.7025   0.84 
 rid      tnametrait1sc 0.6217        
          tnametrait2sc 0.7265   -0.13
 Residual               0.6952        
Number of obs: 1000, groups:  id, 500; rid, 500
Fixed Effects:
tnametrait1sc  tnametrait2sc  
     -0.01571       -0.01471  
```


Fit the standard bi-variate model (no weights/redundant residual random effects):



```r
mod <- relmer(tvalue ~ -1 + tname + (0 + tname|id) + (0 + tname|rid), 
  bdat, relmat = list(id = G), 
  calc.derivs = FALSE)
mod
```

```
Linear mixed model fit by REML ['lmerMod']
Formula: tvalue ~ -1 + tname + (0 + tname | id) + (0 + tname | rid)
   Data: bdat
REML criterion at convergence: 2674.109
Random effects:
 Groups   Name          Std.Dev. Corr 
 id       tnametrait1sc 0.7992        
          tnametrait2sc 0.7041   0.84 
 rid      tnametrait1sc 0.3013        
          tnametrait2sc 0.4833   -0.42
 Residual               0.5410        
Number of obs: 1000, groups:  id, 500; rid, 500
Fixed Effects:
tnametrait1sc  tnametrait2sc  
     -0.01590       -0.01486  
```

Fit uni-variate models for each trait and estimates of heritability:


```r
mod1 <- relmer(trait1 ~ (1|id), dat, relmat = list(id = G)) 
mod1 %>% VarProp %>% .[1, "prop"]
```

```
[1] 0.6097547
```

```r
mod2 <- relmer(trait2 ~ (1|id), dat, relmat = list(id = G)) 
mod2 %>% VarProp %>% .[1, "prop"]
```

```
[1] 0.4845626
```


## Two random-sloped traits + bi-variate model

We simulate two traits under a extreme scenario for bi-variate model
(random-slope data simulation model):

- no additive variance for any trait;
- random-slope genetic covariance.

What estimates of heritability are produced by uni-variate models?

Simulate two traits:


```r
set.seed(1)
simdat <- genCov.simulatePairKin(n_ind = 500, n_perInd = 1, 
  h2.cov = 0.9, h2.add = 0)

dat <- simdat$pheno %>% as_tibble %>%
  mutate(
    trait1sc = as.numeric(scale(trait1)), trait2sc = as.numeric(scale(trait2)),
    id = as.factor(as.character(id)), rid = id)
bdat <- gather(dat, tname, tvalue, trait1sc, trait2sc)

G <- simdat$G
```

Fit the standard bi-variate model:

- in some runs: the overall genetic component of variance is close to zero;
- in other runs: the trait2-specific heritabilty seems different from zero.


```r
mod <- relmer(tvalue ~ -1 + tname + (0 + tname|id) + (0 + tname|rid), 
  bdat, relmat = list(id = G), 
  calc.derivs = FALSE)
mod
```

```
Linear mixed model fit by REML ['lmerMod']
Formula: tvalue ~ -1 + tname + (0 + tname | id) + (0 + tname | rid)
   Data: bdat
REML criterion at convergence: 2831.881
Random effects:
 Groups   Name          Std.Dev.  Corr 
 id       tnametrait1sc 2.681e-06      
          tnametrait2sc 5.167e-06 1.00 
 rid      tnametrait1sc 6.670e-01      
          tnametrait2sc 6.670e-01 -0.36
 Residual               7.450e-01      
Number of obs: 1000, groups:  id, 500; rid, 500
Fixed Effects:
tnametrait1sc  tnametrait2sc  
   -2.914e-13     -5.619e-13  
```

Fit uni-variate models for each trait and estimates of heritability:

- the estimates of heritability for each trait is zero.


```r
mod1 <- relmer(trait1 ~ (1|id), dat, relmat = list(id = G)) 
mod1 %>% VarProp %>% .[1, "prop"]
```

```
[1] 0
```

```r
mod2 <- relmer(trait2 ~ (1|id), dat, relmat = list(id = G)) 
mod2 %>% VarProp %>% .[1, "prop"]
```

```
[1] 1.332829e-15
```

# Population vs. within individual correlation

## Two pleitropic traits + repeats 


```r
set.seed(1)
simdat <- pleiotropy.simulatePairKin(n_ind = 10, n_perInd = 100, 
  h2 = 0.95, overlap = 0.95)

dat <- simdat$pheno %>% as_tibble %>%
  mutate(
    trait1sc = as.numeric(scale(trait1)), trait2sc = as.numeric(scale(trait2)))

ggplot(dat, aes(trait1sc, trait2sc)) + geom_point() + geom_smooth(method = "lm", color = "red") + geom_hline(yintercept = 0, linetype = 3) + geom_abline(linetype = 3) + labs(title = "Population-level corr.")
```

![](figures/bivar_rep_cor-1.png)<!-- -->

```r
ggplot(filter(dat, id %in% 1:4), aes(trait1sc, trait2sc)) + geom_point() + geom_smooth(method = "lm", color = "blue") + geom_hline(yintercept = 0, linetype = 3) + geom_abline(linetype = 3) + facet_wrap(~ id)  + labs(title = "Individual-level corr.")
```

![](figures/bivar_rep_cor-2.png)<!-- -->

## Two random-slope traits + repeats 


```r
set.seed(1)
simdat <- genCov.simulatePairKin(n_ind = 10, n_perInd = 100, 
  h2.cov = 0.95, h2.add = 0)

dat <- simdat$pheno %>% as_tibble %>%
  mutate(
    trait1sc = as.numeric(scale(trait1)), trait2sc = as.numeric(scale(trait2)))

ggplot(dat, aes(trait1sc, trait2sc)) + geom_point() + geom_smooth(method = "lm", color = "red") + geom_hline(yintercept = 0, linetype = 3) + geom_abline(linetype = 3) + labs(title = "Population-level corr.")
```

![](figures/slope_rep_cor-1.png)<!-- -->

```r
ggplot(filter(dat, id %in% 1:4), aes(trait1sc, trait2sc)) + geom_point() + geom_smooth(method = "lm", color = "blue") + geom_hline(yintercept = 0, linetype = 3) + geom_abline(linetype = 3) + facet_wrap(~ id) + labs(title = "Individual-level corr.")
```

![](figures/slope_rep_cor-2.png)<!-- -->


```r
simdat$pheno$g.cov %>% unique %>% head(4) %>% round(2)
```

```
[1]  0.91 -0.24  0.87  0.05
```

