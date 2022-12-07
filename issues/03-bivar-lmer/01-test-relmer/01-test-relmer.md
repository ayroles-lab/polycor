---
title: "Testing relmer function"
author: "Andrey Ziyatdinov"
date: "2019-01-04"
output:
  html_document:
    theme: united
    toc: true
    number_sections: true
    keep_md: true
---






```r
library(magrittr)
library(dplyr)
library(tidyr)

library(MASS)

library(Matrix)
library(lme4)

library(lme4qtl)
library(solarius)
```

# Load `relmer` function


```r
# change to your path here, e.g. based on `nodemame`
path <- switch(Sys.info()[['nodename']],
  "~/git/hemostat/polycor/" )
  
src <- file.path(path, "relmer.R")
source(src)
```

# Test model fitting on examples

## `sleepstudy`: lme4 vs relmer


```r
m1 <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
(m2 <- relmer(Reaction ~ Days + (Days | Subject), sleepstudy))
```

```
Linear mixed model fit by REML ['lmerMod']
Formula: Reaction ~ Days + (Days | Subject)
   Data: sleepstudy
REML criterion at convergence: 1743.628
Random effects:
 Groups   Name        Std.Dev. Corr
 Subject  (Intercept) 24.740       
          Days         5.922   0.07
 Residual             25.592       
Number of obs: 180, groups:  Subject, 18
Fixed Effects:
(Intercept)         Days  
     251.41        10.47  
```

```r
stopifnot(all.equal(logLik(m1), logLik(m2)))
stopifnot(all.equal(fixef(m1), fixef(m2)))
```

## Univariate polygenic: lme4qtl vs relmer


```r
data(dat40)

m1 <- lme4qtl::relmatLmer(trait1 ~ AGE + SEX + (1|ID), dat40, relmat = list(ID = kin2))
(m2 <- relmer(trait1 ~ AGE + SEX + (1|ID), dat40, relmat = list(ID = kin2)))
```

```
Linear mixed model fit by REML ['lmerMod']
Formula: trait1 ~ AGE + SEX + (1 | ID)
   Data: dat40
REML criterion at convergence: 963.3853
Random effects:
 Groups   Name        Std.Dev.
 ID       (Intercept) 2.2988  
 Residual             0.7856  
Number of obs: 224, groups:  ID, 224
Fixed Effects:
(Intercept)          AGE         SEX2  
   7.563248     0.008314    -0.364197  
```

```r
stopifnot(all.equal(logLik(m1), logLik(m2)))
stopifnot(all.equal(fixef(m1), fixef(m2)))
```

## Univariate sex-specifc polygenic: lme4qtl vs relmer


```r
data(dat40)

m1 <- lme4qtl::relmatLmer(trait1 ~ AGE + SEX + (1|ID), dat40, relmat = list(ID = kin2))
(m2 <- relmer(trait1 ~ AGE + SEX + (1|ID), dat40, relmat = list(ID = kin2)))
```

```
Linear mixed model fit by REML ['lmerMod']
Formula: trait1 ~ AGE + SEX + (1 | ID)
   Data: dat40
REML criterion at convergence: 963.3853
Random effects:
 Groups   Name        Std.Dev.
 ID       (Intercept) 2.2988  
 Residual             0.7856  
Number of obs: 224, groups:  ID, 224
Fixed Effects:
(Intercept)          AGE         SEX2  
   7.563248     0.008314    -0.364197  
```

```r
stopifnot(all.equal(logLik(m1), logLik(m2)))
stopifnot(all.equal(fixef(m1), fixef(m2)))
```

## Bivariate polygenic: solarius vs relmer


```r
data(dat40)
# no missing data is a current requirement of `relmer` for nested model
#dat <- na.omit(dat40)
bdat <- within(dat, RID <- ID)
bdat <- tidyr::gather(bdat, tname, trait, trait1, trait2)
wts <- rep(1e10, nrow(bdat))

(m1 <- solarius::solarPolygenic(trait1 + trait2 ~ 1, dat40))
```

```

Call: solarius::solarPolygenic(formula = trait1 + trait2 ~ 1, data = dat40)

File polygenic.out:
	Pedigree:    dat.ped 
	Phenotypes:  dat.phe 
	Trait:       trait1 trait2         Individuals:  234 
 
			 H2r(trait1) is 0.8986090   
	       H2r(trait1) Std. Error:  0.0807720 
 
			 H2r(trait2) is 0.7575300   
	       H2r(trait2) Std. Error:  0.0993358 
 
			 RhoE is 0.3788398   
	       RhoE Std. Error:  0.2825170 
 
			 RhoG is 0.9140881 
	       RhoG Std. Error:  0.0436821 
 
	       Derived Estimate of RhoP is 0.8135769 
 
 
	Loglikelihoods and chi's are in trait1.trait2/polygenic.logs.out 
	Best model is named poly and null0 
	Final models are named poly, spor 
```

```r
(m2 <- relmer(trait ~ (0 + tname|ID) + (0 + dummy(tname)|RID), bdat, relmat = list(ID = kin2)))
```

```
Linear mixed model fit by REML ['lmerMod']
Formula: trait ~ (0 + tname | ID) + (0 + dummy(tname) | RID)
   Data: bdat
REML criterion at convergence: 1737.543
Random effects:
 Groups   Name         Std.Dev. Corr
 ID       tnametrait1  2.3230       
          tnametrait2  2.5550   0.91
 RID      dummy(tname) 0.5624       
 Residual              0.7580       
Number of obs: 438, groups:  ID, 219; RID, 219
Fixed Effects:
(Intercept)  
      7.707  
```

```r
(m3 <- relmer(trait ~ (0 + tname|ID) + (0 + tname|RID), bdat, weights = wts, relmat = list(ID = kin2), calc.derivs = FALSE))
```

```
Linear mixed model fit by REML ['lmerMod']
Formula: trait ~ (0 + tname | ID) + (0 + tname | RID)
   Data: bdat
Weights: wts
REML criterion at convergence: 1737.043
Random effects:
 Groups   Name        Std.Dev. Corr
 ID       tnametrait1 2.2040       
          tnametrait2 2.4519   0.90
 RID      tnametrait1 0.9241       
          tnametrait2 1.0724   0.30
 Residual             1.2306       
Number of obs: 438, groups:  ID, 219; RID, 219
Fixed Effects:
(Intercept)  
      7.694  
```

