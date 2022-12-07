### inc
library(tidyverse)
library(magrittr)

library(Matrix)
library(lme4)

library(devtools)
load_all("~/git/variani/wlm/")

# change to your path here, e.g. based on `nodemame`
path <- switch(Sys.info()[['nodename']],
  "~/git/hemostat/polycor/", 
  "gen-sflaptop" = "/media/sf9/ExtraDrive1/Dropbox/Research_170228/covQTL/bin/polycor/")
  
src <- file.path(path, "relmer.R")
source(src)

### simulate data
N <- 500; M <- 200; h2 <- 0.8
  
Zg <- sapply(1:M, function(i) rbinom(N, 2, 0.5)) # allele freq. = 0.5

col_means <- colMeans(Zg, na.rm = TRUE)
col_freq <- col_means / 2  # col_means = 2 * col_freq
col_sd <- sqrt(2 * col_freq * (1 - col_freq))

Z <- sweep(Zg, 2, col_means, "-")
Z <- sweep(Z, 2, col_sd , "/")

b <- rnorm(M, 0, sqrt(h2/M))
y <- Z %*% b + rnorm(N, 0, sqrt(1 - h2))

ids <- paste0("id", seq(N))
dat <- data.frame(y = y, id = ids)
rownames(Z) <- ids

### (1) fit  LMM by an external function from `wlm` package
#m0 <- lmm1lr(y ~ 1, dat, zmat = Z/sqrt(M))

### (2) fit LMM using `relmat`
G <- tcrossprod(Z/sqrt(M))
m1 <- relmer(y ~ (1|id), dat, relmat = list(id = G))
VarProp(m1)$prop[1] #[1] 0.7839681

### (3) fit LMM using `relfac`
m2 <- relmer(y ~ (1|id), dat, relfac = list(id = Z/sqrt(M)))
VarProp(m2)$prop[1] 

### plot Z'Z matrices
# - !NB! the G = Z'Z matrices in two models are different
varcov(m2, "id", residual = FALSE, scaled = FALSE) %>% .[1:10, 1:10]  %>% image
varcov(m1, "id", residual = FALSE, scaled = FALSE) %>% .[1:10, 1:10]  %>% image


