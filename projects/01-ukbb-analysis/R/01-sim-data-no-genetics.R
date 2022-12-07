### inc
library(MASS)

library(dplyr)
library(magrittr)
library(tidyverse)

library(lme4)

# change to your path here, e.g. based on `nodemame`
path <- switch(Sys.info()[['nodename']],
  "~/git/hemostat/polycor/" )
  
src <- file.path(path, "relmer.R")
source(src)

### par
N <- 1e3
K <- 2
R <- 5

b <- 0.5

var_y <- 1
sd_y <- sqrt(var_y)
var_rep <- 0.5
sd_rep <- sqrt(var_rep)

### data simulations
set.seed(1)

# v1
y1 <- 1 + rnorm(N, sd = sd_y)
y2 <- 2 + b*y1 + rnorm(N, sd = sd_y * sqrt(1 - b^2))

# v2
#M <- matrix(b, K, K)
#diag(M) <- 1
#Y <- mvrnorm(N, rep(0, 2), M, empirical = TRUE)
#y1 <- 1 + Y[, 1] #+ rnorm(N, sd = sd_y)
#y2 <- 2 + Y[, 2] #+ rnorm(N, sd = sd_y)

x1 <- rnorm(N)

dat <- data_frame(id = rep(1:N, each = R), rep = rep(1:R, N),
    x1 = rep(x1, each = R),
    t1 = y1 + rnorm(N*R, sd = var_rep), 
    t2 = y2 + rnorm(N*R, sd = var_rep)) %>%
  mutate(obs = seq(1, n())) %>%
  mutate(id = as.factor(id), rep = as.factor(rep), obs = as.factor(obs))

# check that noise of repeats masks the corr.
cor(y1, y2) # 0.5023334
with(dat, cor(t1, t2)) # 0.2642627

### models
bdat <- gather(dat, tname, tvalue, t1, t2)

bmod <- relmer(tvalue ~ -1 + tname + (0 + tname | id) + (0 + tname | obs), 
  bdat, weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE)
