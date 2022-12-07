### inc
library(dplyr)
library(magrittr)
library(tidyverse)

library(lme4)

library(microbenchmark)
library(peakRAM)

# change to your path here, e.g. based on `nodemame`
path <- switch(Sys.info()[['nodename']],
  "~/git/hemostat/polycor/" )
  
src <- file.path(path, "relmer.R")
source(src)

### data simulation function
sim_dat <- function(n, g = 100)
{
  m <- n/g
  
  data_frame(gr = rep(1:g, each = m)) %>%
    mutate(
      x1 = rnorm(n), x2 = rnorm(n),
      y = rnorm(m)[gr] + rnorm(n),
      gr = as.factor(gr))
}

### run for n = 1000
n <- 500e3
dat <- sim_dat(n)

cpu <- microbenchmark(
  lmer(y ~ (1|gr), dat),
  lmer(y ~ (1|gr) + x1, dat),
  lmer(y ~ (1|gr) + x1 + x2, dat),
  times = 5)
  
autoplot(out)

mem <- peakRAM(
  lmer(y ~ 0 + (1|gr), dat),
  lmer(y ~ (1|gr), dat),
  lmer(y ~ (1|gr) + x1, dat),
  lmer(y ~ (1|gr) + x1 + x2, dat))
