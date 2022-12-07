# https://stats.stackexchange.com/questions/13166/rs-lmer-cheat-sheet
# http://lme4.r-forge.r-project.org/slides/2011-03-16-Amsterdam/2Longitudinal.pdf
# http://rpubs.com/bbolker/3336

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

### par
e1 <- 0.5
e2 <- 0.5
re <- 0.75

sim_ids <- function(N) paste0("ID", 1:N)

sim_bivar <- function(e1, e2, re, N, rep = 1) 
{
  # V = [I]
  e12 <- re*sqrt(e1)*sqrt(e2)
  Sigma_resid <- matrix(c(e1, e12, e12, e2), 2, 2)
  V_resid <- kronecker(Sigma_resid, Diagonal(N))

  V <- V_resid
  
  lapply(seq(1, rep), function(r) {
    y12 <- mvrnorm(1, rep(0, 2*N), V)
    y1 <- y12[seq(1, N)]
    y2 <- y12[N + seq(1, N)]
    
    ids <- sim_ids(N)
    reps <- paste0("rep", r)
    data_frame(id = ids, rid = ids, y1 = y1, y2 = y2, rep = reps)
  }) %>% 
  bind_rows %>%
  mutate(obs = paste0("obs", seq(n())))
}

# Start with no rep.
N <- 1000
dat <- sim_bivar(e1, e2, re, N)
bdat <- gather(dat, tname, tvalue, y1, y2) 

m1 <- relmer(tvalue ~ tname - 1 + (0 + tname | id), data = bdat, weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE)


N <- 500
R <- 10
dat <- sim_bivar(e1, e2, re, N, R)
bdat <- gather(dat, tname, tvalue, y1, y2) %>%
  mutate(tindex = (as.integer(as.factor(tname)) - 1))

m2 <- relmer(tvalue ~ tname - 1 + (0 + tname | id) + (0 + tname | obs), data = bdat, weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE)

m3 <- relmer(tvalue ~ tname - 1 + (0 + tindex || id) + (1 + tname | obs), data = bdat, weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE)

N <- 10
R <- 500
dat <- sim_bivar(e1, e2, re, N, R)
bdat <- gather(dat, tname, tvalue, y1, y2) %>%
  mutate(tindex = (as.integer(as.factor(tname)) - 1))

m4 <- relmer(tvalue ~ tname - 1 + (0 + tname | id) + (0 + tname | obs), data = bdat, weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE)

m5 <- relmer(tvalue ~ tname - 1 + (0 + tindex || id) + (1 + tname | obs), data = bdat, weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE)

stop()


#### nlme
ctrl <- lmeControl(opt = "optim")

m1 <- lme(trait ~ tname - 1, random = ~tname - 1 | id, data = bdat, control = ctrl)
m2 <- lme(trait ~ tname - 1, random = ~tname - 1 | id, weights = varIdent(form = ~1 | trait), data = bdat, control = ctrl)
#m3 <- lme(trait ~ tname - 1, random = ~tname - 1 | id, weights = varIdent(form = ~1 | trait), correlation = corSymm(form = ~1 | id/obs), data =

