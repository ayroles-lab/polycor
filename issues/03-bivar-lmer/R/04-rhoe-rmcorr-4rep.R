### inc 
library(tidyverse)

library(rmcorr)
library(lme4)

# change to your path here, e.g. based on `nodemame`
path <- switch(Sys.info()[['nodename']],
  "~/git/hemostat/polycor/" )
  
src <- file.path(path, "relmer.R")
source(src)

### data
#> ?gilden201
# A dataset containing four repeated measurements of reaction time
# (RT) and accuracy from eleven subjects in a visual search
# experiment. Each measurement is the mean RT and accuracy from a
# block of 288 search trials. blocks of visual search, for eleven
# subjects
     
data(gilden2010, package = "rmcorr")
dat <- as_tibble(gilden2010) %>% 
  mutate(
    sub = as.factor(sub), 
    block = as.factor(block),
    obs = as.factor(1:n()))

bdat <- gather(dat, tname, trait, rt, acc)

### lme4
mod <- relmer(trait ~ tname - 1 + (0 + tname | sub) + (0 + tname | obs), data = bdat, weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE)

Ztlist <- getME(mod, "Ztlist")
lapply(Ztlist, dim)

stop()

### rmcorr
c1 <- rmcorr(participant = Participant, measure1 = Age, measure2 = Volume, dataset = dat)

### lm: thec corr. is the same
lm(scale(Volume) ~ scale(Age), dat) %>% coef
lm(scale(Age) ~ scale(Volume), dat) %>% coef

### lmer: the corr. is different
relmatLmer(scale(Volume) ~ scale(Age) + (1|Participant), dat) %>% fixef
relmatLmer(scale(Age) ~ scale(Volume) + (1|Participant), dat) %>% fixef


