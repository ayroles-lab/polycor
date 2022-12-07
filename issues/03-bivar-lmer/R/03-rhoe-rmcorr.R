# https://stats.stackexchange.com/questions/13166/rs-lmer-cheat-sheet
# http://lme4.r-forge.r-project.org/slides/2011-03-16-Amsterdam/2Longitudinal.pdf
# Multivariate analysis with mixed modeling tools in R: http://rpubs.com/bbolker/3336
# BLUP formula: http://variani.github.io/talks/2014/01-mixed-models-qtl/#19
# Explained rep. meas.: https://www.theanalysisfactor.com/the-intraclass-correlation-coefficient-in-mixed-models/
# Diff. between VarCor and ranef: https://stats.stackexchange.com/questions/153253/random-slope-and-intercept-correlation-not-consistent-in-output-vs-manual-calcu
# Removing random correlation parameter for categorical variable in lmer: https://stat.ethz.ch/pipermail/r-sig-mixed-models/2016q1/024338.html
# Example datasets with ind.-level cor.: https://cran.r-project.org/web/packages/rmcorr/vignettes/Combined_paper_figures.html

# Tutorials on mixed models and lme4
# - https://rpsychologist.com/r-guide-longitudinal-lme-lmer#conditional-growth-model-dropping-intercept-slope-covariance
# - https://stats.stackexchange.com/questions/320978/understanding-and-coding-random-intercept-correlation-lmer

### inc 
library(tidyverse)

library(Matrix)
library(MASS)

library(rmcorr)
library(nlme)
library(lme4)

### data
data(raz2005, package = "rmcorr")
dat <- as_tibble(raz2005) %>% 
  mutate(
    Participant = as.factor(Participant), 
    Time = as.factor(Time),
    Obs = as.factor(1:n()))

bdat <- gather(dat, tname, trait, Age, Volume)

#m <- lmer(trait ~ (0 + tname|rid), bdat)

#### nlme
#ctrl <- lmeControl(opt = "optim")
#m1 <- lme(trait ~ tname - 1, random = ~tname - 1 | Participant, data = bdat, control = ctrl)
#m2 <- lme(trait ~ tname - 1, random = ~tname - 1 | Participant, weights = varIdent(form = ~1 | tname), data = bdat, control = ctrl)
#m3 <- lme(trait ~ tname - 1, random = ~tname - 1 | Participant, weights = varIdent(form = ~1 | tname), correlation = corSymm(form = ~1 | Participant/Obs), data = bdat, control = ctrl)
#m3 <- lme(trait ~ tname - 1, random = ~tname - 1 | id, weights = varIdent(form = ~1 | trait), correlation = corSymm(form = ~1 | id/obs), data =

### lme4
m4 <- relmatLmer(trait ~ tname - 1 + (0 + tname | Participant) + (0 + tname | Obs), data = bdat, weights = rep(1e10, nrow(bdat)))
#m5 <- relmatLmer(trait ~ tname - 1 + (0 + tname | Participant) + (0 + tname | Obs), data = bdat, weights = rep(1e10, nrow(bdat)))

Ztlist <- getME(m4, "Ztlist")

### rmcorr
c1 <- rmcorr(participant = Participant, measure1 = Age, measure2 = Volume, dataset = dat)

### lm: thec corr. is the same
lm(scale(Volume) ~ scale(Age), dat) %>% coef
lm(scale(Age) ~ scale(Volume), dat) %>% coef

### lmer: the corr. is different
relmatLmer(scale(Volume) ~ scale(Age) + (1|Participant), dat) %>% fixef
relmatLmer(scale(Age) ~ scale(Volume) + (1|Participant), dat) %>% fixef


