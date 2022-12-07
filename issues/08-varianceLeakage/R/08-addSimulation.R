library(tidyverse)
library(magrittr)
library(MASS)
library(lme4)
library(lme4qtl)
source('polycor/simFunctions.R')
source('polycor/relmer.R')

#Random slope models fitted to additive data. Are there any leakage from purely additive genetic effects to the "h2 of the slope"?
#I use the simulation function pleiotropy.simulatePairKin with rhog = 0. Both traits have the same additive h2


#### Simulation ####
vals_h2_add <- c(0, .2, .4, .6, .8)
vals_sim <- 1:5
vals_n <- c(100, 500, 1e3)
nrep <- 5
grid <- expand.grid(n = as.integer(vals_n), h2_add = vals_h2_add, sim = vals_sim)
addSimRes <- data.frame(grid, h2.slope = NA, h2.slope_withAdd = NA, h2.add_trait2_slopeMod = NA, h2.add_trait1 = NA, h2.add_trait2 = NA)

ptm <- Sys.time()
for (i in 1:nrow(grid)) {
  n <- grid$n[i]
  h2_add <- grid$h2_add[i]
  sim <- grid$sim[i]
  print(grid[i, ])
  
  simdat <- pleiotropy.simulatePairKin(n_ind = n, n_perInd = nrep, h2 = h2_add, overlap = 0)
  dat <- simdat$pheno %>% mutate(id = as.character(id), rid = id)
  G <- simdat$G
  
  #Random slope
  mod <- relmer(trait2 ~ obs + (0 + trait1|id), data = dat, relmat = list(id = G))
  addSimRes$h2.slope[i] <- VarProp(mod) %>% as.data.frame %>% filter(grp == "id") %$% prop
  
  #Random slope with add
  mod <- relmer(trait2 ~ obs + (1|rid) + (0 + trait1|id), data = dat, relmat = list(id = G, rid = G))
  addSimRes$h2.slope_withAdd[i] <- VarProp(mod) %>% as.data.frame %>% filter(grp == "id") %$% prop
  addSimRes$h2.add_trait2_slopeMod[i] <- VarProp(mod) %>% as.data.frame %>% filter(grp == "rid") %$% prop
  
  #Additive model
  #Trait 1
  mod <- relmer(trait1 ~ obs + (1|id), data = dat, relmat = list(id = G))
  addSimRes$h2.add_trait1[i] <- VarProp(mod) %>% as.data.frame %>% filter(grp == "id") %$% prop
  #Trait 2
  mod <- relmer(trait2 ~ obs + (1|id), data = dat, relmat = list(id = G))
  addSimRes$h2.add_trait2[i] <- VarProp(mod) %>% as.data.frame %>% filter(grp == "id") %$% prop
}
Sys.time() - ptm
save(addSimRes, file = paste('polycor/issues/08-varianceLeakage/results/190128_addSim_nRep', nrep, '.RData', sep = ''))



#### Plot ####
source('polycor/plotFunctions.R')
png('polycor/issues/08-varianceLeakage/figures/addSim_randSlope_5reps_h2slope_v2.png', width = 1200, height = 800)
par(mfrow = c(1,3))
incl <- addSimRes$n == 100
plot.contBoxes(x = addSimRes$h2_add[incl], y = addSimRes[incl, 4:5], legend.text = c('slope', 'slope_add'), legend.title = 'Models', ylab = 'Estimated h2 slope', legend.cex = 3, spacer2 = 2, ylim = c(0,.7))
title('n = 50 x 2', cex.main = 3)
incl <- addSimRes$n == 500
plot.contBoxes(x = addSimRes$h2_add[incl], y = addSimRes[incl, 4:5], legend.text = c('slope', 'slope_add'), legend.title = 'Models', ylab = 'Estimated h2 slope', legend.cex = 3, spacer2 = 2, ylim = c(0,.7))
title('n = 250 x 2', cex.main = 3)
incl <- addSimRes$n == 1000
plot.contBoxes(x = addSimRes$h2_add[incl], y = addSimRes[incl, 4:5], legend.text = c('slope', 'slope_add'), legend.title = 'Models', ylab = 'Estimated h2 slope', legend.cex = 3, spacer2 = 2, ylim = c(0,.7))
title('n = 500 x 2', cex.main = 3)
dev.off()


