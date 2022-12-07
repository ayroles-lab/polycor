library(tidyverse)
library(magrittr)
library(MASS)
library(lme4)
library(lme4qtl)
source('polycor/simFunctions.R')
source('polycor/relmer.R')

#Random slope models fitted to pleiotropic data. Are there any leakage from genetic correlation to the "h2 of the slope"?

#### Simulation ####
vals_h2_add <- .75
vals_sim <- 1:5
vals_n <- c(100, 500, 1e3)
vals_rhog <- c(0.1, 0.5, 0.9)
nrep <- c(2, 5, 10, 20)
grid <- expand.grid(n = as.integer(vals_n), h2_add = vals_h2_add, sim = vals_sim, rhog = vals_rhog, nrep = nrep)
pleioSimRes <- data.frame(grid, h2.slope = NA, h2.slope_withAdd = NA, h2.add_trait2_slopeMod = NA, h2.add_trait1 = NA, h2.add_trait2 = NA,
                          rhog_hat = NA)

ptm <- Sys.time()
for (i in 1:nrow(grid)) {
  n <- grid$n[i]
  h2_add <- grid$h2_add[i]
  sim <- grid$sim[i]
  rhog <- grid$rhog[i]
  nrep <- grid$nrep[i]
  print(grid[i, ])
  
  simdat <- pleiotropy.simulatePairKin(n_ind = n, n_perInd = nrep, h2 = h2_add, overlap = rhog)
  dat <- simdat$pheno %>% mutate(id = as.character(id), rid = id)
  G <- simdat$G
  
  #Random slope
  mod <- relmer(trait2 ~ obs + (0 + trait1|id), data = dat, relmat = list(id = G))
  pleioSimRes$h2.slope[i] <- VarProp(mod) %>% as.data.frame %>% filter(grp == "id") %$% prop
  
  #Random slope with add
  mod <- relmer(trait2 ~ obs + (1|rid) + (0 + trait1|id), data = dat, relmat = list(id = G, rid = G))
  pleioSimRes$h2.slope_withAdd[i] <- VarProp(mod) %>% as.data.frame %>% filter(grp == "id") %$% prop
  pleioSimRes$h2.add_trait2_slopeMod[i] <- VarProp(mod) %>% as.data.frame %>% filter(grp == "rid") %$% prop
  
  #Additive model
  #Trait 1
  mod <- relmer(trait1 ~ obs + (1|id), data = dat, relmat = list(id = G))
  pleioSimRes$h2.add_trait1[i] <- VarProp(mod) %>% as.data.frame %>% filter(grp == "id") %$% prop
  #Trait 2
  mod <- relmer(trait2 ~ obs + (1|id), data = dat, relmat = list(id = G))
  pleioSimRes$h2.add_trait2[i] <- VarProp(mod) %>% as.data.frame %>% filter(grp == "id") %$% prop
  
  
  #Bi-variate model
  bdat <- gather(dat, tname, tvalue, trait1, trait2)
  mod <- relmer(tvalue ~ -1 + tname + obs + (0 + tname | id) + (0 + tname | rid), 
                data = bdat, relmat = list(id = G), calc.derivs = FALSE)
  vf <- VarProp(mod) %>% as.data.frame %>% filter(grp == "id")
  cov <- filter(vf, grp == "id" & var1 == "tnametrait1" & var2 == "tnametrait2") %$% vcov
  pleioSimRes$rhog_hat[i] <- cov / sqrt(vf$vcov[1]*vf$vcov[2])
  
}
Sys.time() - ptm #~32min
save(pleioSimRes, file = 'polycor/issues/08-varianceLeakage/results/190129_pleioSim.RData')


#### Plot ####
source('polycor/plotFunctions.R')

#nRep = 2
pleioSimRes.nrep2 <- pleioSimRes[pleioSimRes$nrep == 2, ]
png('polycor/issues/08-varianceLeakage/figures/pleioSim_randSlope_2reps_h2slope.png', width = 1200, height = 800)
par(mfrow = c(1,3))
incl <- pleioSimRes.nrep2$n == 100
plot.contBoxes(x = pleioSimRes.nrep2$rhog[incl], y = pleioSimRes.nrep2[incl, 6:7], legend.text = c('slope', 'slope_add'), legend.title = 'Models', ylab = 'Estimated h2 slope', xlab = 'Genetic cor', legend.cex = 3, spacer2 = 2, ylim = c(0,.9))
title('n = 100 x 2', cex.main = 3)
incl <- pleioSimRes.nrep2$n == 500
plot.contBoxes(x = pleioSimRes.nrep2$rhog[incl], y = pleioSimRes.nrep2[incl, 6:7], legend.text = c('slope', 'slope_add'), legend.title = 'Models', ylab = 'Estimated h2 slope', xlab = 'Genetic cor', legend.cex = 3, spacer2 = 2, ylim = c(0,.9))
title('n = 500 x 2', cex.main = 3)
incl <- pleioSimRes.nrep2$n == 1000
plot.contBoxes(x = pleioSimRes.nrep2$rhog[incl], y = pleioSimRes.nrep2[incl, 6:7], legend.text = c('slope', 'slope_add'), legend.title = 'Models', ylab = 'Estimated h2 slope', xlab = 'Genetic cor', legend.cex = 3, spacer2 = 2, ylim = c(0,.9))
title('n = 1000 x 2', cex.main = 3)
dev.off()

#nRep = 5
pleioSimRes.nrep5 <- pleioSimRes[pleioSimRes$nrep == 5, ]
png('polycor/issues/08-varianceLeakage/figures/pleioSim_randSlope_5reps_h2slope.png', width = 1200, height = 800)
par(mfrow = c(1,3))
incl <- pleioSimRes.nrep5$n == 100
plot.contBoxes(x = pleioSimRes.nrep5$rhog[incl], y = pleioSimRes.nrep5[incl, 6:7], legend.text = c('slope', 'slope_add'), legend.title = 'Models', ylab = 'Estimated h2 slope', xlab = 'Genetic cor', legend.cex = 3, spacer2 = 2, ylim = c(0,.9))
title('n = 100 x 5', cex.main = 3)
incl <- pleioSimRes.nrep5$n == 500
plot.contBoxes(x = pleioSimRes.nrep5$rhog[incl], y = pleioSimRes.nrep5[incl, 6:7], legend.text = c('slope', 'slope_add'), legend.title = 'Models', ylab = 'Estimated h2 slope', xlab = 'Genetic cor', legend.cex = 3, spacer2 = 2, ylim = c(0,.9))
title('n = 500 x 5', cex.main = 3)
incl <- pleioSimRes.nrep5$n == 1000
plot.contBoxes(x = pleioSimRes.nrep5$rhog[incl], y = pleioSimRes.nrep5[incl, 6:7], legend.text = c('slope', 'slope_add'), legend.title = 'Models', ylab = 'Estimated h2 slope', xlab = 'Genetic cor', legend.cex = 3, spacer2 = 2, ylim = c(0,.9))
title('n = 1000 x 5', cex.main = 3)
dev.off()

#nRep = 10
pleioSimRes.nrep10 <- pleioSimRes[pleioSimRes$nrep == 10, ]
png('polycor/issues/08-varianceLeakage/figures/pleioSim_randSlope_10reps_h2slope.png', width = 1200, height = 800)
par(mfrow = c(1,3))
incl <- pleioSimRes.nrep10$n == 100
plot.contBoxes(x = pleioSimRes.nrep10$rhog[incl], y = pleioSimRes.nrep10[incl, 6:7], legend.text = c('slope', 'slope_add'), legend.title = 'Models', ylab = 'Estimated h2 slope', xlab = 'Genetic cor', legend.cex = 3, spacer2 = 2, ylim = c(0,.9))
title('n = 100 x 10', cex.main = 3)
incl <- pleioSimRes.nrep10$n == 500
plot.contBoxes(x = pleioSimRes.nrep10$rhog[incl], y = pleioSimRes.nrep10[incl, 6:7], legend.text = c('slope', 'slope_add'), legend.title = 'Models', ylab = 'Estimated h2 slope', xlab = 'Genetic cor', legend.cex = 3, spacer2 = 2, ylim = c(0,.9))
title('n = 500 x 10', cex.main = 3)
incl <- pleioSimRes.nrep10$n == 1000
plot.contBoxes(x = pleioSimRes.nrep10$rhog[incl], y = pleioSimRes.nrep10[incl, 6:7], legend.text = c('slope', 'slope_add'), legend.title = 'Models', ylab = 'Estimated h2 slope', xlab = 'Genetic cor', legend.cex = 3, spacer2 = 2, ylim = c(0,.9))
title('n = 1000 x 10', cex.main = 3)
dev.off()

#nRep = 20
pleioSimRes.nrep20 <- pleioSimRes[pleioSimRes$nrep == 20, ]
png('polycor/issues/08-varianceLeakage/figures/pleioSim_randSlope_20reps_h2slope.png', width = 1200, height = 800)
par(mfrow = c(1,3))
incl <- pleioSimRes.nrep20$n == 100
plot.contBoxes(x = pleioSimRes.nrep20$rhog[incl], y = pleioSimRes.nrep20[incl, 6:7], legend.text = c('slope', 'slope_add'), legend.title = 'Models', ylab = 'Estimated h2 slope', xlab = 'Genetic cor', legend.cex = 3, spacer2 = 2, ylim = c(0,.9))
title('n = 100 x 20', cex.main = 3)
incl <- pleioSimRes.nrep20$n == 500
plot.contBoxes(x = pleioSimRes.nrep20$rhog[incl], y = pleioSimRes.nrep20[incl, 6:7], legend.text = c('slope', 'slope_add'), legend.title = 'Models', ylab = 'Estimated h2 slope', xlab = 'Genetic cor', legend.cex = 3, spacer2 = 2, ylim = c(0,.9))
title('n = 500 x 20', cex.main = 3)
incl <- pleioSimRes.nrep20$n == 1000
plot.contBoxes(x = pleioSimRes.nrep20$rhog[incl], y = pleioSimRes.nrep20[incl, 6:7], legend.text = c('slope', 'slope_add'), legend.title = 'Models', ylab = 'Estimated h2 slope', xlab = 'Genetic cor', legend.cex = 3, spacer2 = 2, ylim = c(0,.9))
title('n = 1000 x 20', cex.main = 3)
dev.off()
