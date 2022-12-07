source('polycor/simFunctions.R')
source('polycor/relmer.R')
#Here I'm just making sure that the pleiotropic simulation does not induce any correlation within individuals

#### Simulation ####
n_ind <- 50
n_rep <- 20
simdat <- pleiotropy.simulatePairKin(n_ind = n_ind, n_perInd = n_rep, h2 = .5, overlap = .8)
cor.test(simdat$pheno$trait1, simdat$pheno$trait2) #Population level correlation
#Within individual Correlations
indCors <- data.frame(cor = rep(NA, n_ind), p = NA)
for (i in 1:n_ind) {
  dat <- simdat$pheno[simdat$pheno$id == i, ]
  dat.cor <- cor.test(dat$trait1, dat$trait2)
  indCors$cor[i] <- dat.cor$estimate
  indCors$p[i] <- dat.cor$p.value
}


#### Plot ####
source('polycor/plotFunctions.R')
png(paste('polycor/issues/08-varianceLeakage/figures/pleioSimExample_nInd', n_ind, '_nRep', n_rep, '.png', sep = ''), width = 1200, height = 800)
layout(rbind(1:2, 3:4))
#Pop
plot(simdat$pheno$trait1, simdat$pheno$trait2, xlab = '', ylab = '')
abline(lm(simdat$pheno$trait2 ~ simdat$pheno$trait1), lty = 2, col = 'red', lwd = 4)
title('Wole Popultion', cex.main = 2)
mtext('Trait 2', side = 2, line = 2, cex = 2)
mtext('Trait 1', side = 1, line = 3, cex = 2)
#Ind example
ind <- simdat$pheno$id == 20
plot(simdat$pheno$trait1[ind], simdat$pheno$trait2[ind], xlab = '', ylab = '')
abline(lm(simdat$pheno$trait2[ind] ~ simdat$pheno$trait1[ind]), lty = 2, col = 'red', lwd = 4)
title('Individual 20', cex.main = 2)
mtext('Trait 2', side = 2, line = 2, cex = 2)
mtext('Trait 1', side = 1, line = 3, cex = 2)
#All within ind cors
boxplot(indCors$cor)
mtext('Individual Correlation', side = 2, line = 2, cex = 2)
abline(h = 0, lty = 2)
#All within ind p-vals
boxplot(indCors$p)
mtext('p-val Individual Correlation', side = 2, line = 2, cex = 2)
abline(h = .05, lty = 2)
dev.off()
