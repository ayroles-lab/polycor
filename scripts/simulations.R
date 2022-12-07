##### Simulate a pair of traits #####
source('simFunctions.R')
n_ind <- 200
n_perInd <- 3
H2 <- .5
sim.pair <- genCov.simulatePair(n_ind = n_ind, n_perInd = n_perInd, H2 = H2)

#Calculate kinship
library(Matrix)
test <- apply(sim.pair$geno, 2, var) != 0
G <- scale(sim.pair$geno[,test]) %*% t(scale(sim.pair$geno[,test]))
G <- G/mean(diag(G))
rownames(G) <- colnames(G) <- 1:n_ind #lme4qtl requires rownames for some reason
if(!all(eigen(G)$values > 0)){
  G.posdef <- nearPD(G) #Negative eigenvalues causes lme4qtl to crash. Solution from https://github.com/variani/lme4qtl/issues/1
  G <- G.posdef$mat
}

#Fit models
library(lme4)
library(lme4qtl)
sim.pair_randSlopeModel_kinship <- relmatLmer(gene2.expr ~ gene1.expr + (0 + gene1.expr|id), sim.pair$pheno, relmat = list(id = G))
#h2 estimate, or whatever it should be called in this case
vf <- as.data.frame(VarCorr(sim.pair_randSlopeModel_kinship))[, c("grp", "vcov")]
vf$vcov[1]/sum(vf$vcov) 
#Same model with uncorrelated random effects
sim.pair_randSlopeModel_unCor <- lmer(gene2.expr ~ gene1.expr + (0 + gene1.expr|id), sim.pair$pheno)
#h2 estimate, or whatever it should be called in this case
vf <- as.data.frame(VarCorr(sim.pair_randSlopeModel_unCor))[, c("grp", "vcov")]
vf$vcov[1]/sum(vf$vcov) #About right

#blups VS fixed eff estimates
#According to Andrey https://github.com/variani/lme4qtl/issues/13
R <- chol(G)
blups.kinship <- sim.pair_randSlopeModel_kinship@beta[2] + R %*% lme4::ranef(sim.pair_randSlopeModel_kinship)$id[,1]
blups.uncor <- sim.pair_randSlopeModel_unCor@beta[2] + lme4::ranef(sim.pair_randSlopeModel_unCor)$id[,1]
plot(blups.kinship, blups.uncor)

#Repeated lms
fixedSlopes <- numeric(n_ind)
for(j in 1:n_ind){
  obs <- sim.pair$pheno$id == j
  model <- summary(lm(sim.pair$pheno$gene2.expr[obs] ~ sim.pair$pheno$gene1.expr[obs]))
  fixedSlopes[j] <- model$coefficients[2,1]
}

#Plot
tmp <- c(fixedSlopes, blups.kinship[,1], blups.uncor)
ylim <- c(min(tmp), max(tmp))
plot(unique(sim.pair$pheno$g_cov), fixedSlopes, pch = 19, xlab = '', ylab = '', ylim = ylim)
mtext(text = 'Estimated Slope', side = 2, cex = 2, line = 2.5)
mtext(text = 'True Slope', side = 1, cex = 2, line = 3)
points(unique(sim.pair$pheno$g_cov), blups.kinship, pch = 19, col = 'red', ylim = ylim)
points(unique(sim.pair$pheno$g_cov), blups.uncor, pch = 19, col = 'green', ylim = ylim)
legend('topleft', c('BLUP', 'BLUP uncor', 'Fixed Effect'), pch = 19, col = c('red', 'green', 'black'))


##### Repeated simulations #####
n <- 20 #Number of simulations
n_ind <- 100
n_perInd <- 3
H2 <- .5
h2s.est <- data.frame(intraCor.kin = rep(NA, n), intraCor.unCor = NA)
estimateVStrue <- data.frame(blup_kinship.cor = rep(NA, n), blup_unCor.cor = NA, fxed.cor = NA)
for(i in 1:n){
  simulation <- genCov.simulatePair(n_ind = n_ind, n_perInd = n_perInd, H2 = H2)
  
  #Calculate G matrix, excluding sites with no variation in genotypes
  test <- apply(simulation$geno, 2, var) != 0
  G <- scale(simulation$geno[,test]) %*% t(scale(simulation$geno[,test]))
  G <- G/mean(diag(G))
  rownames(G) <- colnames(G) <- 1:n_ind #lme4qtl requires rownames for some reason
  G.posdef <- nearPD(G) #Negative eigenvalues causes lme4qtl to crash. Solution from https://github.com/variani/lme4qtl/issues/1
  
  #Fit random slope model with kinship
  randInt_kin.model <- relmatLmer(gene2.expr ~ gene1.expr + (0 + gene1.expr|id), simulation$pheno, relmat = list(id = G.posdef$mat))
  #h2 estimate, or whatever it should be called in this case
  vf <- as.data.frame(VarCorr(randInt_kin.model))[, c("grp", "vcov")]
  h2s.est$intraCor.kin[i] <- vf$vcov[1]/sum(vf$vcov) #About right
  R <- chol(G.posdef$mat)
  blups.kin <- randInt_kin.model@beta[2] + R %*% lme4::ranef(randInt_kin.model)$id[,1]
  
  #Fit random slope model without kinship
  randInt_uncor.model <- lmer(gene2.expr ~ gene1.expr + (0 + gene1.expr|id), simulation$pheno)
  #h2 estimate, or whatever it should be called in this case
  vf <- as.data.frame(VarCorr(randInt_uncor.model))[, c("grp", "vcov")]
  h2s.est$intraCor.unCor[i] <- vf$vcov[1]/sum(vf$vcov) #About right
  blups.uncor <- randInt_uncor.model@beta[2] + lme4::ranef(randInt_uncor.model)$id[,1]
  
  #Repeated lms
  fixedSlopes <- numeric(n_ind)
  for(j in 1:n_ind){
    obs <- simulation$pheno$id == j
    model <- summary(lm(simulation$pheno$gene2.expr[obs] ~ simulation$pheno$gene1.expr[obs]))
    fixedSlopes[j] <- model$coefficients[2,1]
  }
  estimateVStrue$blup_kinship.cor[i] <- cor(unique(simulation$pheno$g_cov), blups.kin[,1])
  estimateVStrue$blup_unCor.cor[i] <- cor(unique(simulation$pheno$g_cov), blups.uncor)
  estimateVStrue$fxed.cor[i] <- cor(unique(simulation$pheno$g_cov), fixedSlopes)
}
boxplot(h2s.est, xaxt = 'n')
axis(side = 1, at = 1:2, labels = c('Yes', 'No'), cex.axis = 2)
mtext(text = 'Kinship in model', side = 1, cex = 2, line = 3, font = 2)
mtext(text = 'Estimated \"h2\" of the slope', side = 2, cex = 2, line = 2.5)

plot(h2s.est)

boxplot(estimateVStrue$blup_unCor.cor, estimateVStrue$blup_kinship.cor, estimateVStrue$fxed.cor, xaxt = 'n')
axis(side = 1, at = 1:3, labels = c('BLUP_uncor', 'BLUP_kinship', 'Fixed Effect'), cex.axis = 1.7)
mtext(text = 'Cor(true, estimate)', side = 2, cex = 2, line = 2.5)


##### Simulate "standard" additive genetic architecture #####
#Code used to generate the figure on slide 10
library(hglm)
library(lme4qtl)
n_loci <- 5000
n_ind <- 500
n <- 10
h2s.est_myG <- list()
h2s.true <- seq(.1, .9, .1)

ptm <- Sys.time()
for(i in 1:length(h2s.true)){
  h2s.tmp <- data.frame(h2.lmeqtl = rep(NA, n), h2.hglm = NA)
  
  for(j in 1:n){
    #Draw allele frequencies from a U-shaped distribution
    tmp <- seq(from = .01, to = .99, by = .01)
    prob <- 1/(tmp*(1-tmp))
    prob <- prob/max(prob)
    p <- sample(size = n_loci, x = tmp , prob = prob, replace = T)
    
    #Assign genotypes to according to HW
    geno <- matrix(nrow = n_ind, ncol = n_loci)
    for(k in 1:n_ind){
      geno[k,] <- sapply(X = p, FUN = function(x){sample(x = 0:2, size = 1, prob = c(x^2, 2*x*(1-x), (1-x)^2))})
    }
    
    #Assign phenotypes
    beta <- rnorm(n_loci)
    g <- geno %*% beta
    y <- g + rnorm(n = length(g), mean = 0, sd = sd(g)*((1-h2s.true[i])/h2s.true[i]))
    
    #Kinship
    test <- apply(geno, 2, var) != 0
    G <- scale(geno[,test]) %*% t(scale(geno[,test]))
    G <- G/mean(diag(G))
    rownames(G) <- colnames(G) <- 1:n_ind #lme4qtl requires rownames for some reason
    if(!all(eigen(G)$values > 0)){
      G.posdef <- nearPD(G) #Negative eigenvalues causes lme4qtl to crash. Solution from https://github.com/variani/lme4qtl/issues/1
      G <- G.posdef$mat
    }

    #Fit models
    #hglm
    svd <- svd(G)
    Z <- svd$u %*% diag(sqrt(svd$d))
    y.hglm <- hglm(X = matrix(1, nrow = n_ind, ncol = 1), y = y, Z = Z)
    h2s.tmp$h2.hglm[j] <- y.hglm$varRanef/(y.hglm$varRanef + y.hglm$varFix)
    #lme4qtl
    pheno <- data.frame(id = 1:n_ind, pheno = y)
    y.lme4qtl <- relmatLmer(pheno ~ (1|id), data = pheno, relmat = list(id = G))
    vf <- as.data.frame(VarCorr(y.lme4qtl))[, c("grp", "vcov")]
    h2s.tmp$h2.lmeqtl[j] <- vf$vcov[1]/sum(vf$vcov) #About right
    
    print(paste('h2', h2s.true[i], 'number', j))
  }
  h2s.est_myG[[i]] <- h2s.tmp
}
Sys.time() - ptm #~1.1h

#reshape
tmp <- sapply(h2s.est_myG, function(x){x$h2.lmeqtl})
tmp2 <- sapply(h2s.est_myG, function(x){x$h2.hglm})
h2s.est_myG.frame <- data.frame(h2.true = rep(h2s.true, each = 10), h2.lmeqtl = c(tmp), h2.hglm = c(tmp2))

#plot
library(wesanderson)
cols <- wes_palette('Darjeeling2', 2)
spacer1 <- .01
spacer2 <- 3
xPos <- c(h2s.est_myG.frame$h2.true - spacer1, h2s.est_myG.frame$h2.true + spacer1)
box <- boxplot(c(h2s.est_myG.frame$h2.lmeqtl, h2s.est_myG.frame$h2.hglm) ~ xPos, xaxt = 'n', yaxt = 'n', col = cols)
tmp <- seq(1, length(box$names)*spacer2/2, spacer2)
at = sort(c(tmp, tmp+1))
png('../results/2018-12-05_compare_hglm_lmeqtl_etc/normalAddSim_comp_hglm_lme4qtl.png', width = 1200, height = 800)
boxplot(c(h2s.est_myG.frame$h2.lmeqtl, h2s.est_myG.frame$h2.hglm) ~ xPos, xaxt = 'n', at = at, yaxt = 'n', col = cols)
grid(10,10)
par(new = T)
boxplot(c(h2s.est_myG.frame$h2.lmeqtl, h2s.est_myG.frame$h2.hglm) ~ xPos, xaxt = 'n', at = at, yaxt = 'n', col = cols)
axis(side = 1, at = seq(1.5, length(box$names)*spacer2/2 - .5, spacer2), labels = h2s.true, las = 2, cex.axis = 2)
axis(side = 2, at = seq(0,1,.1), cex.axis = 2)
mtext(text = 'Simulated h2', side = 1, line = 4, cex = 2, font = 2)
mtext(text = 'Estimated h2', side = 2, line = 2.5, cex = 2, font = 2)
legend('topleft', c('lme4qtl', 'hglm'), col = cols, title = 'h2 estimates', pch = 15, cex = 2.5)
dev.off()



