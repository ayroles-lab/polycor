#This wrapper script simulates random slope data, and polygenic additive data, for a range of h2s and multiple repeats per h2.
#The actual simulation is done in the function genCov.simulatePair (simFunctions.R)
library(lme4)
library(lme4qtl)
library(magrittr)

#Change to your paths 
source('lme4qtl/R/utils.lib.R')  #To get the function relmatRanef
source('polycor/simFunctions.R') #To get the simulation function genCov.simulatePair

#### Simulation parameters ####
n_loci <- 1000    #Number of loci. All will be independent (no LD) and causal
n_ind <- 50       #Number of individuals
n_perInd <- 2     #Observations per individual
n <- 2            #Number of simulations per h2
freqDistr <- 'U'  #Shape of allele frequency distribution. Options: 'U' (U-shaped) or 'uni' (uniform)
Fst <- .3         #Controls the amount of stratification in the population (together with nrPops)
nrPops <- 10      #Number of sub-populations

#Data structures to store results
h2s.est <- list()
genEffects <- list()

h2s.true <- seq(.1, .9, .1) #The true h2s to simulate
for(i in 1:length(h2s.true)){
  h2 <- h2s.true[i]
  
  #Data structures to store results
  h2s.tmp <- data.frame(h2_slope.lmeqtl = rep(NA, n), h2_slope.lmer = NA, h2_add.lmeqtl = NA, h2_add.lmer = NA,
                        h2_slope_addData.lmeqtl = rep(NA, n), h2_slope_addData.lmer = NA, h2_add_polyCorData.lmeqtl = NA, h2_add_polyCorData.lmer = NA)
  genEffects.tmp <- list()
  
  for(j in 1:n){
    #### Simulate Data ####
    #The function genCov.simulatePair returns a list with four elements: 
    #   pheno       - Simulated phenotypes generated under three different scenarios (see the function for details). Genetic effects 
    #   geno.raw    - Unscaled genotypes (coded as 0,1,2)
    #   geno.scaled - Scaled genotypes 
    #   G           - GRM
    simPair <- genCov.simulatePair(n_loci = n_loci, n_ind = n_ind, n_perInd = n_perInd, h2 = h2, freqDistr = freqDistr, nrPops = nrPops, Fst = Fst)
    
    #Data structure to store results
    genEffects.tmp2 <- data.frame(g.add = rep(NA, n_ind*nrPops), g.cov = NA, g.add_blupLme4qtl = NA, g.add_blupLmer = NA, 
                                  g.cov_blupLme4qtl = NA, g.cov_blupLmer = NA, g.add_fixed = NA, g.cov_fixed = NA,
                                  g.add_polyCorData_blupLme4qtl = NA, g.add_polyCorData_blupLmer = NA,
                                  g.cov_addData_blupLme4qtl = NA, g.cov_addData_blupLmer = NA)
    
    
    #### Fit models ####
    #Phenotype: Polygenic correlation (random slope)
    dat <- data.frame(id = simPair$pheno$id, y = simPair$pheno$gene2.expr, x = simPair$pheno$gene1.expr)
    #Model:     Random slope (model specification matching data)
    m_slope.lme4qtl <- relmatLmer(y ~ x + (0 + x|id), data = dat, relmat = list(id = simPair$G))
    m_slope.lmer <- lmer(y ~ x + (0 + x|id), data = dat)
    #Model:     Random intercept (wrong model specification)
    m_add_polyCorData.lme4qtl <- relmatLmer(y ~ (1|id), data = dat, relmat = list(id = simPair$G))
    m_add_polyCorData.lmer <- lmer(y ~ (1|id), data = dat)
    
    #Phenotype: Polygenic additive architecture (NO random slope) 
    dat <- data.frame(id = simPair$pheno$id, y = simPair$pheno$gene2.expr_onlyAdd, x = simPair$pheno$gene1.expr)
    #Model:     Random intercept (model specification matching data)
    m_add.lme4qtl <- relmatLmer(y ~ (1|id), data = dat, relmat = list(id = simPair$G))
    m_add.lmer <- lmer(y ~ (1|id), data = dat)
    #Model:     Random slope (wrong model specification)
    m_slope_addData.lme4qtl <- relmatLmer(y ~ x + (0 + x|id), data = dat, relmat = list(id = simPair$G))
    m_slope_addData.lmer <- lmer(y ~ x + (0 + x|id), data = dat)
    
    #Repeated lms
    fixedSlopes <- numeric(n_ind*nrPops)
    fixedInt <- numeric(n_ind*nrPops)
    for(k in 1:length(fixedInt)){
      obs <- simPair$pheno$id == k
      model <- summary(lm(simPair$pheno$gene2.expr[obs] ~ simPair$pheno$gene1.expr[obs]))
      fixedSlopes[k] <- model$coefficients[2,1]
      fixedInt[k] <- mean(simPair$pheno$gene2.expr_onlyAdd[obs])
    }
    
    
    #### h2 estimates ####
    #Model specification matching data
    vf <- as.data.frame(VarCorr(m_slope.lme4qtl))[, c("grp", "vcov")]
    h2s.tmp$h2_slope.lmeqtl[j] <- vf$vcov[1]/sum(vf$vcov)
    vf <- as.data.frame(VarCorr(m_slope.lmer))[, c("grp", "vcov")]
    h2s.tmp$h2_slope.lmer[j] <- vf$vcov[1]/sum(vf$vcov)
    vf <- as.data.frame(VarCorr(m_add.lme4qtl))[, c("grp", "vcov")]
    h2s.tmp$h2_add.lmeqtl[j] <- vf$vcov[1]/sum(vf$vcov)
    vf <- as.data.frame(VarCorr(m_add.lmer))[, c("grp", "vcov")]
    h2s.tmp$h2_add.lmer[j] <- vf$vcov[1]/sum(vf$vcov)
    #Wrong model specification 
    vf <- as.data.frame(VarCorr(m_add_polyCorData.lme4qtl))[, c("grp", "vcov")]
    h2s.tmp$h2_add_polyCorData.lmeqtll[j] <- vf$vcov[1]/sum(vf$vcov)
    vf <- as.data.frame(VarCorr(m_add_polyCorData.lmer))[, c("grp", "vcov")]
    h2s.tmp$h2_add_polyCorData.lmer[j] <- vf$vcov[1]/sum(vf$vcov)
    vf <- as.data.frame(VarCorr(m_slope_addData.lme4qtl))[, c("grp", "vcov")]
    h2s.tmp$h2_slope_addData.lmeqtl[j] <- vf$vcov[1]/sum(vf$vcov)
    vf <- as.data.frame(VarCorr(m_slope_addData.lmer))[, c("grp", "vcov")]
    h2s.tmp$h2_slope_addData.lmer[j] <- vf$vcov[1]/sum(vf$vcov)
    
    
    #### Estimates (BLUPs and fixed eff.) & true genetic effects ####
    #BLUPs
    #Model specification matching data
    genEffects.tmp2$g.add_blupLmer <- m_add.lmer@beta + lme4::ranef(m_add.lmer)$id[,1]
    genEffects.tmp2$g.add_blupLme4qtl <- m_add.lme4qtl@beta + relmatRanef(m_add.lme4qtl, 'id')
    genEffects.tmp2$g.cov_blupLmer <- m_slope.lmer@beta[2] + lme4::ranef(m_slope.lmer)$id[,1]
    genEffects.tmp2$g.cov_blupLme4qtl <- m_slope.lme4qtl@beta[2] + relmatRanef(m_slope.lme4qtl, 'id')
    #Wrong model specification
    genEffects.tmp2$g.add_polyCorData_blupLme4qtl <- m_add_polyCorData.lme4qtl@beta + relmatRanef(m_add_polyCorData.lme4qtl, 'id')
    genEffects.tmp2$g.add_polyCorData_blupLmer <- m_add_polyCorData.lmer@beta + lme4::ranef(m_add_polyCorData.lmer)$id[,1]
    genEffects.tmp2$g.cov_addData_blupLme4qtl <- m_slope_addData.lme4qtl@beta + relmatRanef(m_slope_addData.lme4qtl, 'id')
    genEffects.tmp2$g.cov_addData_blupLmer <- m_slope_addData.lmer@beta + lme4::ranef(m_slope_addData.lmer)$id[,1]
    
    #Fixed effect estimates
    genEffects.tmp2$g.add_fixed <- fixedInt
    genEffects.tmp2$g.cov_fixed <- fixedSlopes
    
    #True effects
    genEffects.tmp2$g.add <- unique(simPair$pheno$g.add)
    genEffects.tmp2$g.cov <- unique(simPair$pheno$g.cov)
    
    #Put in storage
    genEffects.tmp[[j]] <- genEffects.tmp2
    
    print(paste('h2', h2s.true[i], 'number', j))
  }
  
  #Put in storage
  h2s.est[[i]] <- h2s.tmp
  genEffects[[i]] <- genEffects.tmp
}
save(h2s.est, genEffects, file = paste('181211_simPair_nLoci', n_loci, '_nInd', n_ind*nrPops, '_nPerInd', n_perInd, '_freqDistr', freqDistr, '_nSubPops', nrPops, '_Fst', sub(pattern = '\\.', replacement = '', x = as.character(Fst)), '_nSims', n, '.RData', sep = ''))
