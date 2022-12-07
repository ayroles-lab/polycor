genCov.simulatePair <- function(n_loci = 200, n_ind = 50, n_perInd = 3, 
                                h2, trait1 = NULL, geno = NULL, freqDistr = 'U', returnGRM = T, nrPops = NULL, Fst = NULL){
  #In this simulation, expression of gene 1 is only given by E and cor(expr gene 1, expr gene 2) is a polygenic trait
  #I simulate three scenarios:
  # 1. gene2 ~ gene1*g.cov
  # 2. gene2 ~ gene1*g.cov + g.add
  # 3. gene2 ~ g.add
  #Where:
  # g.cov = polygenic effect on cor(gene1, gene2)
  # g.add = polygenic effect on gene2
  
  #The argument h2 isn't heritability strictly speaking (I think), but it does control the signal/noise ratio of the genetic effects
  if(h2 < 0 | h2 > 1)
    stop('h2 must be in the range [0,1]')
  
  
  #### Simulate Genotypes ####
  if(is.null(geno)){
    if(!is.null(nrPops) & !is.null(Fst)){ #Simulate stratified population
      geno <- sim_pop(N = n_ind, M = n_loci, Fst = Fst, nrPops = nrPops)
      n_ind <- n_ind*nrPops
    }
    else{ #Simulate unrelated individuals
      #Allele frequency distribution
      if(freqDistr == 'U' | freqDistr == 'u'){
        #Draw alleles from from a U-shaped freq distribution
        tmp <- seq(from = .01, to = .99, by = .01)
        prob <- 1/(tmp*(1-tmp))
        prob <- prob/max(prob)
        p <- sample(size = n_loci, x = tmp , prob = prob, replace = T)
      }
      else if(freqDistr == 'uni' | freqDistr == 'uniform')
        p <- runif(n_loci, 0.1, 0.9)
      else if(freqDistr > 0 & freqDistr < 1)
        p <- rep(freqDistr, n_loci)
      else
        stop(paste('Did not recognize freqDistr:', freqDistr))
      
      #Assign genotypes according to HW
      #This step could be elaborated to increase relatedness between certain individuals
      geno <- matrix(nrow = n_ind, ncol = n_loci)
      for(i in 1:n_ind){
        #X[i,] <- sapply(X = p, FUN = function(x){sample(x = -1:1, size = 1, prob = c(x^2, 2*x*(1-x), (1-x)^2))})
        haplo1 <- as.numeric(runif(n_loci) < p)
        haplo2 <- as.numeric(runif(n_loci) < p)
        geno[i,] <- haplo1 + haplo2
      }
    }
    
    #Remove non-polymorphic loci
    fixed <- apply(geno, 2, var) == 0
    if(any(fixed)){
      geno <- geno[, !fixed]
      n_loci <- sum(!fixed)
    }
  }
  
  # center/scale
  col_means <- colMeans(geno, na.rm = TRUE)
  col_freq <- col_means / 2  # col_means = 2 * col_freq
  col_sd <- sqrt(2 * col_freq * (1 - col_freq))
  Z <- sweep(geno, 2, col_means, "-")
  Z <- sweep(Z, 2, col_sd , "/")
  
  if(returnGRM){
    Zg <- Z / sqrt(n_loci)
    G <- tcrossprod(Zg) # the same as tcrossprod(Z) / M
    rownames(G) <- colnames(G) <- 1:n_ind 
  }
  else
    G <- NA
  
  #### Simulate phenotypes ####
  mu1 <- 4
  mu2 <- 0
  data.repObs <- data.frame(id = sort(rep(c(1:n_ind), n_perInd)), 
                            obs = rep(1:n_perInd, n_ind), 
                            gene1.expr = NA, 
                            gene2.expr = NA,
                            gene2.expr_withAdd = NA,
                            gene2.expr_onlyAdd = NA,
                            g.add = NA,
                            g.cov = NA)
  
  #Genetic effect per locus
  b.add <- rnorm(n_loci, 0, sqrt(h2/n_loci)) 
  b.cov <- rnorm(n_loci, 0, sqrt(h2/n_loci)) 
  #Adds up to genetic effect per individual. Saving in order to compare to blup estimates
  g.add <- Z %*% b.add
  g.cov <- Z %*% b.cov
  
  data.repObs$g.add <- g.add[data.repObs$id]
  data.repObs$g.cov <- g.cov[data.repObs$id]
  #Expression
  if(is.null(trait1))
    data.repObs$gene1.expr <- mu1 + rnorm(n_ind*n_perInd) #Draw a random trait1
  else
    data.repObs$gene1.expr <- trait1 #Use the provided trait1
  
  # data.repObs$gene2.expr <- mu2 + data.repObs$gene1.expr*data.repObs$g.cov + rnorm(n = n_ind*n_perInd, mean = 0, sd = sd(g.cov)*((1-h2)/h2))
  # data.repObs$gene2.expr_withAdd <- mu2 + data.repObs$gene1.expr*data.repObs$g.cov + data.repObs$g.add + rnorm(n = n_ind*n_perInd, mean = 0, sd = (sd(g.add) + sd(g.cov))*((1-h2)/h2))
  # data.repObs$gene2.expr_onlyAdd <- mu2 + data.repObs$g.add + rnorm(n = n_ind*n_perInd, mean = 0, sd = sd(g.add)*((1-h2)/h2))
  
  data.repObs$gene2.expr <- mu2 + data.repObs$gene1.expr*data.repObs$g.cov + rnorm(n = n_ind*n_perInd, mean = 0, sd = sqrt(1 - h2))
  data.repObs$gene2.expr_withAdd <- mu2 + data.repObs$gene1.expr*data.repObs$g.cov + data.repObs$g.add + rnorm(n = n_ind*n_perInd, mean = 0, sd = sqrt(1 - h2))
  data.repObs$gene2.expr_onlyAdd <- mu2 + data.repObs$g.add + rnorm(n = n_ind*n_perInd, mean = 0, sd = sqrt(1 - h2))
  
  return(list(pheno = data.repObs, geno.raw = geno, geno.scaled = Z, G = G))
}


pleiotropy.simulatePair <- function(n_loci = 200, n_ind = 50, n_perInd = 3, h2, trait1 = NULL, geno = NULL, freqDistr = 'U', returnGRM = T, nrPops = NULL, Fst = NULL, overlap){
  if(h2 < 0 | h2 > 1)
    stop('h2 must be in the range [0,1]')
  if(overlap < 0 | overlap > 1)
    stop('overlap must be in the range [0,1]')
  
  #### Simulate Genotypes ####
  if(is.null(geno)){
    if(!is.null(nrPops) & !is.null(Fst)){ #Simulate stratified population
      geno <- sim_pop(N = n_ind, M = n_loci, Fst = Fst, nrPops = nrPops)
      n_ind <- n_ind*nrPops
    }
    else{ #Simulate unrelated individuals
      #Allele frequency distribution
      if(freqDistr == 'U' | freqDistr == 'u'){
        #Draw alleles from from a U-shaped freq distribution
        tmp <- seq(from = .01, to = .99, by = .01)
        prob <- 1/(tmp*(1-tmp))
        prob <- prob/max(prob)
        p <- sample(size = n_loci, x = tmp , prob = prob, replace = T)
      }
      else if(freqDistr == 'uni' | freqDistr == 'uniform')
        p <- runif(n_loci, 0.1, 0.9)
      else if(freqDistr > 0 & freqDistr < 1)
        p <- rep(freqDistr, n_loci)
      else
        stop(paste('Did not recognize freqDistr:', freqDistr))
      
      #Assign genotypes according to HW
      #This step could be elaborated to increase relatedness between certain individuals
      geno <- matrix(nrow = n_ind, ncol = n_loci)
      for(i in 1:n_ind){
        #X[i,] <- sapply(X = p, FUN = function(x){sample(x = -1:1, size = 1, prob = c(x^2, 2*x*(1-x), (1-x)^2))})
        haplo1 <- as.numeric(runif(n_loci) < p)
        haplo2 <- as.numeric(runif(n_loci) < p)
        geno[i,] <- haplo1 + haplo2
      }
    }
    
    #Remove non-polymorphic loci
    fixed <- apply(geno, 2, var) == 0
    if(any(fixed)){
      geno <- geno[, !fixed]
      n_loci <- sum(!fixed)
    }
  }
  
  # center/scale
  col_means <- colMeans(geno, na.rm = TRUE)
  col_freq <- col_means / 2  # col_means = 2 * col_freq
  col_sd <- sqrt(2 * col_freq * (1 - col_freq))
  Z <- sweep(geno, 2, col_means, "-")
  Z <- sweep(Z, 2, col_sd , "/")
  
  if(returnGRM){
    Zg <- Z / sqrt(n_loci)
    G <- tcrossprod(Zg) # the same as tcrossprod(Z) / M
    rownames(G) <- colnames(G) <- 1:n_ind 
  }
  else
    G <- NA
  
  #### Simulate phenotypes ####
  mu1 <- 1
  mu2 <- 1
  data.repObs <- data.frame(id = sort(rep(c(1:n_ind), n_perInd)), 
                            obs = rep(1:n_perInd, n_ind), 
                            trait1 = NA, 
                            trait2 = NA,
                            g.trait1 = NA,
                            g.trait2 = NA)
  
  #Genetic effect per locus. Amount of overlap in causative loci determined by the overlap argument
  b.trait1 <- b.trait2 <- rnorm(n_loci, 0, sqrt(h2/n_loci))
  mid <- floor(median(1:n_loci))
  b.trait1[1:(mid - round(overlap*n_loci/2))] <- 0
  b.trait2[(mid + round(overlap*n_loci/2)):n_loci] <- 0
  
  #Adds up to genetic effect per individual. Saving in order to compare to blup estimates
  g.trait1 <- Z %*% b.trait1
  g.trait2 <- Z %*% b.trait2
  
  data.repObs$g.trait1 <- g.trait1[data.repObs$id]
  data.repObs$g.trait2 <- g.trait2[data.repObs$id]
  
  data.repObs$trait1 <- mu1 + data.repObs$g.trait1 + rnorm(n = n_ind*n_perInd, mean = 0, sd = sqrt(1 - h2))
  data.repObs$trait2 <- mu1 + data.repObs$g.trait2 + rnorm(n = n_ind*n_perInd, mean = 0, sd = sqrt(1 - h2))
  return(list(pheno = data.repObs, geno.raw = geno, geno.scaled = Z, G = G))
}

#--------------------------------------------------
# Simulation of two traits (related ind.; kinship)
#--------------------------------------------------

simulateKinship <- function(N, L = 5)
{
  if(!require(Matrix))
    stop('Could not load package Matrix')
  
  stopifnot(!(N %% L)) # family size = L
  Nl <- N / L
  
  # kinship of a single family
  K0 <- matrix(0.5, L, L)
  diag(K0) <- 1
  K0[1, 2] <- 0
  K0[2, 1] <- 0
  
  K <- kronecker(Diagonal(Nl), K0)
  
  ids <- seq(N)
  rownames(K) <- ids
  colnames(K) <- ids
  
  return(K)
}

# Partially copied from: https://github.com/hemostat/polycor/blob/master/issues/03-bivar-lmer/runs/02-sim-rand-noise/R/01-run-gen-rep.R
pleiotropy.simulatePairKin <- function(
  n_ind = 250, n_perInd = 4, # 250 x 4 samples seems ok to recover rhog 
    # see also: https://github.com/hemostat/polycor/issues/3#issuecomment-453800742
  h2 = 0.75, # the same h2 for trait1 and trait2
  overlap = 0.5 # genetic correlation (rhog)
    # here we assume there is no residual correlation (rhoe)
)
{
  K <- simulateKinship(n_ind)
  
  #### Simulate phenotypes ####
  mu1 <- 1
  mu2 <- 2
  data.repObs <- data.frame(id = rep(c(1:n_ind), each = n_perInd), 
    rep = rep(1:n_perInd, n_ind), obs = seq(1, n_ind*n_perInd),
    trait1 = NA, trait2 = NA,
    g.trait1 = NA, g.trait2 = NA)
  
  # Genetic part of two genetically correlated traits
  h12 <- overlap*h2
  g.sigma <- matrix(c(h2, h12, h12, h2), 2, 2)
  v.sigma <- kronecker(g.sigma, K)

  g.trait12 <- mvrnorm(1, rep(0, 2*n_ind), v.sigma)

  #Adds up to genetic effect per individual. Saving in order to compare to blup estimates
  g.trait1 <- g.trait12[seq(1, n_ind)]
  g.trait2 <- g.trait12[n_ind + seq(1, n_ind)]
  
  data.repObs$g.trait1 <- g.trait1[data.repObs$id]
  data.repObs$g.trait2 <- g.trait2[data.repObs$id]

  data.repObs$r.trait1 <- rnorm(n = n_ind*n_perInd, mean = 0, sd = sqrt(1 - h2))
  data.repObs$r.trait2 <- rnorm(n = n_ind*n_perInd, mean = 0, sd = sqrt(1 - h2))
    
  data.repObs$trait1 <- mu1 + data.repObs$g.trait1 + data.repObs$r.trait1
  data.repObs$trait2 <- mu2 + data.repObs$g.trait2 + data.repObs$r.trait2
  
  return(list(pheno = data.repObs, G = K))
}


genCov.simulatePairKin <- function(n_ind = 250, n_perInd = 4, 
  h2.cov = 0.4, h2.add = 0.4,
  mu1 = 0, mu2 = 2)
{
  # In this simulation: 
  #  - expression of trait1 is only given by E and 
  #  - cor(trait1, trait2) is a polygenic trait
  # Three scenarios:
  # 1. trait2 ~ trait1 * g.cov
  # 2. trait2 ~ trait1 * g.cov + g.add
  # 3. trait2 ~ g.add
  # where:
  # g.cov = polygenic effect on cor(trait1, trait2)
  # g.add = polygenic effect on trait2
  
  stopifnot(mu1 == 0) # not necessary if trait1 is scaled before model fitting
  stopifnot(!(h2.cov + h2.add > 1))
  
  K <- simulateKinship(n_ind)
  
  #### Simulate phenotypes ####
  dat <- data.frame(id = rep(c(1:n_ind), each = n_perInd), 
    rep = rep(1:n_perInd, n_ind), obs = seq(1, n_ind*n_perInd),
    trait1 = NA, 
    trait2 = NA, trait2_withAdd = NA, trait2_onlyAdd = NA,
    g.add = NA, g.cov = NA)
  
  # Genetic effect per locus
  if(h2.add) {
    g.add <- mvrnorm(1, rep(0, n_ind), h2.add * K)
  } else {
    g.add <- rep(0, n_ind)
  }
  if(h2.cov) {
    g.cov <- mvrnorm(1, rep(0, n_ind), h2.cov * K)
  } else {
    g.cov <- rep(0, n_ind)
  }
  
  dat$g.add <- g.add[dat$id]
  dat$g.cov <- g.cov[dat$id]
  
  # Trait 1
  dat$trait1 <- mu1 + rnorm(n_ind*n_perInd) # Draw a random trait1
  
  # Trait 2 under three scenarios
  dat$trait2 <- mu2 + dat$g.cov * dat$trait1 + 
    rnorm(n_ind*n_perInd, mean = 0, sd = sqrt(1 - h2.cov))
    
  dat$trait2_withAdd <- mu2 + dat$g.cov * dat$trait1 + dat$g.add + 
    rnorm(n_ind*n_perInd, mean = 0, sd = sqrt(1 - h2.cov - h2.add))

  dat$trait2_onlyAdd <- mu2 + dat$g.add + 
    rnorm(n_ind*n_perInd, mean = 0, sd = sqrt(1 - h2.add))
  
  return(list(pheno = dat, G = K))
}

genBivarCov.simulatePairKin <- function(n_ind = 250, n_perInd = 4, 
  h2.cov = 0.4, h2.add = 0.4,
  h2.add.t1 = 0.4,
  rhog = 0,
  mu1 = 0, mu2 = 2,
  postproc = FALSE)
{
  # In this simulation: 
  #  - expression of trait1 is only given by G and E
  #  - trait1 and trait2 can also come from a bi-variate model
  #  - cor(trait1, trait2) is a polygenic trait
  # Simulation:
  #   trait1 ~ g.add.t1 + g.bivar + e.t1
  #   trait2 ~ g.cov * trait1 + g.add.t2 + g.bivar + e.t2
  # where:
  # g.cov = polygenic effect on cor(trait1, trait2)
  # g.add = additive polygenic effect on trait1/trait2
  # g.bivar = pleiotropuc effect on both trait1 and trait2
  #
  # var(g.cov) = h2.cov * Kinship
  # var(g.add.t2) = h2.add * Kinship
  # var(g.add.t1) = h2.add.t1 * Kinship
  # var(g.bivar) = h2.bivar * Kinship = rhog * sqrt(h2.add) * sqrt(h2.add.t1) * Kinship
  
  stopifnot(mu1 == 0) # not necessary if trait1 is scaled before model fitting
  
  # need non-zero variances; at least variances close to zero, but positive
  stopifnot(h2.cov > 0)
  stopifnot(h2.add > 0)
  stopifnot(h2.add.t1 > 0)
  
  h2.bivar <- rhog * sqrt(h2.add) * sqrt(h2.add.t1)
  stopifnot(!(h2.cov + h2.add + h2.bivar > 1))
  
  K <- simulateKinship(n_ind)
  
  #### Simulate phenotypes ####
  dat <- data.frame(id = rep(c(1:n_ind), each = n_perInd), 
    rep = rep(1:n_perInd, n_ind), obs = seq(1, n_ind*n_perInd),
    trait1 = 0.0, trait2 = 0.0, 
    g.add = 0.0, g.cov = 0.0,
    g.add.t1 = 0.0, g.bivar = 0.0)
  
  # Genetic effect per locus
  if(rhog) {
    g.sigma <- matrix(c(h2.add.t1, h2.bivar, h2.bivar, h2.add), 2, 2)
    v.sigma <- kronecker(g.sigma, K)

    g.add.t12 <- mvrnorm(1, rep(0, 2*n_ind), v.sigma)
  
    g.add.t1 <- g.add.t12[seq(1, n_ind)]
    g.add <- g.add.t12[n_ind + seq(1, n_ind)]
  } else {
    g.add.t1 <- mvrnorm(1, rep(0, n_ind), h2.add.t1 * K)
    g.add <- mvrnorm(1, rep(0, n_ind), h2.add * K)
  }
  
  g.cov <- mvrnorm(1, rep(0, n_ind), h2.cov * K)
  
  dat$g.add.t1 <- g.add.t1[dat$id]
  dat$g.add <- g.add[dat$id]
  dat$g.cov <- g.cov[dat$id]
  
  dat$trait1 <- mu1 + dat$g.add.t1 + dat$g.bivar +
    rnorm(n_ind*n_perInd, sd = sqrt(1 - h2.add.t1 - h2.bivar))
  
  dat$trait2 <- mu2 + dat$g.cov * dat$trait1 + 
    dat$g.add + dat$g.bivar +
    rnorm(n_ind*n_perInd, sd = sqrt(1 - h2.add - h2.bivar - h2.cov))

  if(postproc) {
    dat <- within(dat, {
      id <- as.character(id)
      rid <- id
      trait1 <- as.numeric(scale(trait1))
      trait2 <- as.numeric(scale(trait2))
    })  
  }

  return(list(pheno = dat, G = K))
}

genBivarCov.simulatePair_outbred <- function(n_loci = 200, n_ind = 250, n_perInd = 4, 
                                        h2.add.t1 = 0.4, h2.add.t2 = 0.4, # The heritabilities of the two traits. If rhog > 0, h2.add.t2 will approach h2.add.t1
                                        h2.cov = 0.4,                     # The heritability of the intra-individual correlation/causality t1 -> t2
                                        rhog = 0,                         # Genetic correlation (pleiotropy). Here, it is defined as the fraction of overlapping causal loci. Hence rhog is proportional, but not equal to, genetic correlation
                                        mu1 = 0, mu2 = 0,                 # The intercepts/population means for the traits
                                        beta_12 = 0,                      # The population mean of the effect t1 -> t2
                                        geno = NULL,                      # If null, the genotype is simulated
                                        freqDistr = 'U',                  # Distribution to draw allele frequencies from 
                                        returnGRM = T, 
                                        returnGeno = F,
                                        nrPops = NULL, Fst = NULL,        # Controls the degree of genetic stratification in the simulated population
                                        postproc = FALSE,
                                        e.cov = 0,
                                        trait1.binary = F, 
                                        trait1.fixed = F)                 #trait1 the same across all individuals. For instance equal time points when trait2 is measured. This should affect var(trait2)
{
  #Same scenario as in the function genBivarCov.simulatePairKin() but simulating genotypes in outbred populations, rather than nuclear families
  #The GRM is SNP based and hence more dense than in the *Kin functions
  #NOTE: rhog is the fraction of overlapping causal loci, which is proportional but not equal to the genetic correlation. When rhog > 0, h2.add.t2 and h2.cov should also be proportional but not equal to the true parameter values
  # Simulation:
  #   trait1 ~ g.add.t1 + e.t1
  #   trait2 ~ g.add.t2 + g.cov * trait1 + e.t2
  # Where:
  #   g.cov = polygenic effect on cor(trait1, trait2)
  #   g.add = additive polygenic effect on trait1/trait2. g.add per individual is the sum of all the SNP effects, where the fraction of overlapping causal SNPs is determined by rhog

  stopifnot(h2.add.t1 >= 0 & h2.add.t1 <= 1)
  totVar.t2 <- h2.cov + h2.add.t2 + e.cov
  stopifnot(totVar.t2 >= 0 & totVar.t2 <= 1)

  #### Simulate Genotypes ####
  if(is.null(geno)){
    if(!is.null(nrPops) & !is.null(Fst)){ #Simulate stratified population
      stopifnot(!n_ind %% nrPops)
      geno <- sim_pop(N = n_ind/nrPops, M = n_loci, Fst = Fst, nrPops = nrPops)
    }
    else{ #Simulate unrelated individuals
      #Allele frequency distribution
      if(freqDistr == 'U' | freqDistr == 'u'){
        #Draw alleles from from a U-shaped freq distribution
        tmp <- seq(from = .01, to = .99, by = .01)
        prob <- 1/(tmp*(1-tmp))
        prob <- prob/max(prob)
        p <- sample(size = n_loci, x = tmp , prob = prob, replace = T)
      }
      else if(freqDistr == 'uni' | freqDistr == 'uniform')
        p <- runif(n_loci, 0.1, 0.9)
      else if(freqDistr > 0 & freqDistr < 1)
        p <- rep(freqDistr, n_loci)
      else
        stop(paste('Did not recognize freqDistr:', freqDistr))
      
      #Assign genotypes according to HW
      geno <- matrix(nrow = n_ind, ncol = n_loci)
      for(i in 1:n_ind){
        # geno[i,] <- sapply(X = p, FUN = function(x){sample(x = -1:1, size = 1, prob = c(x^2, 2*x*(1-x), (1-x)^2))}) #Same thing as below but slower
        haplo1 <- as.numeric(runif(n_loci) < p)
        haplo2 <- as.numeric(runif(n_loci) < p)
        geno[i,] <- haplo1 + haplo2
      }
    }
    
    #Remove non-polymorphic loci
    fixed <- apply(geno, 2, var) == 0
    if(any(fixed)){
      geno <- geno[, !fixed]
      n_loci <- sum(!fixed)
    }
  }
  
  # Center & scale
  col_means <- colMeans(geno, na.rm = TRUE)
  col_sd <- apply(geno, 2, sd, na.rm = TRUE)
  # col_freq <- col_means / 2  # col_means = 2 * col_freq
  # col_sd <- sqrt(2 * col_freq * (1 - col_freq))
  Z <- sweep(geno, 2, col_means, "-")
  Z <- sweep(Z, 2, col_sd , "/")
  
  if(returnGRM){
    Zg <- Z / sqrt(n_loci)
    G <- tcrossprod(Zg) # the same as tcrossprod(Z) / M
    rownames(G) <- colnames(G) <- 1:n_ind 
  }
  else
    G <- NA
  
  
  #### Simulate phenotypes ####
  dat <- data.frame(id = rep(c(1:n_ind), each = n_perInd),
                    rep = rep(1:n_perInd, n_ind), obs = seq(1, n_ind*n_perInd),
                    trait1 = 0.0, trait2 = 0.0,
                    g.add.t1 = 0.0, g.add.t2 = 0.0,
                    g.cov = 0.0)
  
  #Additive genetic effect per locus. Amount of overlap in causative loci determined by rhog. This is what induces the pleiotropy / genetic correlation
  b.trait1 <- rnorm(n_loci, 0, sqrt(h2.add.t1/n_loci))
  b.trait2 <- rnorm(n_loci, 0, sqrt(h2.add.t2/n_loci))
  if(rhog > 0)
    b.trait2[1:round(rhog*n_loci)] <- b.trait1[1:round(rhog*n_loci)]
  
  #Adds up to genetic effect per individual. Saving in order to compare to blup estimates
  g.add.t1 <- Z %*% b.trait1
  g.add.t2 <- Z %*% b.trait2
  dat$g.add.t1 <- g.add.t1[dat$id]
  dat$g.add.t2 <- g.add.t2[dat$id]
  
  #Genetic effect per locus on individual correlation
  b.cov <- rnorm(n_loci, 0, sqrt(h2.cov/n_loci)) 
  #Adds up to genetic effect per individual. Saving in order to compare to blup estimates
  g.cov <- Z %*% b.cov
  dat$g.cov <- g.cov[dat$id]
  
  #Non genetic effect on individual slope
  dat$res.cov <- rep(rnorm(n_ind, sd = sqrt(e.cov)), each = n_perInd)
  
  #Individual phenotypes are the sum of the genetic effects + environmental noise
  if(trait1.fixed){
    t <- 1:n_perInd - mean(1:n_perInd)
    dat$trait1 <- rep(t, times = n_ind)
  }
  else{
    y1 <- mu1 + dat$g.add.t1 + rnorm(n = n_ind*n_perInd, mean = 0, sd = sqrt(1 - h2.add.t1))
    if(trait1.binary)
      dat$trait1 <- as.numeric(cut(y1, breaks = 2)) - 1 #Discretize trait1
    else
      dat$trait1 <- y1
  }
  
  dat$trait2 <- mu2 + dat$g.add.t2 + (beta_12 + dat$g.cov + dat$res.cov) * dat$trait1 + rnorm(n = n_ind*n_perInd, mean = 0, sd = sqrt(1 - h2.add.t2 - h2.cov - e.cov))
  
  if(postproc) {
    dat <- within(dat, {
      id <- as.character(id)
      rid <- id
      #trait1 <- as.numeric(scale(trait1))
      #trait2 <- as.numeric(scale(trait2))
    })  
  }
  
  if(!returnGeno)
    geno <- NULL
  
  return(list(pheno = dat, G = G, geno = geno))
}

#--------------------------------------
# Simulation of multiple traits (>2)
#--------------------------------------
#Note: These functions are still quite experimental

genCov.simulateMultiTrait <- function(n_loci = 200, n_ind = 50, n_perInd = 5, n_traits = 10, h2){
  if(!require(Matrix))
    stop('Couldn\'t load package Matrix')
  if(!require(MASS))
    stop('Couldn\'t load package MASS')
  
  #Structure to store the data
  data.repObs <- data.frame(id = sort(rep(c(1:n_ind), n_perInd)), 
                            obs = rep(1:n_perInd, n_ind),
                            matrix(nrow = n_ind*n_perInd, ncol = n_traits))
  colnames(data.repObs)[3:ncol(data.repObs)] <- paste('trait', 1:n_traits, sep = '')
  
  #Draw alleles from from a U-shaped freq distribution
  tmp <- seq(from = .01, to = .99, by = .01)
  prob <- 1/(tmp*(1-tmp))
  prob <- prob/max(prob)
  p <- sample(size = n_loci, x = tmp , prob = prob, replace = T)
  
  #Assign genotypes according to HW
  #This step could be elaborated to increase relatedness between certain individuals
  geno <- matrix(nrow = n_ind, ncol = n_loci)
  for(i in 1:n_ind){
    #X[i,] <- sapply(X = p, FUN = function(x){sample(x = -1:1, size = 1, prob = c(x^2, 2*x*(1-x), (1-x)^2))})
    haplo1 <- as.numeric(runif(n_loci) < p)
    haplo2 <- as.numeric(runif(n_loci) < p)
    geno[i,] <- haplo1 + haplo2
  }
  
  #Calculate G matrix, excluding sites with no variation in genotypes
  test <- apply(geno, 2, var) != 0
  G <- scale(geno[,test]) %*% t(scale(geno[,test]))
  G <- G/mean(diag(G))
  rownames(G) <- colnames(G) <- 1:n_ind #lme4qtl requires rownames for some reason
  if(!all(eigen(G)$values > 0)){
    G.posdef <- nearPD(G) #Negative eigenvalues causes lme4qtl to crash. Solution from https://github.com/variani/lme4qtl/issues/1
    G <- G.posdef$mat
  }
  
  #Create the base trait covariance matrix. Is it sensible to just draw random cov like this?
  A <- matrix(runif(n_traits^2)*2-1, ncol=n_traits) 
  cov.base <- t(A) %*% A
  
  #Create individual trait cov-matrices
  #Genetic effect on the individual cov-matrix. The effects are correlated according to kinship
  #Here, I'm drawing one term per trait pair and individual. Another way might be to have just one term per individual that is added to the cov-matrix
  #g is a matrix with one col per individual. There are as many rows as there are trait pairs. 
  g <- mvrnorm(n = choose(n_traits, 2), mu = rep(0, n_ind), Sigma = G)
  
  g.colSDs <- apply(g, 2, sd)
  j <- 1
  indCovMat <- list()
  for(i in 1:n_ind){
    #Create individual cov-matrix by adding noise and genetic effect to cov.base
    cov.ind <- cov.base
    e <- rnorm(n = nrow(g), mean = 0, sd = g.colSDs[i]*((1-h2)/h2))
    cov.ind[upper.tri(cov.ind)] <- cov.ind[upper.tri(cov.ind)] + g[,i]
    cov.ind[lower.tri(cov.ind)] <- t(cov.ind)[lower.tri(cov.ind)] #Make symetric
    #Check if it is positive definite
    if(!all(eigen(cov.ind)$values > 0)){ #I hope this wont create any problems
      print(paste('Covariance matrix for individual', i, 'is not PD. Fixing it with nearPD'))
      cov.ind_pd <- nearPD(cov.ind)
      data.repObs[j:(j+n_perInd-1), 3:ncol(data.repObs)] <- mvrnorm(n = n_perInd, mu = rep(0, n_traits), Sigma = cov.ind_pd$mat) + e
    }else
      data.repObs[j:(j+n_perInd-1), 3:ncol(data.repObs)] <- mvrnorm(n = n_perInd, mu = rep(0, n_traits), Sigma = cov.ind) + e
    
    indCovMat[[i]] <- cov.ind
    j <- j + n_perInd
  }
  
  return(list(data = data.repObs, indCovMat = indCovMat, geno = geno, G = G))
}

build_Rmatrix <- function(simRes, reshapePheno = T){
  #Convenience function to build the block structured cov-matrix. Input is the output from the genCov.simulateMultiTrait() function
  
  #If reshapePheno = T, a vector is also returned with the phenotypes reshaped into one long vector, and a vector with ids per ind*trait
  #The vectors and the R matrix should all add up correctly to be used in lme4qtl. All the indices are a mindfuck but I think I got it right ...
  n_ind <- length(unique(simRes$data$id))
  n_traits <- ncol(simRes$data) - 2
  n_perInd <- length(unique(simRes$data$obs))
  if(reshapePheno){
    y <- numeric(n_ind*n_perInd*n_traits)
    id <- character(n_ind*n_traits)
  }
  else{
    y <- NA
    id <- NA
  }
  
  #There should be some elegant way of doing this using matrix multiplications
  R <- matrix(0, nrow = n_ind*n_traits, ncol = n_ind*n_traits)
  j <- 1
  m <- 1
  for(i in 1:n_ind){
    data <- simRes$data[simRes$data$id == i, 3:ncol(simRes$data)]
    cov.ind_i <- var(data)
    
    if(reshapePheno){
      y[m:(m + n_perInd*n_traits - 1)] <- unlist(data)
      id[m:(m + n_perInd*n_traits - 1)] <- rep(paste('ind', i, '_trait', 1:n_traits, sep = ''), each = n_perInd)
      m <- m + n_perInd*n_traits
    }
    
    R[j:(j+n_traits-1), j:(j+n_traits-1)] <- cov.ind_i #The individual cov-matrix on the diagonal
    if(i < n_ind){
      #Fill the other blocks on this "row"
      l <- j+n_traits
      for(k in (i+1):n_ind){
        data <- simRes$data[simRes$data$id == k, 3:ncol(simRes$data)]
        cov.ind_k <- var(data)
        
        R[j:(j+n_traits-1), l:(l+n_traits-1)] <- cov.ind_i*cov.ind_k*simRes$G[i,k]
        l <- l+n_traits
      }
    }
    j <- j+n_traits
  }
  R[lower.tri(R)] <- t(R)[lower.tri(R)] #Fill in the lower triangle
  
  return(list(R = R, y = y, id = id))
}

genCov.simulateMultiTrait_v2 <- function(n_loci = 200, n_ind = 50, n_perInd = 5, n_traits = 10, h2, netCaseProbs = NULL, returnG = T){
  #A completely different approach. Here, I'm repeatedly calling genCov.simulatePair to simulate trait pairs with a genetic effect on the covariance
  #I think this is a more explicit way of building the correlation network, although computationaly much less efficient
  if(!require(Matrix))
    stop('Couldn\'t load package Matrix')
  
  #Structure to store the data
  data.repObs <- data.frame(id = sort(rep(c(1:n_ind), n_perInd)), 
                            obs = rep(1:n_perInd, n_ind),
                            matrix(nrow = n_ind*n_perInd, ncol = n_traits))
  colnames(data.repObs)[3:ncol(data.repObs)] <- paste('trait', 1:n_traits, sep = '')
  geneticEffects <- data.frame(id = sort(rep(c(1:n_ind), n_perInd)), 
                               obs = rep(1:n_perInd, n_ind),
                               matrix(nrow = n_ind*n_perInd, ncol = n_traits-1))
  
  #Draw alleles from from a U-shaped freq distribution
  tmp <- seq(from = .01, to = .99, by = .01)
  prob <- 1/(tmp*(1-tmp))
  prob <- prob/max(prob)
  p <- sample(size = n_loci, x = tmp , prob = prob, replace = T)
  
  #Assign genotypes according to HW
  #This step could be elaborated to increase relatedness between certain individuals
  geno <- matrix(nrow = n_ind, ncol = n_loci)
  for(i in 1:n_ind){
    #X[i,] <- sapply(X = p, FUN = function(x){sample(x = -1:1, size = 1, prob = c(x^2, 2*x*(1-x), (1-x)^2))})
    haplo1 <- as.numeric(runif(n_loci) < p)
    haplo2 <- as.numeric(runif(n_loci) < p)
    geno[i,] <- haplo1 + haplo2
  }
  
  if(returnG){
    test <- apply(geno, 2, var) != 0
    G <- scale(geno[,test]) %*% t(scale(geno[,test]))
    G <- G/mean(diag(G))
    rownames(G) <- colnames(G) <- 1:n_ind #lme4qtl requires rownames for some reason
    if(!all(eigen(G)$values > 0)){
      G.posdef <- nearPD(G) #Negative eigenvalues causes lme4qtl to crash. Solution from https://github.com/variani/lme4qtl/issues/1
      G <- G.posdef$mat
    }
  }
  else
    G <- NA
  
  y1 <- data.repObs$trait1 <- rnorm(n_ind*n_perInd)
  traitA.name <- 1
  traitNr <- 2
  for(i in 4:(n_traits + 2)){
    #Every pair of traits is generated by  
    #trait2 ~ trait1*g.cov , where g.cov is a polygenic effect
    #Ad hoc way of building the network: 
    #case1: trait1 = trait2 in the previous round. Creates a causal chain A -> B -> C ...
    #case2: trait1 = trait1 in the previous round. One trait directly affecting many A -> B, A -> C ...
    #case3: trait1 ~ E, randomly draw a new trait1. Creates separate "modules" A -> B, C -> D ...
    case <- sample(1:3, size = 1, prob = netCaseProbs)

    if(case == 1 & exists('simPair')){
      simPair <- genCov.simulatePair(n_loci = n_loci, n_ind = n_ind, n_perInd = n_perInd, h2 = h2, trait1 = data.repObs[, i-1], geno = geno)
      traitA.name <- traitNr-1
    }
    else if(case == 2){
      simPair <- genCov.simulatePair(n_loci = n_loci, n_ind = n_ind, n_perInd = n_perInd, h2 = h2, trait1 = y1, geno = geno)
    }
    else{
      y1 <- rnorm(n_ind*n_perInd)
      traitA.name <- 'e'
      simPair <- genCov.simulatePair(n_loci = n_loci, n_ind = n_ind, n_perInd = n_perInd, h2 = h2, trait1 = y1, geno = geno)
    }
    data.repObs[, i] <- simPair$pheno$gene2.expr
    
    geneticEffects[, i-1] <- simPair$pheno$g.cov
    names(geneticEffects)[i-1] <- paste(traitA.name, '->', traitNr, sep = '')
    traitNr <- traitNr+1
  }
  
  return(list(data = data.repObs, geneticEffects = geneticEffects, geno = geno, G = G))
}


sim_pop <- function(N = 200, M = 1000, Fst = 0.1, maf_max = 0.5, maf_min = 0.05, nrPops = 2){ 
  #Adopted from https://variani.github.io/bigcov/vignettes/popstrat.html
  # set.seed(seed)
  maf_values <- runif(M, maf_min, maf_max)
  
  gdat <- matrix(ncol = M, nrow = N*nrPops)
  k <- 1
  for(j in 1:nrPops){
    freq.pop <- sapply(1:M, function(i) rbeta(1, 
                                              maf_values[i] * (1 - Fst) / Fst, 
                                              (1 - maf_values[i]) * (1 - Fst) / Fst))
    
    gdat.pop <- sapply(1:M, function(i) sample(c(0, 1, 2), N, replace = TRUE,
                                               prob = c(((1 - freq.pop[i])^2), (2 * freq.pop[i] * (1 - freq.pop[i])), (freq.pop[i]^2))))
    gdat[k:(k+N-1), ] <- gdat.pop
    k <- k + N
  }
  return(gdat)
}

