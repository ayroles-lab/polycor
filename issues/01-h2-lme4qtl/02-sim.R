library(hglm)
library(lme4qtl)
library(lme4)
library(Matrix)

M <- 1000 # n_loci
N <- 500 # n_ind
n <- 10 #Number of simulations per true h2
h2s.true <- seq(.1, .9, .1)
h2s.est <- list()

ptm <- Sys.time()
for(i in 1:length(h2s.true)){
  h2 <- h2s.true[i]
  h2s.tmp <- data.frame(h2.lmeqtl = rep(NA, n), h2.hglm = NA)
  
  for(j in 1:n){
    ### simulate genotypes
    freqs <- runif(M, 0.1, 0.9)
    
    # genotypes of independent individuals
    X <- sapply(1:M, function(i) {
      p <- freqs[i]
      q <- 1 - p
      rbinom(N, 2, p)
    })
    
    ### compute GRM
    # The GRM composed of these M (scaled) genotypes computed as:
    # GRM_jk = (1/M) Sum_i^M ((x_ij - 2p_i)(x_ik - 2p_i) / 2p_i(1-p_i)
    
    # center/scale 
    col_means <- colMeans(X, na.rm = TRUE)
    col_freq <- col_means / 2  # col_means = 2 * col_freq
    col_sd <- sqrt(2 * col_freq * (1 - col_freq))
    Z <- sweep(X, 2, col_means, "-")
    Z <- sweep(Z, 2, col_sd , "/")
    Zg <- Z / sqrt(M)
    
    # GRM 
    G <- tcrossprod(Zg) # the same as tcrossprod(Z) / M
    
    ### simulate phenotype
    b <- rnorm(M, 0, sqrt(h2/M)) 
    y <- Z %*% b + rnorm(N, 0, sqrt(1 - h2))
    
    ### prepare data for model fitting
    ids <- as.character(1:N)
    dat <- data.frame(id = ids, y = y)
    rownames(G) <- ids
    colnames(G) <- ids
    
    # svd
    svd <- svd(G)
    W <- svd$u %*% diag(sqrt(svd$d))
    
    ### fit models
    #lme4qtl
    npd <- nearPD(G)
    Gnpd <- npd$mat
    rownames(Gnpd) <- ids
    colnames(Gnpd) <- ids
    m2 <- relmatLmer(y ~ (1|id), data = dat, relmat = list(id = Gnpd))
    vf <- as.data.frame(VarCorr(m2))[, c("grp", "vcov")]
    h2s.tmp$h2.lmeqtl[j] <- vf$vcov[1]/sum(vf$vcov) #About right
    
    #hglm
    X0 <- matrix(1, nrow = N, ncol = 1)
    m3 <- hglm(X = X0, y = y, Z = W)
    h2s.tmp$h2.hglm[j] <- m3$varRanef/(m3$varRanef + m3$varFix)
    
    print(paste('h2', h2s.true[i], 'number', j))
  }
  h2s.est[[i]] <- h2s.tmp
}
Sys.time() - ptm #~42m
save(h2s.est, file = '02-sim.RData')

#reshape
tmp <- sapply(h2s.est, function(x){x$h2.lmeqtl})
tmp2 <- sapply(h2s.est, function(x){x$h2.hglm})
h2s.est.frame <- data.frame(h2.true = rep(h2s.true, each = 10), h2.lmeqtl = c(tmp), h2.hglm = c(tmp2))
#plot
library(wesanderson)
cols <- wes_palette('Darjeeling2', 2)
spacer1 <- .01
spacer2 <- 3
xPos <- c(h2s.est.frame$h2.true - spacer1, h2s.est.frame$h2.true + spacer1)
box <- boxplot(c(h2s.est.frame$h2.lmeqtl, h2s.est.frame$h2.hglm) ~ xPos, xaxt = 'n', yaxt = 'n', col = cols)
tmp <- seq(1, length(box$names)*spacer2/2, spacer2)
at = sort(c(tmp, tmp+1))
png('normalAddSim_comp_hglm_lme4qtl.png', width = 1200, height = 800)
boxplot(c(h2s.est.frame$h2.lmeqtl, h2s.est.frame$h2.hglm) ~ xPos, xaxt = 'n', at = at, yaxt = 'n', col = cols)
grid(10,10)
par(new = T)
boxplot(c(h2s.est.frame$h2.lmeqtl, h2s.est.frame$h2.hglm) ~ xPos, xaxt = 'n', at = at, yaxt = 'n', col = cols)
axis(side = 1, at = seq(1.5, length(box$names)*spacer2/2 - .5, spacer2), labels = h2s.true, las = 2, cex.axis = 2)
axis(side = 2, at = seq(0,1,.1), cex.axis = 2)
mtext(text = 'Simulated h2', side = 1, line = 4, cex = 2, font = 2)
mtext(text = 'Estimated h2', side = 2, line = 2.5, cex = 2, font = 2)
legend('topleft', c('lme4qtl', 'hglm'), col = cols, title = 'h2 estimates', pch = 15, cex = 2.5)
dev.off()
