library(devtools)
load_all("~/git/variani/lme4qtl")
#library(lme4qtl)

library(Matrix)

M <- 1000 # n_loci
N <- 500 # n_ind
h2 <- 0.5

### reproducibility
set.seed(1)

### simulate genotypes
freqs <- rep(0.5, M) # fix all freq to 0.5 

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
#svd <- svd(G)
#W <- svd$u %*% diag(sqrt(svd$d))

# evd
evd <- eigen(G, symmetric = TRUE)
#R <- diag(sqrt(out$values)) %*% t(out$vectors)
R <- relfac.evd(G, tol = 1e-10) # lme4 requires that class(R) = "dgeMatrix"

G2 <- nearPD(G)$mat
  
### fit models
m1 <- relmatLmer(y ~ (1|id), data = dat, relmat = list(id = G))
