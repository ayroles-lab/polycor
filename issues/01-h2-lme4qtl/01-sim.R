library(hglm)
library(lme4qtl)

library(Matrix)

M <- 1000 # n_loci
N <- 500 # n_ind
h2 <- 0.9

### reproducibility
set.seed(1)

### simulate genotypes
freqs <- rep(0.5, M) # fix all freq to 0.5 
#freqs <- runif(M, 0.1, 0.9)

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

# evd
#out <- eigen(G, symmetric = TRUE)
#R <- diag(sqrt(out$values)) %*% t(out$vectors)
#R <- Matrix::Matrix(R)
  
### fit models
m1 <- relmatLmer(y ~ (1|id), data = dat, relmat = list(id = G))

#> VarProp(m1)
#       grp        var1 var2       vcov     sdcor      prop
#1       id (Intercept) <NA> 0.87055745 0.9330367 0.8980196
#2 Residual        <NA> <NA> 0.09886177 0.3144229 0.1019804

# observation 1: `m1` recovers the true value of h2 = 0.9

npd <- nearPD(G)
Gnpd <- npd$mat
rownames(Gnpd) <- ids
colnames(Gnpd) <- ids
m2 <- relmatLmer(y ~ (1|id), data = dat, relmat = list(id = Gnpd))

#> VarProp(m2)
#       grp        var1 var2       vcov     sdcor      prop
#1       id (Intercept) <NA> 0.87055743 0.9330367 0.8980196
#2 Residual        <NA> <NA> 0.09886179 0.3144229 0.1019804

# observation 2: `m2` performs as well as `m1`
# do we need nearPD? can we just copy upper-triangle part to lower-triangle part

X0 <- matrix(1, nrow = N, ncol = 1)
m3 <- hglm(X = X0, y = y, Z = W)
  
h2_m3 <- with(m3, varRanef / (varRanef + varFix))
#> h2_m3
#[1] 0.8924215

# observation 3: `m3` performs as well as `m1`
