#' A light-weight version of lme4qtl::relmatLmer
#'
#' This function overcomes issues when fitting nested models, 
#' such as gene-env or multi-trait models.
#' In a nutshel, it is all about (i) updating Zt* = R*Z, where R'R = A;
#' and (ii) control under indexing or multiple instances of the same Zt matrix.
#'
#' Note that missing data should be handled before calling the function,
#' e.g. no missing values for two traits in a bi-variate model.
#' It is ok for doing simulations; todo when coming to real data analysis.
#'
#' @export
relmer <- function(formula, data = NULL, REML = TRUE,
  control = lmerControl(), start = NULL,
  verbose = 0L, subset, weights, na.action, offset,
  contrasts = NULL, 
  # relmat-specific arguments  
  relmat = list(), relfac = list(),
  check.nobs = "ignore", calc.derivs = TRUE,
  method.relfac = "auto",
  update = 2, debug = 0)
{
  stopifnot(require("Matrix"))
  stopifnot(require("lme4"))
  stopifnot(require("lme4qtl"))

  mc <- mcout <- match.call()

  if(check.nobs == "ignore") {
    control$checkControl$check.nobs.vs.rankZ <- "ignore"
    control$checkControl$check.nobs.vs.nlev <- "ignore"
    control$checkControl$check.nobs.vs.nRE <- "ignore"
  }
  control$calc.derivs <- calc.derivs 
  
  mc$control <- control
  
  ### Step 1: parse formula
  mc[[1]] <- quote(lme4::lFormula)

  args <- c("relmat", "relfac", "check.nobs", "calc.derivs", "method.relfac", "update", "debug")
  for(arg in args) {
    if(arg %in% names(mc)) {
      mc[[arg]] <- NULL
    }
  }

  lmod <- eval(mc, parent.frame(1L))
  
  ### Step 1.5 (relmat/relfac-specific)
  stopifnot(is.list(relmat), length(names(relmat)) == length(relmat))
  stopifnot(is.list(relfac), length(names(relfac)) == length(relfac))
  
  stopifnot(!any(names(relmat) %in% names(relfac)))
  stopifnot(!any(names(relfac) %in% names(relmat)))
  
  relnms <- c(names(relmat), names(relfac))
  relsrc <- c(rep("relmat", length(relmat)), rep("relfac", length(relfac)))
  
  flist <- lmod$reTrms[["flist"]]
  fnmns <- names(flist) 

  if(debug) {
    cat(" - fnmns:", fnmns, "\n")
    cat(" - relnms:", relnms, "\n")
  }
  # process `relmat/relfac` to update `Zt`
  if(any(fnmns %in% relnms)) {
    Ztlist <- lmod$reTrms[["Ztlist"]]
    
    ind <- which(fnmns %in% relnms)
    if(debug) {
      cat(" - ind:", ind, "\n")
    }
    for(i in ind) {
      fn <- fnmns[i] # "ID"
      if(debug) {
        cat("  -- fn:", fn, "\n")
      }
      obs <- rownames(lmod$fr) # just rownames like "1" "2" "3", ...
      zn <- lmod$fr[, fn] # ID values: 101 102 103 104
      zn.unique <- levels(zn)

      Zt <- Ztlist[[i]]
      rownames.Zt <- rownames(Zt) # rownames are (repeated) ID values: ID101 ID102 ID103...
      colnames.Zt <- colnames(Zt) # values are meaningless: "1" "2" "3", ..., but correspond to `zn`
      nrow.Zt <- nrow(Zt)
      ncol.Zt <- ncol(Zt)
      if(debug) {
        cat("  -- nrow.Zt:", nrow.Zt, "\n")
        cat("  -- ncol.Zt:", ncol.Zt, "\n")
      }
      
      # update `colnames.Zt`: match `obs` & `zn`
      adat <- data.frame(obs = obs, zn = zn) # annotation `dat`
      zdat <- data.frame(obs = colnames.Zt)
      mdat <- merge(zdat, adat, by = "obs", all.x = TRUE) # merged `dat`
      colnames.Zt <- mdat$zn
      
      # compute `R`: A = R' R 
      if(relsrc[i] == "relmat") {
        reln <- rownames(relmat[[fn]])
        stopifnot(!is.null(reln))
        stopifnot(all(zn.unique %in% reln))

        A <- relmat[[fn]][zn.unique, zn.unique]
        A <- Matrix::Matrix(A, sparse = TRUE)

        R <- relfac(A, method.relfac = method.relfac) # A = R' R 
        dim.R <- nrow(R) # nrow(R) = ncol(R)
      } else if(relsrc[i] == "relfac") {
        reln <- rownames(relfac[[fn]])
        stopifnot(!is.null(reln))
        stopifnot(all(zn.unique %in% reln))
        
        R0 <- relfac[[fn]][zn.unique, ] # nrow(R0) != ncol(R0)
        stopifnot(ncol(R0) <= nrow(R0))
        if(ncol(R0) < nrow(R0)) {
          Rpad <- Matrix(0, nrow = nrow(R0), ncol = nrow(R0) - ncol(R0), sparse = TRUE)
          R <- cbind(R0, Rpad)
        }
        dim.R <- nrow(R) # nrow(R) != ncol(R)
      } else {
        stop("relsrc[i]")
      }
          
      #----------
      # update v1
      #----------
      if(update == 1) {
      
      rows <- sapply(zn.unique, function(x) head(which(rownames.Zt == x), 1))
      cols <- sapply(zn.unique, function(x) head(which(colnames.Zt == x), 1))
      Zt.unique <- Zt[rows, cols]
      
      # new `Zt.unique` matrix, named as `Wt.unique`
      Wt.unique <- R %*% Zt.unique # Z* = Z L, then Z*' = L' Z' = R Z'
      
      rows <- sapply(rownames.Zt, function(x) which(zn.unique == x))
      cols <- sapply(colnames.Zt, function(x) which(zn.unique == x)) 
      Wt <- Wt.unique[rows, cols]
      
      } # end of update v1
      #-----------
      # update v2
      #-----------
      if(update == 2) {

      if(nrow.Zt == dim.R & ncol.Zt == dim.R) {
        Wt <- R %*% Zt
      } else {
        stopifnot(!(nrow.Zt %% dim.R))
        stopifnot(!(ncol.Zt %% dim.R))
        nrep.row <- nrow.Zt / dim.R
        nrep.col <- ncol.Zt / dim.R
  
        ind.rows <- lapply(zn.unique, function(x) which(rownames.Zt == x)) 
        ind.cols <- lapply(zn.unique, function(x) which(colnames.Zt == x)) 
        stopifnot(all(sapply(ind.rows, length) == nrep.row))
        stopifnot(all(sapply(ind.cols, length) == nrep.col))
  
        Wt <- Zt
        for(ir in seq(1, nrep.row)) {
          for(ic in seq(1, nrep.col)) {
            rows <- sapply(ind.rows, function(x) x[ir])
            cols <- sapply(ind.cols, function(x) x[ic])
            Wt[rows, cols] <- R %*% Wt[rows, cols]
          }
        }
      }
      
      } # end of update v2
      
      # store updated matrices 
      Ztlist[[i]] <- Wt
    }
    
    # update `lmod`
    lmod$reTrms[["Ztlist"]] <- Ztlist
    lmod$reTrms[["Zt"]] <- do.call(rBind, Ztlist)
  }
  
  ### Step 2: deviance function
  mcout$formula <- lmod$formula
  
  devfun <- do.call(mkLmerDevfun, 
    c(lmod, list(start = start, verbose = verbose, control = control)))

  ### Step 3: optimize
  opt <- optimizeLmer(devfun,
    optimizer = control$optimizer,
    restart_edge = control$restart_edge,
    boundary.tol = control$boundary.tol,
    control = control$optCtrl,
    verbose = verbose,
    start = start,
    calc.derivs = control$calc.derivs,
    use.last.params = control$use.last.params)
        
  cc <- lme4:::checkConv(attr(opt,"derivs"), opt$par,
    ctrl = control$checkConv,
    lbound = environment(devfun)$lower)  

  ### Step 4: make a model object
  mod <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr, mcout, lme4conv = cc)
  
  ### addition checks
  if(sd(mod@frame[, 1]) < sigma(mod)) {
    # the variance of response > residual variance
    warning("sd(y) > s2")
  }
  return(mod)
}

# Examples of use:
if(FALSE) {

library(tidyr)
library(lme4)
library(lme4qtl)

data(dat40)
# no missing data is a current requirement of `relmer` for nested model
dat <- na.omit(dat40)

# sex-specificity model of gene-env. interaction for a single trait
m1 <- lme4qtl::relmatLmer(trait1 ~ AGE + SEX + (0 + SEX | ID), dat, relmat = list(ID = kin2))
m2 <- relmer(trait1 ~ AGE + SEX + (0 + SEX | ID), dat, relmat = list(ID = kin2))
stopifnot(all.equal(getME(m1, "Ztlist")[[1]]@x, getME(m2, "Ztlist")[[1]]@x))
stopifnot(all.equal(logLikNum(m1), logLikNum(m2)))

# genetic correlation model for two traits
bdat <- within(dat, RID <- ID)
bdat <- tidyr::gather(bdat, tname, trait, trait1, trait2)

m3 <- relmer(trait ~ (0 + tname|ID) + (0 + dummy(tname)|RID), bdat, relmat = list(ID = kin2))
#> VarCorr(m3)
# Groups   Name         Std.Dev. Corr
# ID       tnametrait1  2.32297
#          tnametrait2  2.55500  0.913
# RID      dummy(tname) 0.56241
# Residual              0.75802
 
wts <- rep(1e10, nrow(bdat))
m4 <- relmer(trait ~ (0 + tname|ID) + (0 + tname|RID), bdat, weights = wts, relmat = list(ID = kin2), calc.derivs = FALSE)
#> VarCorr(m4)
# Groups   Name        Std.Dev. Corr
# ID       tnametrait1 2.20402
#          tnametrait2 2.45187  0.900
# RID      tnametrait1 0.92409
#          tnametrait2 1.07239  0.299
# Residual             1.23058
 
# SOLAR models
library(solarius)
solarPolygenic(trait1 ~ 1, dat)$vcf[1, 2] # h2 = 0.8652342
solarPolygenic(trait2 ~ 1, dat)$vcf[1, 2] # h2 = 0.8561782
solarPolygenic(trait1 + trait2 ~ 1, dat)$vcf # rhog = 0.9124713; rhoe = 0.2963357
}
