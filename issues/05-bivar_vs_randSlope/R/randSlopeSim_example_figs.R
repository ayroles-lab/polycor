#Plotting function
plot.contBoxes <- function(x, y, col = NULL, spacer1 = .01, spacer2 = 3, xlab = 'Simulated h2', ylab = 'Estimated h2', legend.text = NULL, legend.title = NULL, legend.pos = 'topleft', legend.cex = 1, ...){
  if(is.null(col)){
    library(wesanderson)
    cols <- wes_palette('Darjeeling2', ncol(y))
  }
  
  #Set up grouping variable
  tmp <- seq(from = 0, length.out = ncol(y), by = spacer1)
  tmp <- tmp - mean(tmp)
  grouping <- rep(x, 2) + rep(tmp, each = length(x))
  #Set up between "group" spacing
  box <- boxplot(unlist(y) ~ grouping, xaxt = 'n', yaxt = 'n', col = cols, plot = F)
  tmp <- rep(seq(0, length.out = length(unique(x)), by = 2), each = ncol(y))
  xPos <- rep(1:ncol(y), length(unique(x))) + tmp*spacer2
  #pos for x-axis
  start <- mean(xPos[1:ncol(y)])
  end <- mean(xPos[(length(xPos) - ncol(y) + 1):length(xPos)])
  xlab.pos <- seq(start, end, length.out = length(unique(x)))
  
  #plot
  boxplot(unlist(y) ~ grouping, xaxt = 'n', at = xPos, yaxt = 'n', col = cols, ...)
  grid(10,10)
  par(new = T)
  boxplot(unlist(y) ~ grouping, xaxt = 'n', at = xPos, yaxt = 'n', col = cols, ...)
  axis(side = 1, at = xlab.pos, labels = unique(x), las = 2, cex.axis = 2)
  axis(side = 2, at = unique(x), labels = unique(x), cex.axis = 2)
  mtext(text = xlab, side = 1, line = 4, cex = 2, font = 2)
  mtext(text = ylab, side = 2, line = 2.5, cex = 2, font = 2)
  # lines(x = c(0,25.5), y = c(0, .9), lty = 2, lwd = 3)
  if(!is.null(legend.text))
    legend(legend.pos, legend.text, col = cols, title = legend.title, pch = 15, cex = legend.cex)
}


load('181211_simPair_nLoci1000_nInd500_nPerInd2_freqDistrU_nSubPops10_Fst03.RData')
#Reshape my ad hoc data structure from the simulations
tmp <- sapply(h2s.est, function(x){x$h2_add.lmeqtl})
tmp2 <- sapply(h2s.est, function(x){x$h2_add.lmer})
tmp3 <- sapply(h2s.est, function(x){x$h2_slope.lmeqtl})
tmp4 <- sapply(h2s.est, function(x){x$h2_slope.lmer})
tmp5 <- sapply(h2s.est, function(x){x$h2_slope_addData.lmeqtl})
tmp6 <- sapply(h2s.est, function(x){x$h2_slope_addData.lmer})
tmp7 <- sapply(h2s.est, function(x){x$h2_add_polyCorData.lmeqtll}) #misspelled
tmp8 <- sapply(h2s.est, function(x){x$h2_add_polyCorData.lmer})
h2s.true <- seq(.1, .9, .1)
h2s.est.frame <- data.frame(h2.true = rep(h2s.true, each = 10), h2_add.lmeqtl = c(tmp), h2_add.lmer = c(tmp2), h2_slope.lmeqtl = c(tmp3), h2_slope.lmer = c(tmp4), 
                            h2_slope_addData.lmeqtl = c(tmp5), h2_slope_addData.lmer = c(tmp6), h2_add_polyCorData.lmeqtl = c(tmp7), h2_add_polyCorData.lmer = c(tmp8))
#plot additive simulation results
png('issues/05-bivar_vs_randSlope/figures/additiveSim_stratified_nLoci1000_nInd500_nPerInd2.png', width = 1200, height = 800)
plot.contBoxes(x = h2s.est.frame$h2.true, y = h2s.est.frame[, c("h2_add.lmeqtl", "h2_add.lmer", "h2_slope_addData.lmeqtl", "h2_slope_addData.lmer")], legend.text = c('lme4qtl_add', 'lmer_add', 'lme4qtl_cov', 'lmer_cov'), legend.cex = 2.5)
title('Additive Simulation', cex.main = 2)
mtext('Stratified populations', side = 3, cex = 1.5)
dev.off()
#plot poly cor simulation results
png('issues/05-bivar_vs_randSlope/figures/polyCorSim_stratified_nLoci1000_nInd500_nPerInd2.png', width = 1200, height = 800)
plot.contBoxes(x = h2s.est.frame$h2.true, y = h2s.est.frame[, c("h2_slope.lmeqtl", "h2_slope.lmer", "h2_add_polyCorData.lmeqtl", "h2_add_polyCorData.lmer")], legend.text = c('lme4qtl_cov', 'lmer_cov', 'lme4qtl_add', 'lmer_add'), legend.pos = 'bottomright', legend.cex = 2.5)
title('Polygenic Cor Simulation', cex.main = 2)
mtext('Stratified populations', side = 3, cex = 1.5)
dev.off()
