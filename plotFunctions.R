plot.contBoxes <- function(x, y, col = NULL, spacer1 = .01, spacer2 = 3, 
                           xlab = 'Simulated h2', ylab = 'Estimated h2', 
                           legend.text = NULL, legend.title = NULL, legend.pos = 'topleft', 
                           cex.legend = 1, cex.axis = 2, cex.lab = 2, grid = F, ...){
  if(is.null(col)){
    if(!require(wesanderson))
      stop('Could not load package wesanderson')
    cols <- wes_palette('Darjeeling2', ncol(y))
  }
  
  #Deal with NaN:s
  y <- data.frame(
    apply(y, 2, FUN = function(x){
      x[is.nan(x)] <- NA
      x
    })
  )
  
  #Set up grouping variable
  tmp <- seq(from = 0, length.out = ncol(y), by = spacer1)
  tmp <- tmp - mean(tmp)
  grouping <- rep(x, ncol(y)) + rep(tmp, each = length(x))
  #Set up between "group" spacing
  box <- boxplot(unlist(y) ~ grouping, xaxt = 'n', yaxt = 'n', col = cols, plot = F)
  tmp <- rep(seq(0, length.out = length(unique(x)), by = 2), each = ncol(y))
  xPos <- rep(1:ncol(y), length(unique(x))) + tmp*spacer2
  #pos for x-axis
  start <- mean(xPos[1:ncol(y)])
  end <- mean(xPos[(length(xPos) - ncol(y) + 1):length(xPos)])
  xlab.pos <- seq(start, end, length.out = length(unique(x)))
  
  #plot
  boxplot(unlist(y) ~ grouping, xaxt = 'n', at = xPos, yaxt = 'n', col = cols, na.action = na.pass, ...)
  if(grid){
    grid(10,10)
    par(new = T)
    boxplot(unlist(y) ~ grouping, xaxt = 'n', at = xPos, yaxt = 'n', col = cols, na.action = na.pass, ...)
  }
  axis(side = 1, at = xlab.pos, labels = unique(x), las = 2, cex.axis = cex.axis)
  # if(par('usr')[4] < max(x))
    axis(side = 2, cex.axis = cex.axis)
  # else
  #   axis(side = 2, at = unique(x), labels = unique(x), cex.axis = cex.axis)
  mtext(text = xlab, side = 1, line = 4, cex = cex.lab, font = 2)
  mtext(text = ylab, side = 2, line = 2.5, cex = cex.lab, font = 2)
  # lines(x = c(0,25.5), y = c(0, .9), lty = 2, lwd = 3)
  if(!is.null(legend.text))
    legend(legend.pos, legend.text, col = cols, title = legend.title, pch = 15, cex = cex.legend)
}
