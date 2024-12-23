#' Plot Average Gain/Loss Profiles from Cohort CNA data
#' @param bins
#' @param mar Base graphics margin; argument \code{mar} to \code{par()}.
#' @export
oneD <- function(bins, mar = c(1, 2.5, 1, 2.5)){

   chrsize <- chrSizes()
   offset <- cumsum(chrsize)
   offset <- c(0, offset[seq(length(chrsize)-1)])
   names(offset) <- names(chrsize)
   x <- bins
   x2 <- cbind(x, data.frame(
     x.start = offset[as.character(x$chromosome)] + x$start,
     x.end = offset[as.character(x$chromosome)] + x$end))
   xmax <- max(x2$x.end)
   x2$x.start <- x2$x.start/xmax
   x2$x.end <- x2$x.end/xmax

   par(mar = mar, lwd = 0.5)
   plot(NA, xlim = c(0, 1), ylim = c(-1, 1), xlab='',ylab = '% gain / loss',
        yaxt = 'n',xaxt = 'n' ,mgp = c(1.4,1,0), cex.lab = .8,tck = 0, bty = 'l')
   axis(side = 2, at = seq(-1, 1, 0.5), las = 2, tck = -0.01, cex.axis = 0.6,
        mgp = c(2, 0.5, 0), lwd = 0.0, lwd.ticks = 0.5,
        label = seq(-1, 1, 0.5)*100)
   rect(xleft = offset[seq(2, 22, 2)]/xmax,
        xright = c(offset[seq(3, 21, 2)]/xmax, 1),
        ytop = 1, ybottom = -1, border = NA, col = 'gray95')

   par(lwd = 1.0)
   for(i in seq(nrow(x2))){
     rect(xleft = x2[i, 'x.start'], xright = x2[i, 'x.end'],
          ytop = x2[i,'famp'], ybottom = 0,col = 'red', border = NA)
     rect(xleft = x2[i, 'x.start'], xright = x2[i, 'x.end'],
          ytop = 0, ybottom = -x2[i, 'floss'], col='blue', border = NA)
   }
   for(i in seq_along(offset))
     segments(x0 = offset[i]/xmax, x1 = offset[i]/xmax, y0 = -1, y1 = 1,
              lty = 1,lwd = 0.2)
   x1 <- (offset[22] + chrsize[22])/xmax
   segments(x0 = x1, x1 = x1, y0 = -1,y1 = 1,lty = 1, lwd = 0.2)
   xchr <- (offset[1:21] + offset[2:22])/2/xmax
   xchr <- c(xchr, (2*offset['22'] + chrsize['22'])/2/xmax)
   text(x = xchr[seq(1, 21, 2)], y = -0.8, cex = 0.6,
        label = names(offset)[seq(1, 21, 2)], xpd = NA)
   text(x = xchr[seq(2, 22, 2)], y = -0.9, cex = 0.6,
        label = names(offset)[seq(2, 22, 2)], xpd = NA)
}
