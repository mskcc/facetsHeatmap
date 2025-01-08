#' Plot Average Gain/Loss Profiles from Cohort CNA data
#' @param bins
#' @param mar Base graphics margin; argument \code{mar} to \code{par()}.
#' @param show.loh Show LOH fractions
#' @export
oneD <- function(bins, mar = c(1.5, 2.5, 1, 2.6), cex = 0.5,
                 show.loh = TRUE){

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

   famp.cut <- 0.2*(floor(max(c(x2$fgain + x2$famp, x2$floh)) / 0.2) + 1)
   floss.cut <- 0.2*(floor(max(x2$floss + x2$fdel) / 0.2) + 1)

   par(mar = mar, lwd = 0.3)
   plot(NA, xlim = c(0, 1), ylim = c(-floss.cut, famp.cut), xlab='',
        ylab = 'Freq. (%)', yaxt = 'n',xaxt = 'n' ,mgp = c(1.4,1,0),
        cex.lab = cex*1.3,tck = 0, bty = 'n')
   tmp <- seq(-floss.cut, famp.cut, 0.2)
   if(length(tmp) > 6)
     tmp <- seq(-floss.cut, famp.cut, 0.4)
   axis(side = 2, las = 2, tck = -0.005, cex.axis = cex,
        mgp = c(2, 0.5, 0), lwd = 0.5, lwd.ticks = 0.5,
        at = tmp,
        label = tmp*100)
   rect(xleft = offset[seq(2, 22, 2)]/xmax,
        xright = c(offset[seq(3, 21, 2)]/xmax, 1),
        ytop = 1, ybottom = -1, border = NA, col = 'gray95')

   par(lwd = 1.0)
   for(i in seq(nrow(x2))){
     rect(xleft = x2[i, 'x.start'], xright = x2[i, 'x.end'],
          ytop = x2[i,'fgain'], ybottom = 0,col = 'orange', border = NA)
     if(x2[i,'famp'] > 0)
       rect(xleft = x2[i, 'x.start'], xright = x2[i, 'x.end'],
            ytop = x2[i,'famp'] + x2[i,'fgain'], ybottom = x2[i,'fgain'],
            col = 'red', border = NA)
     rect(xleft = x2[i, 'x.start'], xright = x2[i, 'x.end'],
          ytop = 0, ybottom = -x2[i, 'floss'], col='cornflowerblue',
          border = NA)
     if(x2[i,'fdel'] > 0)
       rect(xleft = x2[i, 'x.start'], xright = x2[i, 'x.end'],
            ytop = -x2[i,'floss'], ybottom = -x2[i,'floss'] - x2[i,'fdel'],
            col = 'blue', border = NA)
   }
   if(show.loh){
     for(i in seq(nrow(x2))){
       segments(x0 = x2[i, 'x.start'], x1 = x2[i, 'x.end'],
          y0 = x2[i, 'floh'], y1 = x2[i, 'floh'], col = 'forestgreen', lwd = 1)
       if(i == nrow(x2)) break()
       segments(x0 = x2[i, 'x.end'], x1 = x2[i+1, 'x.start'],
          y0 = x2[i, 'floh'], y1 = x2[i+1, 'floh'], col = 'forestgreen', lwd = 1)
     }
   }

   for(i in seq_along(offset))
     segments(x0 = offset[i]/xmax, x1 = offset[i]/xmax, y0 = -1, y1 = 1,
              lty = 1,lwd = 0.2)
   x1 <- (offset[22] + chrsize[22])/xmax
   segments(x0 = x1, x1 = x1, y0 = -1,y1 = 1,lty = 1, lwd = 0.2)
   xchr <- (offset[1:21] + offset[2:22])/2/xmax
   xchr <- c(xchr, (2*offset['22'] + chrsize['22'])/2/xmax)
   text(x = xchr[seq(1, 21, 2)], y = -0.6, cex = .9*cex,
        label = names(offset)[seq(1, 21, 2)], xpd = NA)
   text(x = xchr[seq(2, 22, 2)], y = -0.7, cex = .9*cex,
        label = names(offset)[seq(2, 22, 2)], xpd = NA)
   legend(x = 1.01, y = famp.cut, bty='n',
          fill = c('red','orange', 'cornflowerblue', 'blue'),
          border = c('red', 'orange', 'cornflowerblue', 'blue'),
          legend = c('Amp', 'Gain', 'Loss','Del'), cex = cex*0.7, xpd = NA)
   if(show.loh)
     legend(x = 1.01, y = 0.05, bty='n',
          lty = 1, col = 'forestgreen', lwd = 2, legend = 'LOH', cex = cex*0.7,
          xpd = NA)
}
