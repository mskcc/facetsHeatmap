#' Generates Table of Mean 1D CNV Profile
#' @param cncf Input cncf file
#' @param bin.size Size of bin in mb
#' @param clustering.prep Generate chromosome level features for clustering
#' @param progress.bar Display progress bar
#' @return List of components \code{bins} (chromosome, start, and end of bins),
#'         and \code{matrix} (sample vs. bin matrix of copy numbers).
#' @export
#' @examples
#' z <- mat2D(cncf = system.file('extdata','example.cncf', package = 'facetsHeatmap'))
#' head(z)
#'
mat2D <- function(cncf, bin.size = 10, clustering.prep = TRUE,
                  progress.bar = TRUE){

  chrsize <- chrSizes()
  dx <- bin.size * 1e6

  xtot <- sum(chrsize)
  ntot <- round(xtot/dx) + 1
  nchr <- round(chrsize/dx)
# if(sum(nchr)!=ntot) stop('error')
  xchr <- cumsum(nchr)

  chrno <- rep(0,ntot)  # chromosome no. each bin belongs to
  chrno[seq(1,xchr[1])] <- '1'
  for(i in seq(2,length(nchr)))
    chrno[seq(xchr[i-1]+1,xchr[i])] <- as.character(i)

  bins <- data.frame(id=seq(length(chrno)),chromosome=chrno,start=0,end=0)
  bins <- bins[bins$chromosome>0,]
  for(k in seq(as.integer(names(chrsize)))){
    nb <- sum(bins$chromosome==k)
    nk <- floor(chrsize[k]/sum(nb))
    start <- seq(1,chrsize[k]-nk,length.out=nb)
    end <- c(start[-1]-1,chrsize[k])
    bins[bins$chromosome==k,'start'] <- start
    bins[bins$chromosome==k,'end'] <- end
  }

  if(!file.exists(cncf)) stop(paste0(cncf,' does not exist'))
  cncf <- read.table(cncf,header=TRUE,sep='\t')
  sid <- unique(cncf$sid)
  nsamp <- length(unique(cncf$sid))
  bins <- cbind(bins,famp=0, floss=0)

  m <- nrow(bins)
  mat <- matrix(2, nrow = nsamp, ncol = m)
  rownames(mat) <- sid
  colnames(mat) <- bins$id
  if(progress.bar) pb <- txtProgressBar(style=3)

  for(i in seq(m)){
    famp <- floss <- 0
    for(k in seq_along(sid)){
      z <- cncf[cncf$sid==sid[k] & cncf$chrom==bins[i,'chromosome'],]
      flag <- !(z$loc.end < bins[i,'start']) &
        !(bins[i,'end'] < z$loc.start) # some overlap between the two
#       if(sum(z[flag,'tcn.em']>2) > sum(z[flag,'tcn.em']<2))
#         famp <- famp + 1
#       else if(sum(z[flag,'tcn.em']>2) < sum(z[flag,'tcn.em']<2))
#         floss <- floss + 1
      if(any(flag)){
        if(all(z[flag,'tcn.em'] > 2)) famp <- famp + 1
        if(all(z[flag,'tcn.em'] < 2)) floss <- floss + 1
        mat[k, i] <- mean(z[flag, 'tcn.em'])
      }
    }
    bins[i,'famp'] <- famp/nsamp
    bins[i,'floss'] <- floss/nsamp
    if(progress.bar) setTxtProgressBar(pb, i/m)
  }
  if(progress.bar) close(pb)

  if(!clustering.prep)
    w <- list(bins=bins, matrix = mat)
  else{
    m <- length(chrsize)
    dat <- dat2 <- matrix(2, nrow = nsamp, ncol = m)
    rownames(dat) <- rownames(dat2) <- sid
    colnames(dat) <- colnames(dat2) <- names(chrsize)
    for(k in seq_along(sid)){
      for(i in seq(m)){
        z <- cncf[cncf$sid==sid[k] & cncf$chrom==i,]
        dat[k, i] <- round(median(z[, 'tcn.em']))
        dat2[k, i] <- round(median(z[, 'lcn.em'], na.rm = TRUE))
      }
    }
    w <- list(bins=bins, matrix = mat, dat.tcn = dat, dat.lcn = dat2)
  }

  return(w)
}
