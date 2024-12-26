#' Generate heatmap of Copy-Number Profiles
#' @param z List of components \code{bins}, \code{matrix}, \code{dat.tcn},
#'          and \code{dat.lcn}
#' @param clusterZ Cluster samples using abbreviated CNA features
#' @param K Number of clusters
#' @param reorder.hc Reorder samples by hierarchical clustering after grouping by
#'        cluster
#' @param tcn.max Maximum total copy number to display
#' @param relative Use relative copy numbers
#' @param tcn.min If \code{Relative}, lower bound of relative tcn loss to show
#' @param mar Base graphics margin; argument \code{mar} to \code{par()}
#' @param legend.title Tile of legend
#' @export
showHeatmap <- function(z, clusterZ = TRUE, K = 2, reorder.hc = TRUE,
                        tcn.max = 6, relative = FALSE,
                        tcn.min = -2,
                        margin=c(0.1, 0.5, 0.5,1.3),
                        mar = c(1, 2.5, 1, 2.6),
                        cex = 0.6,
                        legend.title = 'tcn'){

  sid <- rownames(z$matrix)
  nsamp <- length(sid)
  if(clusterZ){
    dat <- clusterZ(z, K = K)
    cid <- dat[sid,'cluster']
  } else{
    cid <- rep(1, nsamp)
  }
  names(cid) <- sid

  sid1 <- NULL
  for(k in seq(K)){
    sidk <- sid[cid == k]
    if(reorder.hc){
      hc <- hclust(dist(z$matrix[sidk,]))$order
      sidk <- sidk[hc]
    }
    sid1 <- c(sid1, sidk)
  }
  sid <- sid1
  cid <- cid[sid]
  x <- z$matrix[sid, ]
  if(relative)
    x <- x - rowMeans(x)
  bins <- z$bins
  xb <- bins[!duplicated(bins[, 2]), 'id']
  xb <- c(xb, nrow(bins))

  dat <- cbind(data.frame(sid = sid), x)
  w <- data.frame(tidyr::pivot_longer(dat, cols = colnames(x)))
  w$sid <- factor(w$sid, levels = sid)
  w$name <- factor(w$name, levels = colnames(x))

  wmax <- tcn.max
  w$value[w$value > wmax] <- wmax
  wmin <- ifelse(relative, tcn.min, 0)
  w$value[w$value < wmin] <- wmin

  if(relative){
    np1 <- ceiling(wmax*10)
    np2 <- ceiling(abs(wmin)*10)
    legend.title <- expression(Delta*' tcn')
  } else{
    np1 <- (ceiling(wmax) - 2)*10
    np2 <- 20
  }
  pal1 <- colorRampPalette(RColorBrewer::brewer.pal(9,'OrRd'))(np1)
  pal2 <- colorRampPalette(RColorBrewer::brewer.pal(9,'Blues'))(np2)
  pal <- c(rev(pal1),'white',pal2)

  a <- cumsum(as.numeric(table(cid)))
  label <- paste('C', seq(K))

  x0 <- 0
  ticks <- c(rev(seq(0,wmin, by = -2)), seq(2, wmax, by = 2))
  p2d <- ggplot2::ggplot(w, ggplot2::aes(x = name, y = sid, color = value,
                                         fill = value)) +
    ggplot2::xlab('') +
    ggplot2::ylab('') +
    ggplot2::labs(fill = legend.title) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colors = rev(pal),
                  na.value = 'white', breaks = ticks) +
    ggplot2::scale_color_gradientn(colors = rev(pal),
                  na.value = 'white', breaks = ticks) +
    ggplot2::geom_vline(xintercept = c(xb[seq(22)] - 0.5, nrow(bins) + 0.5),
                        linewidth = 0.2) +
    ggplot2:: guides(color = 'none') +
    ggplot2::annotate('text', x = x0 +
        (xb[seq(1, 21, 2)] + xb[seq(2, 22, 2)])/2,
                      y = -0.008*nsamp, label = seq(1, 21, 2), size = cex*4.5) +
    ggplot2::annotate('text', x = x0 +
        (xb[seq(2, 22, 2)] + xb[seq(3, 23, 2)])/2,
                      y = -0.015*nsamp, label = seq(2, 22, 2), size = cex*4.5) +
    ggplot2::coord_cartesian(xlim = c(x0, ncol(dat)), ylim = c(1, nrow(x)),
                             clip = 'off') +
    ggplot2::annotate('text', x = -1, y = nsamp - c(0,a[seq(1,length(a)-1)]) - 1,
                      label = label, size = cex*5,
                      hjust = 1, color = 'black') +
    ggplot2::geom_hline(yintercept = c(nsamp + 0.5, nsamp - a + 0.5),
      linewidth = 0.3) +
    ggplot2::scale_y_discrete(limits = rev(sid)) +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          plot.margin = ggplot2::margin(margin, unit = 'cm'),
          legend.key.size = ggplot2::unit(0.3, 'cm'),
          panel.grid = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          legend.text = ggplot2::element_text(size = cex*15),
          legend.title = ggplot2::element_text(size = cex*17))

  dummy <- function(){
    oneD(z$bins, mar = mar, cex = cex)
  }
  plt2 <- cowplot::plot_grid(dummy, p2d, nrow = 2, ncol = 1,
                             rel_heights = c(0.1, 0.9))
  return(plt2)
}
