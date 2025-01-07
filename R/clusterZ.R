#' Cluster Samples based on Abbreviated Copy-Number Profiles
#' @param z List of components \code{bins}, \code{matrix}, \code{dat.tcn},
#'          and \code{dat.lcn}
#' @param K (Range of) number of clusters to consider
#' @return If \code{K} is of length more than 1, ggplot object of \code{K}
#'        vs. BIC. If \code{K} is a single value, data frame of cluster IDs
#'        for samples with posterior prob.
#' @export
clusterZ <- function(z, K = 2, plotBIC = TRUE){

   A <- data.frame(apply(z$dat.tcn, 1:2, function(x){if(is.na(x)) return(1);
                                                     as.numeric(x > 2) + 1}))
   L <- data.frame(apply(z$dat.lcn, 1:2, function(x){if(is.na(x)) return(1);
                                                     as.numeric(x > 0) + 1}))
   colnames(L) <- paste0('Y', seq(ncol(L)))
   x <- cbind(A, L)
   f <- cbind(
     X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,
     X13,X14,X15,X16,X17,X18,X19,X20,X21,X22,
     Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,Y11,Y12,
     Y13,Y14,Y15,Y16,Y17,Y18,Y19,Y20,Y21,Y22) ~ 1

   nsamp <- NROW(z$matrix)
   K <- min(c(K, nsamp)) # cannot cluster into more than nsamp dimensions
   fit <- list()
   bic <- NULL
   for(k in K){
     fit[[as.character(k)]] <- poLCA::poLCA(f, data=x, nclass=k, verbose = FALSE)
     bic <- rbind(bic,data.frame(k=k,bic=fit[[as.character(k)]]$bic))
   }
   plt <- NULL
   if(length(K) > 1){
     if(plotBIC){
       plt <- ggplot2::ggplot(bic, ggplot2::aes(x = k, y = bic)) +
               ggplot2::theme_linedraw() +
               ggplot2::geom_path() +
               ggplot2::geom_point() +
               ggplot2::ylab('BIC') +
               ggplot2::xlab('K')
       return(plt)
     } else{
       return(bic)
     }
   }

   fk <- fit[[1]]
   pred <- fk$predclass

   dat <- data.frame(posterior = fk$posterior, cluster = pred)
   rownames(dat) <- rownames(fk$x)

   ave.tcn <- rep(2, K)
   for(k in seq(K))
     ave.tcn[k] <- mean(z$matrix[rownames(dat)[dat$cluster==k],])
   idx <- match(seq(K), order(ave.tcn, decreasing = FALSE))
   dat$cluster <- idx[dat$cluster]

   return(dat)
}
