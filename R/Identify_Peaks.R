#' Identify Peaks
#'
#' Given an input data frame (or vector) of sequential observations, calculate the peaks
#'
#' @param x A vector or data frame of observations. If a data frame, the rows should correspond
#'          to the sequence of observations, and each column is an independent experiement.
#'
#' @param method A character string denoting which method to use. Valid options are: PeakSeg and PoT
#' @examples
#' lambdas <- rep( c(10,40,10,40,10), times=c(20,10,30,8,30) )
#' data <- rpois( length(lambdas), lambda = lambdas )
#' identify_peaks( x=data )
#' @export
identify_peaks <- function(x, method='PoT', penalty=NULL, threshold=NULL){
  if(method == 'PoT' & is.vector(x) ){
    if( is.null(threshold) ){
      clusters <- kmeans(x, 2)
      index <- which.min( clusters$centers )
      threshold <- (max(x[clusters$cluster == index]) + min(x[clusters$cluster != index])) / 2
    }
    out <- ifelse( x < threshold, 0, 1 )
  }else if( method == 'PeakSeg' & is.vector(x) ){
    model <- PeakSegOptimal::PeakSegFPOP(x, penalty = penalty)
    ends <- model$ends.vec
    ends <- ends[ ends > 0 ] %>% rev() %>% c( length(x) )
    means <- model$mean.vec
    means <- means[ means != Inf ] %>% rev()
    n <- c(ends[1], diff(ends) )

    clusters <- kmeans(means, 2)
    index <- which.min( clusters$centers )
    threshold <- (max(x[clusters$cluster == index]) + min(x[clusters$cluster != index])) / 2

    grp <- ifelse( means > threshold, 1, 0)
    out <- rep( grp, times=n)
  }
  temp <- data.frame(count= as.integer(x)) %>%
    mutate(chromStart = as.integer(1:n() -1),
           chromEnd   = as.integer(chromStart +1 ))
  fit <- PeakSegOptimal::PeakSegPDPAchrom(temp, as.integer(20))

  max.feasible.peaks <- data.table(fit$loss)[feasible==TRUE, max(peaks)]
  show.segments <- data.table(fit$segments)[peaks <= max.feasible.peaks+2]
  show.changes <- show.segments[, data.table(
    position=chromStart[-1]+0.5,
    diff=diff(mean)
  ), by=list(peaks)]
  show.changes[, constraint := ifelse(diff==0, "equality", "inequality")]

  temp %>%
    mutate( x= 1:n()) %>%
    ggplot(aes(x=x, y=count))+
    geom_point()+
    geom_segment(aes(
      chromStart+0.5, mean,
      xend=chromEnd+0.5, yend=mean),
      color="green",
      data=show.segments)+
    scale_linetype_manual(values=c(equality="solid", inequality="dotted"))+
    geom_vline(aes(
      xintercept=position, linetype=constraint),
      color="green",
      data=show.changes)
  #pdf("figure-pepseq-example-mean.pdf", 12, 8)
  print(gg)
  #dev.off()


  # label peaks sequentially
  out <- rle(out)
  index <- which( out$values > 0 )
  out$values[index] <- 1:length(index)
  out <- inverse.rle(out)

  foo <- data.frame( y=x, x = 1:length(x), Peak=factor(out))
  ggplot(foo, aes(x=x, y=y, color=Peak)) + geom_point()

  return(out)
}


