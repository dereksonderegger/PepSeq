#' Identify Peaks
#'
#' Given (x,y) pairs of datapoints, calculate peaks.
#'
#' @param x A vector x-values.
#' @param y A vector of y-values.
#' @param method A character string denoting which method to use. Valid options are: PeakSeg and PoT
#' @param param A parameter controling how many peaks are detected. For PoT, it is the treshold.
#' @return A data frame with columns `Peak`, `Start`, and `End`.
#' @export
identify_peaks <- function(x, y, method='PeakSeg', param=NA){
  # if(length(method) > 1){
  #   method = unique(method)
  # }
  # if( length(param) > 1){
  #   param = unique(param)
  # }

  if(method == 'PoT' ){
    foo <-  data.frame(x=x, y=y) %>% arrange(x) %>% mutate(z=row_number() )
    if( is.null(param) | is.na(param) ){
      message('Using default threshold')
      param <- foo %>%
        mutate( q=rank(y) / n(), logq = log(q)   ) %>%
        filter( q >= .85 ) %>%
        summarize(y = min(y)) %>% pull(y)
    }
    out <-
      foo %>%
      filter( y >= param ) %>%
      mutate( delta = diff(c(1, z))) %>%
      mutate( Peak = cumsum(delta) - 1:n() + 1 ) %>%
      group_by(Peak) %>% summarize( Start=min(x), End=max(x) ) %>%
      mutate( Peak = row_number() )

  }else if( method == 'PeakSeg' & is.vector(x) ){
    if( length(unique(y)) <= 2 ){
      out <- data.frame(Peak=NULL, Start=NULL, End=NULL) # not enough unique values
    }else{
      if( is.na(param) ){ param <- 10 }
      foo <- data.frame(count= as.integer(pmax(y,0))) %>%
        mutate(chromStart = as.integer(x), chromEnd   = as.integer(x+1 ))

      # fit <- PeakSegOptimal::PeakSegPDPAchrom(temp, as.integer(20))   # find sequence of best 1 peak, 2 peaks, 3 peaks, etc
      fit <- PeakSegOptimal::PeakSegFPOPchrom(foo, as.integer(param))    # Find the best number of peaks
      out <-
        fit$segments %>%
        filter(status == 'peak') %>%
        rename(Start=chromStart, End=chromEnd) %>%
        mutate(Peak = 1:n()) %>%
        select(Peak, Start, End)
    }
  }

  return(out)
}


