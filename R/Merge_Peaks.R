#' Given a data frame of peaks, merge peaks that are close together.
#'
#' @param x A data frame of peaks wich columns `Peak`, `Start`, and `End`.
#' @param min.distance The minimum distance allowed between two peaks.
#' @return A dataframe with as many or fewer peaks.
#'
#' @export
merge_peaks <- function(x, min.distance=5){
  x <- x %>% arrange(Start) %>% mutate(Peak2 = 1:n())
  n <- nrow(x)
  index <- which(x$Start[2:n] - x$End[1:(n-1)] < min.distance )
  out <- x
  for( i in index ){
    out <- out %>% mutate( Peak2 = ifelse( Peak2 == i,  i+1, Peak2 ) )
  }
  out <-
    out %>% group_by(Peak2) %>%
    summarize(Start = min(Start), End=max(End)) %>%
    left_join( select(x, -Start, -End) ) %>%
    select(Peak, Start, End, everything() ) %>%
    select( -Peak2 )
  return(out)
}
