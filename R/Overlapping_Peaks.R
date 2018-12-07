#' Calculate Overlapping Peaks
#'
#' Given two sets of peaks, figure out which peaks overlap and how much
#'
#' @examples
#' A <- data.frame(Start=c(9,20,50,100), End=c(12,25,55,199))
#' B <- data.frame(Start=c(3,24,60,120), End=c( 7,28,63,130))
#' Overlapping_Peaks(A, B)
#' @export
Overlapping_Peaks <- function(A, B, buffer=0, group_by=''){
  A <- A %>% mutate( Start = Start - buffer,
                     End   = End   + buffer)
  B <- B %>% mutate( Start = Start - buffer,
                     End   = End   + buffer)

  colnames(A) <- paste('A_', colnames(A), sep='')
  colnames(B) <- paste('B_', colnames(B), sep='')

  out <-
    merge(A, B) %>%
    mutate( Overlap = pmin(A_End, B_End) - pmax( A_Start, B_Start ) + 1 ) %>%
    filter( Overlap > 0 ) %>%
    group_by(A_Start, A_End) %>% arrange(A_Start, A_End, desc(Overlap)) %>%
    slice(1)

  return(out)
}


