#' Calculate Overlapping Peaks
#'
#' Given two sets of peaks, figure out which peaks overlap and how much
#'
#' @examples
#' A <- data.frame(Start=c(9,20,50,100), End=c(12,25,55,199))
#' B <- data.frame(Start=c(3,24,60,120), End=c( 7,28,63,130))
#' Overlapping_Peaks(A, B)
#'
#' A <- data.frame(
#'   Start = c( 9,20,10,25),
#'   End   = c(12,25,15,35),
#'   Grp   = c('A','A','B','B') )
#' B <- data.frame(
#'   Start = c( 3,24,10,20),
#'   End   = c( 7,28,13,25),
#'   Grp   = c('A','A','B','B') )
#' Overlapping_Peaks(A, B, group_by='Grp')

#' @export
Overlapping_Peaks <- function(A, B, buffer=0, group_by=''){
  A <- A %>% mutate( Start = Start - buffer,
                     End   = End   + buffer)
  B <- B %>% mutate( Start = Start - buffer,
                     End   = End   + buffer)

  A <- A %>% rename(A_Start = Start, A_End = End)
  B <- B %>% rename(B_Start = Start, B_End = End)

  out <-
    merge(A, B) %>% group_by_(group_by) %>%
    mutate( Overlap = pmin(A_End, B_End) - pmax( A_Start, B_Start ) + 1 ) %>%
    filter( Overlap > 0 ) %>%
    group_by_(group_by, 'A_Start', 'A_End') %>% arrange(A_Start, A_End, desc(Overlap)) %>%
    slice(1)

  return( out )
}


