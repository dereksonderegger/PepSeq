#' Remove some percentage or number of observations.
#'
#' Often it is desirable to remove a certain number or percentage of a vector.
#' This function will accept either an input number \code{size} or percentage to be
#' removed.  These can be removed from either the start or end of the vector, or
#' the vector can be sorted first and then the items can be removed from the start
#' or end.
#'
#' @param x A vector of numeric values
#' @param proportion The proportion of values to be removed. For example, if we
#'                   wanted to remove 15 percent, the proportion entered would be 0.15.
#' @param size The number of observations to be removed.
#' @param side Where should the values be removed from. \code{start} and \code{end}
#'             remove the observations from the start and end of the vector, while
#'             \code{top} and \code{bottom} will cause the vector to first be sorted
#'             and then remove either the largest or smallest observations.
#' @return A vector that is shorter than the input.
#'
#' @examples
#' set.seed(8675309)
#' input <- rpois(10, lambda=4)
#' input
#' shrink( input, proportion=0.1, side='start')
#' shrink( input, proportion=0.1, side='end')
#' shrink( input, proportion=0.1, side='top')
#' shrink( input, proportion=0.1, side='bottom')
#'
#' @export
shrink <- function(x, proportion=NULL, size=NULL, side='top', na.rm=TRUE){
  if( !is.null(proportion) & !is.null(size) ){
    stop('Either proportion or size can be input, but not both!')
  }

  if( na.rm == TRUE ){
    x <- x[!is.na(x)]
  }

  if( side %in% c('top', 'bottom') ){
    x <- sort(x)
  }

  N <- length(x)
  if( is.null(size) ){
    size <- round( N * proportion )
  }

  if(side %in% c('top','end') ){
    out <- x[1:(N-size)]
  }else{
    out <- x[(size+1):N]
  }

  return(out)
}
