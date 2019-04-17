#' Standardize counts relative to uncleaved counts.
#'
#' There are two different ways to handles standardization of
#' cleaved vs uncleaved. For each of the cleaved and uncleaved counts,
#' we first convert to the proportion of reads, and then either look at the
#' difference in proportion (addititive) or the ratio of the proportions
#' (multiplicative). In the addititive case, we might want to normalize
#' the differences based on the differing amounts in the reference library.
#' @param cleaved The cleaved raw counts
#' @param uncleaved The uncleaved raw counts
#' @param ref The raw counts of the reference library. Ideally these should be
#'            identical, but the library likely isn't equally weighted
#' @param type Either `addititive`, `multiplicative`, or `complex`.
#' @examples
#' df <- data.frame( cleaved   = c(20, 10,5),
#'                   uncleaved = c(5, 5, 2),
#'                   ref       = c(3, 3, 2) )
#' with(df, standardize( cleaved, uncleaved, ref ) )
#' with(df, standardize( cleaved, uncleaved, type='multiplicative' ) )
#' @export
standardize <- function(cleaved, uncleaved, ref=NULL, type='additive'){
  if( is.null(ref) ){
    ref = rep(1/length(cleaved), length(cleaved))
  }
  # else{
  # Do some error checking making sure the reference numbers aren't too
  # small.  Zeros would be a huge problem for the reference.
  # }
  if( type == 'additive' ){
    out <- (cleaved/sum(cleaved, na.rm = TRUE) - uncleaved/sum(uncleaved, na.rm=TRUE)) / ( ref / sum(ref) )
  }else if(type == 'multiplicative'){
    out <- (cleaved/sum(cleaved, na.rm=TRUE)) / (uncleaved/sum(uncleaved, na.rm=TRUE))
  }else if(type == 'complex'){
    background_cleaved   <- background_rate(cleaved)
    background_uncleaved <- background_rate(uncleaved)

    cleaved_p   = cleaved   / background_cleaved
    cleaved_p   = cleaved_p / max(cleaved_p)
    uncleaved_p = uncleaved / background_uncleaved
    uncleaved_p = uncleaved_p / max(uncleaved_p)
    diff_p      = cleaved - uncleaved
    diff_p      = diff_p / max(diff_p)
    out = cleaved_p + uncleaved_p + abs(diff_p)
  }else{
    stop("type must be either 'additive', 'multiplicative', or 'complex'")
  }
  return(out)
}

