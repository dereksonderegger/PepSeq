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
#' @param scale Should we scale to have the maximum between 100 and 1000?
#' @examples
#' df <- data.frame( cleaved   = c(20, 10,5),
#'                   uncleaved = c(5, 5, 2),
#'                   ref       = c(3, 3, 2) )
#' with(df, standardize( cleaved, uncleaved, ref ) )
#' with(df, standardize( cleaved, uncleaved, type='multiplicative' ) )
#' @export
standardize <- function(cleaved, uncleaved, ref=NULL, type='additive', scale=TRUE){

  # Make a data frame of the input stuff
  df <- data.frame( cleaved=cleaved, uncleaved=uncleaved )

  # Make the reference group default correct
  if( is.null(ref) ){
    df$ref = 1
  }else{
    df$ref = ref
  }

  # Do some error checking making sure the reference numbers aren't too
  # small.  Zeros would be a huge problem for the reference.

  # standardize for the read depth
  df <- df %>%
    mutate( cleaved = cleaved / sum(cleaved, na.rm=TRUE),
            uncleaved = uncleaved / sum(uncleaved, na.rm=TRUE) )

  # Now make the standardization.
  if( type == 'additive' ){
    df <- df %>% mutate( signal = (cleaved - uncleaved) / ref )
  }else if(type == 'multiplicative'){
    df <- df %>% mutate( signal = cleaved / uncleaved )
  }else if(type == 'complex'){
    # This is for some experimental work
    df <- df %>% mutate( signal = cleaved - uncleaved )
  }else{
    stop("type must be either 'additive', 'multiplicative', 'complex', or 'none' ")
  }

  # I want the maximum signal to live between 100 and 1000
  scale <- (1 / max(df$signal, na.rm=TRUE)) %>%
    log10() %>% ceiling() %>% (function(x){10^x})
  df$signal <- df$signal * scale * 100


  return(df$signal)

}




#' Ceiling to nearest
