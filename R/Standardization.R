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
standardize <- function(cleaved, uncleaved, ref=NULL, type='additive'){

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
    mutate(   cleaved =   cleaved / sum(   cleaved, na.rm=TRUE ),
            uncleaved = uncleaved / sum( uncleaved, na.rm=TRUE ) )

  # Now set the background rates for cleaved and uncleaved to be the same.
  # background_scale <-
  #   (df %>% pull(cleaved)   %>% shrink( proportion=.5, side='top') %>% mean() ) /
  #   (df %>% pull(uncleaved) %>% shrink( proportion=.5, side='top') %>% mean() )
  # df <- df %>%
  #   mutate( cleaved = cleaved / background_scale )


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
  # if( scale == TRUE ){
  #   scale <- df %>% select(signal) %>% drop_na() %>% filter( signal < Inf) %>% pull(signal)
  #   scale <- (1 / max(scale)) %>% log10() %>% ceiling() %>% (function(x){10^x})
  #   df$signal <- df$signal * scale * 100
  # }

  return(df$signal)

}




#' Standardize the cleaved and uncleaved as well as create the signal
#'
#' There are several cases where we want to both standardize the cleaved/uncleaved
#' as well as calculate the signal terms
#'
#' @param df A data frame with columns cleaved and uncleaved. If the data frame is
#'           already grouped, then all the standardization occurs within a group.
#' @param type Either `addititive`, `multiplicative`, or `complex`.
#' @param scale Should we scale the signal to have the maximum between 100 and 1000?
#' @param trim_proportion In the rescaling, what percent of the large values should be
#'                        removed to get to a background rate.
#' @return A data frame with columns new columns cleaved_Z, uncleaved_Z, and signal. The rows correspond
#'         to the rows in the input data.frame.
#' @export
full_standardize <- function(df, type='additive', scale=TRUE, trim_proportion=0.25){

  # standardize for the read depth
  df <- df %>%
    mutate(   cleaved_Z =   cleaved / sum(  cleaved, na.rm=TRUE),
            uncleaved_Z = uncleaved / sum(uncleaved, na.rm=TRUE) )

  # Now set the background rates for cleaved and uncleaved to be the same.
  background_scale <- df %>% summarize(
      cleaved_background   =   cleaved_Z %>% shrink( proportion=trim_proportion, side='top') %>% mean(na.rm=TRUE),
      uncleaved_background = uncleaved_Z %>% shrink( proportion=trim_proportion, side='top') %>% mean(na.rm=TRUE) ) %>%
    mutate( background_scale = uncleaved_background / cleaved_background ) %>%
    mutate( background_scale = ifelse(    background_scale == 0,      1, background_scale ),
            background_scale = ifelse( is.nan(background_scale),      1, background_scale ),
            background_scale = ifelse( is.na(background_scale),       1, background_scale ),
            background_scale = ifelse( is.infinite(background_scale), 1, background_scale ))

  df <- df %>% left_join(background_scale, by=group_vars(df) ) %>%
    mutate( cleaved_Z = cleaved_Z * background_scale ) %>%
    select( -cleaved_background, -uncleaved_background, -background_scale )



  # Now make the standardization.
  if( type == 'additive' ){
    df <- df %>% mutate( signal = (cleaved_Z - uncleaved_Z)  )
  }else if(type == 'multiplicative'){
    df <- df %>% mutate( signal = cleaved_Z / uncleaved_Z )
  }else if(type == 'complex'){
    # This is for some experimental work
    df <- df %>% mutate( signal = cleaved_Z - uncleaved_Z )
  }else{
    stop("type must be either 'additive', 'multiplicative', 'complex', or 'none' ")
  }

  # I want the maximum signal to live between 100 and 1000
  if( scale == TRUE ){
    scale_df <-
      df %>% select(group_cols(), signal) %>%
      drop_na() %>% filter( signal < Inf) %>%
      summarize(scale = 1/max(signal) ) %>%
      mutate( scale = (scale %>% log10() %>% ceiling()) )
    df <- df %>% left_join(scale_df, by=group_vars(.) ) %>%
      mutate( signal = signal * 10^scale * 100 ) %>%
      select( -scale )

  }

  return(df)
}

