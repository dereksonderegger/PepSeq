#' Identify Peaks
#'
#' Given (x,y) pairs of datapoints, calculate peaks. This is the package interal function that does the
#' heavy calculation. The function `identify_peaks()` is the user friendly version of the function
#' in which the user can just pass the imported data.
#'
#' @param x A vector x-values.
#' @param y A vector of y-values.
#' @param method A character string denoting which method to use. Valid options are: PeakSeg and PoT
#' @param param A parameter controling how many peaks are detected. For PoT, it is the treshold.
#' @return A data frame with columns `Peak`, `Start`, and `End`.
#' @export
identify_peaks_aux <- function(x, y, method='PoT', param=NA, min_peak_length=1, merge_peak_gap=2){

  if(method == 'PoT' ){
    foo <-  data.frame(x=x, y=y) %>% arrange(x) %>% mutate(z=row_number() )
    if( is.null(param) | any(is.na(param)) ){
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
        mutate( Peak = row_number() ) %>%
        select(Peak, Start, End)
    }
  }

  # out <- data.frame( Start = c(10, 40, 60, 72, 82, 100, 120, 140, 147),
  #                    End   = c(20, 45, 70, 80, 90, 110, 120, 145, 155) ) %>%
  #   mutate( Peak = 1:n() )
  # out <- tibble( Start = 1, End=1, Peak=1 ) %>% filter(Peak > 10)

  # merge peaks separated by less than merge_peak_gap
  out <- out %>%
    mutate( gap = c(NA, Start[-1] - End[-n()]),
            merge = ifelse(gap <= merge_peak_gap, TRUE, FALSE),
            merge = ifelse( is.na(merge), FALSE, merge ),
            delta = 1,
            Peak = cumsum( delta*!merge ) ) %>%
    group_by(Peak) %>%
    summarize( Start = min(Start), End = max(End) )

  # Remove all peaks where the peaks are smaller than min_peak_length
  out <- out %>% ungroup() %>%
    mutate( width = End - Start ) %>%
    filter( width >= min_peak_length ) %>%
    mutate( Peak = row_number() ) %>%
    select(Peak, Start, End)

  return(out)
}


#' Identify Peaks
#'
#' Given an input data set, calculate peaks.
#'
#' @param data A data frame of values
#' @param x The column name of the x-axis. Defaults to `position`.
#' @param y The column name of the y-axis. Defaults to `signal`.
#' @param protein The column name of groups that break up the x-axis. Defaults to `protein_ID`.
#' @param Group The column name of the grouping variable for different experiments. Defaults to `Group`.
#' @param method A character string denoting which method to use. Valid options are: PeakSeg and PoT
#' @param param A parameter controling how many peaks are detected. For PoT, it is the treshold.
#' @param min_peak_length The smallest length allowed for a peak
#' @param merge_peak_gap If two peaks are separated for less than this amount, they are merged into
#'                       a single peak.
#'
#' @return A data frame with columns for `Group`, `protein_ID`, `Peak`, `Start`, and `End`.
#' @export
identify_peaks <- function(data, # x='position', y='signal', protein='protein_ID', Group='Group',
                               method='PoT', param=NA, min_peak_length=2, merge_peak_gap=2){

  data <- data %>% mutate( z = row_number() )

  if(method == 'PoT' ){
    # Add a column for the input parameter
    if( is.null(param) | any(is.na(param)) ){
      message('Using default thresholds:')
      thresholds <- data %>%
        group_by( Group, protein_ID ) %>%
        mutate( q=rank(signal) / n(), logq = log(q)   ) %>%
        filter( q >= .95 ) %>%
        summarize(param = min(signal))
      message(thresholds)
      data <- left_join(data, thresholds, by=c('Group', 'protein_ID'))
    }else{
      data <- data %>% mutate(param = param)
    }

    # Now do the filtering
    out <- data %>%
      group_by( Group, protein_ID ) %>%
      filter( signal >= param ) %>%
      mutate( delta = diff(c(1, z))) %>%
      mutate( Peak = cumsum(delta) - 1:n() + 1 ) %>%
      group_by(Group, protein_ID, Peak) %>% summarize( Start=min(position), End=max(position) ) %>%
      mutate( Peak = row_number() )
  }else if( method == 'Z-score' ){
    # Z-score method of Larmen et al.
    # input parameter is the percentage of data to mask
    data %>%
      group_by( Group ) %>%
      mutate( r = rank(signal, na.last = NA),
              r = r / max(r) ) %>%
      filter( r <= param ) %>%
      summarise( alpha = 3, beta= 4)

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
        mutate( Peak = row_number() ) %>%
        select(Peak, Start, End)
    }
  }

  # out <- data.frame( Start = c(10, 40, 60, 72, 82, 100, 120, 140, 147),
  #                    End   = c(20, 45, 70, 80, 90, 110, 120, 145, 155) ) %>%
  #   mutate( Peak = 1:n() )
  # out <- tibble( Start = 1, End=1, Peak=1 ) %>% filter(Peak > 10)

  # merge peaks separated by less than merge_peak_gap
  out <- out %>%
    group_by( Group, protein_ID ) %>%
    mutate( gap = c(NA, Start[-1] - End[-n()]),
            merge = ifelse(gap <= merge_peak_gap, TRUE, FALSE),
            merge = ifelse( is.na(merge), FALSE, merge ),
            delta = 1,
            Peak = cumsum( delta*!merge ) ) %>%
    group_by(Group, protein_ID, Peak) %>%
    summarize( Start = min(Start), End = max(End) )

  # Remove all peaks where the peaks are smaller than min_peak_length
  out <- out %>%
    group_by(Group, protein_ID) %>%
    mutate( width = End - Start ) %>%
    filter( width >= min_peak_length ) %>%
    mutate( Peak = row_number() ) %>%
    select(Group, protein_ID, Peak, Start, End)

  return(out)

}


#' Multivariate standardization using Cleaved/Uncleaved responses
#'
#' Given multiple sequences of (x,y) pairs, create a vector idenifying the peaks.
#'
#' @param data A data frame contain an index, cleaved and uncleaved columns.
#' @param index Which column in the data set corresponds to the index
#' @param cleaved Which column in the data set corresponds to the cleaved
#' @param uncleaved Which column in the data set corresponds to the uncleaved values
#' @export
identify_peaks2 <- function(data, index='index', cleaved='Cleaved', uncleaved='Uncleaved'){
  data <- data %>% rename_('cleaved' = cleaved, 'uncleaved'=uncleaved)
  background_cleaved   <- background_rate(data$cleaved)
  background_uncleaved <- background_rate(data$uncleaved)

  data %>%
    mutate( cleaved_p   = cleaved   / background_cleaved,
            cleaved_p   = cleaved_p / max(cleaved_p),
            uncleaved_p = uncleaved / background_uncleaved,
            uncleaved_p = uncleaved_p / max(uncleaved_p),
            diff_p      = cleaved - uncleaved,
            diff_p      = diff_p / max(diff_p) ) %>%
    mutate( signal = cleaved_p + uncleaved_p + abs(diff_p) ) %>%
    pull(signal) %>% return()
}

