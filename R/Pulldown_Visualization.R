#' Create a scatterplot representing a pull-down
#'
#' The pulldown script created by John Altin creates a .csv file
#' with several columns that gives the number of reads for various
#' peptide sequences. This script reads the .csv file and creates
#' a graph using ggplot2 and saves the graph to an output file. The
#' function also invisibly returns the dataset used to create the
#' graph.
#'
#' Critically this function assumes the .csv file has columns for
#' `protein_ID` and `position` and that no other columns other than the read
#' counts start with a `X`.
#'
#' @param input The input data frame read in using the
#'              `import_pulldown()` function.
#'
#' @param output_file The file to save the resulting image. In can be a full
#'                    path to the new file or a relative path from the current
#'                    working directory.
#'
#' @param height The height (in inches) of the resulting file. If NULL, the
#'               default is the number of Group-by-Treatment-by-Rep combinations
#'               plus 3.
#'
#' @param width The width (in inches) of the resulting file. If NULL, the default
#'              is the number of sequences / 1000.
#'
#' @param ymin The minimum value of the response to be shown.  If NULL, all the data is shown.
#'
#' @param ymax The maximum value of the response to be shown.  If NULL, all the data is shown.
#'
#' @param peaks A logical flag denoting if the peaks should be highlighted
#'
#' @param peak_method Which method should be used for peak detection
#'
#' @param peak_param A parameter which controls the peak finding method. For PoT, it is the threshold.
#'
#' @param scales Contols if the y-scale should be 'fixed' across the rows or if the y-scales should be 'free'
#'              to vary across the different rows.
#'
#' @export
plot_pulldown <- function( input, output_file='pulldown.pdf',
                           height=NULL, width=NULL,
                           ymin=NULL, ymax=NULL,
                           peaks=TRUE, peak_method='PoT', peak_param=NA,
                           scales = 'fixed'){

  df = input

  if( is.character(input) ){
    stop('input should be a data frame. Use import_pulldown() to load the data first!')
  }

  if( peaks & any(is.na(peak_param)) ){
    stop('peaks is TRUE, but peak_param is NA and there is no default value')
  }

  # figure out peaks
  if( peaks == TRUE ){
    Peak_Params <- df %>% group_by(Group) %>% count() %>% ungroup() %>% mutate(peak_param=peak_param) %>% select(-n)
    Peaks <- df %>%
      #mutate( signal = signal %>% pmax(signal,0) %>% as.integer() ) %>%
      left_join(Peak_Params, by='Group') %>%
      group_by(Group, protein_ID) %>%
      do( {identify_peaks( .$position, .$signal, method='PoT', .$peak_param[1] ) }) %>%
      left_join( df, by=c('protein_ID', 'Group', 'Start'='position')) %>%
      rename( Start.index = index ) %>% select( Group, protein_ID, Peak, Start, End, Start.index) %>%
      left_join( df, by=c('protein_ID', 'Group', 'End'='position')) %>%
      rename( End.index = index ) %>% select( Group, protein_ID, Peak, Start, End, Start.index, End.index)
    Peaks %>% group_by(Group) %>% count() %>% rename(`Number of Peaks` = n) %>% print()
  }else{
    Peaks <- data.frame(Group=NULL, protein_ID=NULL, Peak=NULL, Start=NULL, End=NULL, Start.index=NULL, End.index=NULL)
  }

  n <- max(df$index)                # For the width of the graph
  p <- length(unique(df$Group))     # For the height of the graph

  if( is.null(height) ){ height=2*p+3 }           # Default values for height/width
  if( is.null(width) ){ width = n/100 }

  if( is.null(ymin) ){ ymin=min(df$signal,na.rm=TRUE) }
  if( is.null(ymax) ){ ymax=max(df$signal, na.rm=TRUE) }

  P <-
    ggplot(df) +
    geom_point(data=df, aes(x=position, y=signal), size=.2) +
    #facet_grid( Group ~ protein_ID, scales='free_x', space='free_x') +
    facet_grid( Group ~ protein_ID, scale='free', space='free_x')

  if( peaks == TRUE & nrow(Peaks) >= 1 ){
    P <- P + geom_rect( data=Peaks, alpha=0.4, aes(xmin=Start, xmax=End, ymin=-Inf, ymax=Inf), fill='salmon')
  }

  if( scales == 'fixed' ){
    P <- P + coord_cartesian(ylim = c(ymin, ymax))
  }

  ggsave(plot=P, filename=output_file, width = width, height=height, limitsize=FALSE)

  return(invisible(P))  # return the plot (invisibly!)
}



# file = system.file('extdata', 'PepSeq_vs_NN.csv', package='PepSeq')
# output_file='pulldown.pdf'
# height=NULL; width=NULL
# ymin=NULL; ymax=NULL
# standardization_method='additive'
# protein_column='Protein_ID'
# position_column='Start_Loc'
# read_indicator=3:6
# peaks=TRUE; peak_method='PoT'; peak_param=c(.5, 95, 10)


# plot_pulldown( file = system.file('extdata', 'PepSeq_vs_NN.csv', package='PepSeq'),
#                output_file='pulldown.pdf',
#                read_indicator = 3:6,
#                protein_column = 'Protein_ID', position_column='Start_Loc',
#                peak_param=c(.5, 95, 10),
#                scales='free')
