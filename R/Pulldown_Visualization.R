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
#' @param file The input .csv file. It can be a full path to a file or a
#'                   relative path from the current working directory.
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
#' @param read_indicator An argument that identifies what columns are responses
#'                       to be plotted. This could be either a vector of integers
#'                       or a character string that is at the beginning of all of the
#'                       column names of the response data.
#'
#' @param standardization_method The method by which the cleaved and uncleaved read
#'                        counts are combined. Valid choices are 'additive' or
#'                        'multiplicative'. The default is additive.
#'
#' @param protein_column A character string indicating which column denotes the protein.
#'
#' @param position_column A character string indicating which column corresponds to the position within a protein.
#'
#' @export
plot_pulldown <- function( file, output_file='pulldown.pdf',
                           height=NULL, width=NULL,
                           ymin=NULL, ymax=NULL,
                           read_indicator = 'X',
                           standardization_method='additive',
                           protein_column='protein_ID',
                           position_column='position'){
  # import the data
  df <- import_pulldown(file, standardization_method, read_indicator, protein_column, position_column)

  # figure out peaks
  # df <- df %>%
  #   group_by(Group) %>%
  #   mutate( Peak = as.integer(identify_peaks(Value, method='PoT')) ) %>%
  #   mutate( ribbon_ymin = ifelse( Peak == 0, Value, -Inf),
  #           ribbon_ymax = ifelse( Peak == 0, Value,  Inf) )

  # combos <- expand.grid(Group=levels(df$Group),                                   # figure out how many sequences
  #                       Trt=levels(df$Treatment),                                 # and treatment combinations we have
  #                       Rep=levels(df$Rep))                                       # so we know how tall/wide the resulting image
  n <- max(df$index)                                                                # will approximately be
  p <- length(unique(df$Group))

  if( is.null(height) ){ height=p+3 }                                             # Default values for height/width
  if( is.null(width) ){ width = n/100 }

  if( is.null(ymin) ){ ymin=min(df$signal,na.rm=TRUE) }
  if( is.null(ymax) ){ ymax=max(df$signal, na.rm=TRUE) }

  P <- df %>%
    ggplot(., aes(x=position, y=signal)) +
    # geom_ribbon( alpha=0.4, aes(ymin=ribbon_ymin, ymax=ribbon_ymax), fill='salmon') +
    geom_point(size=.2) +
    facet_grid( Group ~ protein_ID, scales='free_x', space='free_x') +
    coord_cartesian(ylim = c(ymin, ymax))

  ggsave(plot=P, filename=output_file, width = width, height=height, limitsize=FALSE)

  invisible(df)  # return the data (invisibly!)
}

