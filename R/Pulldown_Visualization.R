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
#' @param input_file The input .csv file. It can be a full path to a file or a
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
#' @param standize.method The method by which the cleaved and uncleaved read
#'                        counts are combined. Valid choices are 'additive' or
#'                        'multiplicative'. The default is additive.
#'
#' @export
plot_pulldown <- function( input_file, output_file='pulldown.pdf',
                           height=NULL, width=NULL,
                           standardize.method='additive' ){

  df <- read.csv(input_file) %>%
    select( protein_ID, position, starts_with('X')) %>%                           # count column names all start with an X
    gather(key='Type', value='Value', -protein_ID, -position) %>%                 # convert to a _long_ orientation
    drop_na() %>%                                                                 # Get rid of missing values
    separate(Type, c('Group','Cleave','Treatment','Rep'), sep=fixed('_')) %>%     # Break column name into component information
    spread(key=Cleave, value=Value) %>%                                           # Spread CleaveState/Value into Cleaved/Uncleaved
    mutate( signal = PepSeq::standardize(.$Cleaved, .$Uncleaved,                  ##
                                          type = standardize.method)) %>%         ## standardize to combine cleaved/uncleaved values
    mutate( Group = factor(Group),                                                # clean-up in preparation for ggplot facets
            Treatment = factor(Treatment),                                        #
            Rep = factor(Rep) )

  combos <- expand.grid(Group=levels(df$Group),                                   # figure out how many sequences
                        Trt=levels(df$Treatment),                                 # and treatment combinations we have
                        Rep=levels(df$Rep))                                       # so we know how tall/wide the resulting image
  n <- nrow(df) / nrow(combos)                                                    # will approximately be
  p <- nrow(combos)

  if( is.null(height) ){ height=p+3 }                                             # Default values for height/width
  if( is.null(width) ){ width = n/100 }

  df %>%
    ggplot(., aes(x=position, y=signal)) +
    geom_point(size=.2) +
    facet_grid( Group*Treatment*Rep ~ protein_ID, scales='free', space='free_x') +
    ggsave(output_file, width = width, height=height, limitsize=FALSE)

  invisible(df)  # return the data (invisibly!)
}

