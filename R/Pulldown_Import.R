#' Reads an input .csv file and organizes it for visualization
#'
#' This function is intended to be internal to the package and not user
#' facing.
#' @param file The input .csv file. It can be a full path to a file or a
#'             relative path from the current working directory.
#' @param standardization_method The method by which the cleaved and uncleaved read
#'                        counts are combined. Valid choices are 'additive' or
#'                        'multiplicative'. The default is additive.
#' @param read_indicator An argument that identifies what columns are responses
#'                       to be plotted. This could be either a vector of integers
#'                       or a character string (or vector of strings) that is at the beginning of all of the
#'                       column names of the response data.
#' @param protein_column A character string indicating which column denotes the protein.
#' @param position_column A character string indicating which column corresponds to the position within a protein.
#' @export
import_pulldown <- function( file, standardization_method = 'additive',
                             read_indicator = 'X', protein_column = 'protein_ID',
                             position_column = 'position' ){

  df <- read.csv(file) %>%
    rename_( protein_ID = protein_column,
             position   = position_column)
  if( is.character( read_indicator )){
    df <- df %>% select( protein_ID, position, starts_with(read_indicator )  )    # specify data columns by name
  }else{
    df <- df %>% select( protein_ID, position,             read_indicator )       # specify data columns by locaton
  }
  df <- df %>%
    arrange( protein_ID, position ) %>%                                           # No guarentee user hasn't sent in mixed up rows
    mutate(index = 1:n()) %>%                                                     # create a peptide_ID column
    gather(key='Type', value='Value', -protein_ID, -position, -index) %>%         # convert to a _long_ orientation
    drop_na() %>%                                                                 # Get rid of missing values
    mutate( protein_ID = factor(protein_ID),
            protein_ID = fct_reorder(protein_ID, index) )
  df1 <- df %>%
    filter(!str_detect(Type, fixed('cleave', ignore_case = TRUE))) %>%
    rename(Group = Type, signal = Value) %>%
    mutate(cleaved=NA, uncleaved=NA) %>%
    select(protein_ID, position, index, Group, cleaved, uncleaved, signal)
  df2 <-
    df %>% filter( str_detect(Type, fixed('cleave', ignore_case = TRUE))) %>%
    mutate( Type = stringi::stri_reverse(Type) ) %>%                              # reversing because cleaved/uncleaved is at the end
    separate(Type, c('Cleave', 'Group'), sep=fixed('_'), extra='merge') %>%       # split on the first '_'
    mutate(Cleave = stringi::stri_reverse(Cleave),                                # undo my reversing operation
           Group  = stringi::stri_reverse(Group)) %>%
    arrange(Group, protein_ID, position, Cleave) %>%
    mutate(Cleave = str_to_lower(Cleave)) %>%                                     # get rid of non-consistent capitalizations
    spread(key=Cleave, value=Value) %>%
    mutate( signal = PepSeq::standardize(.$cleaved, .$uncleaved,                  ##
                                         type = standardization_method)) #%>%      ## standardize to combine cleaved/uncleaved values
    #select(-cleaved, -uncleaved) %>%
    #spread(key=Group, value=signal)

  out <- rbind( df1, df2 )

  return(out)
}



# df <- import_pulldown( file=system.file('extdata', 'example_counts.csv', package='PepSeq'),
#                        read_indicator = 'Y', protein_column = 'protein_ID',
#                        position_column = 'position' )
# df <- import_pulldown(
#   file=system.file('extdata', 'PepSeq_vs_NN.csv', package='PepSeq'),
#   read_indicator = 3:6,
#   protein_column = 'Protein_ID', position_column = 'position')

