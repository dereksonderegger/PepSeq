#' Reads an input .csv file and organizes it for visualization
#'
#' This function is intended to be internal to the package and not user
#' facing.
#' @param file The input .csv file. It can be a full path to a file or a
#'             relative path from the current working directory.
#' @param standardize_method The method by which the cleaved and uncleaved read
#'                        counts are combined. Valid choices are 'additive' or
#'                        'multiplicative'. The default is additive.
#' @param read_indicator An argument that identifies what columns are responses
#'                       to be plotted. This could be either a vector of integers
#'                       or a character string that is at the beginning of all of the
#'                       column names of the response data.
#' @param protein_column A character string indicating which column denotes the protein.
#' @param position_column A character string indicating which column corresponds to the position within a protein.
#' @export
import_pulldown <- function( file,
        standardize_method = 'additive',
        read_indicator = 'X',
        protein_column = 'protein_ID',
        position_column = 'position' ){

  df <- read.csv(file) %>%
    mutate_( protein_ID = protein_column,
             position   = position_column)
  if( is.character( read_indicator )){
    df <- df %>% select( protein_ID, position, starts_with(read_indicator )  )
  }else{
    df <- df %>% select( protein_ID, position,             read_indicator )
  }
  df <- df %>%
    mutate(index = 1:n()) %>%                                                     # create a peptide_ID column
    gather(key='Type', value='Value', -protein_ID, -position, -index) %>%         # convert to a _long_ orientation
    drop_na() %>%                                                                 # Get rid of missing values
    mutate( protein_ID = factor(protein_ID),
            protein_ID = fct_reorder(protein_ID, index) )
  df1 <- df %>%
    filter(!str_detect(Type, fixed('cleave', ignore_case = TRUE))) %>%
    spread(key=Type, value=Value)
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
                                         type = standardize_method)) %>%          ## standardize to combine cleaved/uncleaved values
    select(-cleaved, -uncleaved) %>%
    spread(key=Group, value=signal)

  df <- full_join(df1, df2, by=c('protein_ID', 'position','index')) %>%
    gather('Group', 'signal', -protein_ID, -position, -index)

  return(df)
}
