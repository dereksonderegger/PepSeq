#' Reads an input pulldown .csv file and organizes it for visualization
#'
#' This function is intended to be how users import data from a pulldown.
#'
#' @param file The input .csv, .xls, or .xlsx file. It can be a full path to a file or a
#'             relative path from the current working directory.
#'
#' @param standardization_method The method by which the cleaved and uncleaved read
#'                        counts are combined. Valid choices are 'additive' or
#'                        'multiplicative'. The default is additive.
#'
#' @param scale Should we rescale the signal result so that the maximum signal
#'              is between 100 and 1000.
#'
#' @param read_indicator An argument that identifies what columns are responses
#'                       to be plotted. This could be either a vector of integers
#'                       or a character string (or vector of strings) that is at the beginning of all of the
#'                       column names of the response data.
#'
#' @param protein_column A character string indicating which column denotes the protein.
#'
#' @param position_column A character string indicating which column corresponds to the position within a protein.
#'
#'
#' @export
import_pulldown <- function( file, standardization_method = 'additive', scale=TRUE,
                             read_indicator = 'X',
                             protein_column = 'protein_ID', position_column = 'Peptide.start',
                             Cleaved_Type_Indicators = c('Cl_','Un_') ){

  # Read in the data: from either a .xls or .csv file
  if( is.character(file) ){
    if( str_detect(file, fixed('.csv'))){
      df <- read.csv(file)
    }else{
      df <- readxl::read_xls(file)
    }
  }

  # and get rid of the extraneous columns
  if( is.character( read_indicator )){
    df <- df %>%
      select( protein_column, position_column, starts_with(read_indicator ) )    # specify data columns by name
  }else{
    df <- df %>%                                                                  # specify data columns by locaton
      select( protein_column, position_column, read_indicator  )
  }

  # make the protein and position names consistent
  df <- df %>%
    rename(protein_ID = protein_column, position = position_column)


  # Turn this into a long dataframe
  df <- df %>%
    arrange( protein_ID, position ) %>%                                           # No guarentee user hasn't sent in mixed up rows
    mutate( index = 1:n() ) %>%
    gather(key='Type', value='Value', -protein_ID, -position, -index) %>%         # convert to a _long_ orientation
    drop_na() %>%                                                                 # Get rid of missing values
    mutate( protein_ID = factor(protein_ID),
            protein_ID = fct_reorder(protein_ID, index) )


  # Add a column denoting if the observation is from a cleaved or uncleaved observation.
  Index <- which( str_detect( df$Type, fixed(Cleaved_Type_Indicators[1])) |
                  str_detect( df$Type, fixed(Cleaved_Type_Indicators[2])) )


  df1 <- df[Index, ] %>%
    mutate( Cleave = ifelse( str_detect(.$Type, fixed(Cleaved_Type_Indicators[1])), 'cleaved', 'uncleaved') ) %>%
    mutate( Group = str_remove(.$Type,  fixed(Cleaved_Type_Indicators[1]) )) %>%
    mutate( Group = str_remove(.$Group, fixed(Cleaved_Type_Indicators[2]) )) %>%
    select( index, protein_ID, position, Group, Cleave, Value) %>%
    group_by( index, protein_ID, position, Group ) %>% spread(key=Cleave, value=Value) %>%
    group_by(Group)  %>%
    mutate( signal = PepSeq::standardize(cleaved, uncleaved,                ##
                                         type = standardization_method,     ## standardize to combine cleaved/uncleaved values
                                         scale = scale ) )                  ##


  df2 <- df[-Index, ] %>%
    mutate( cleaved = NA, uncleaved = NA ) %>%
    rename(signal = Value, Group = Type ) %>%
    select( index, protein_ID, position, Group, cleaved, uncleaved, signal)


  out <- bind_rows( df1, df2 )
  return(out)

}




