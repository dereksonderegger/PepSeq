#' Import the pulldown metadata
#'
#' Internally, the PepSeq package keeps track of peptides via the protein column and
#' the position within the protein.  However, that isn't consistent with other programs
#' so we need a way to associate the protein/position information with the library_member
#' code, the actual peptide sequence, as well as the Accession number
#'
#' @param file The input .csv, .xls, or .xlsx file. It can be a full path to a file or a
#'             relative path from the current working directory.
#'
#' @param protein_column A character string indicating which column denotes the protein.
#' @param position_column A character string indicating which column corresponds to the position within a protein.
#'
#' @param meta_columns A vector of character strings indicating which columns contain metadata that should be stored.
#'
#' @export
import_pulldown_metadata <- function( file,
    protein_column = 'protein_ID', position_column = 'Peptide.start',
    meta_columns = c('library_member', 'Peptide.seq', 'Accession') ){


  # Read in the data: from either a .xls or .csv file
  if( is.character(file) ){
    if( str_detect(file, fixed('.csv'))){
      df <- read.csv(file)
    }else{
      df <- readxl::read_xls(file)
    }
  }

  # make the protein and position names consistent
  df <- df %>%
    rename(protein_ID = protein_column, position = position_column) %>%
    select( protein_ID, position, meta_columns )

  return(df)

}
