#' Create a Shiny application showing a scatterplot representing a pull-down
#'
#' The pulldown script created by John Altin creates a .csv file
#' with several columns that gives the number of reads for various
#' peptide sequences. This function reads the .csv file and creates
#' a Shiny application that allows the user to explore the pulldown.
#'
#' The input .csv file should have a columns denoting the protein name and location,
#' (by default we assume they are labeled `protein_ID` and `position`) and then
#' a bunch of response columns.  If the response column has a "XXX_cleaved" and
#' "XXX_uncleaved" pair, then we standardize the cleaved/uncleaved pair into a single
#' response.  If the column has no pair, then we don't do any standardization with it.
#'
#' @param input_file The input .csv file. It can be a full path to a file or a
#'                   relative path from the current working directory.
#' @param standardize.method The method by which the cleaved and uncleaved read
#'                        counts are combined. Valid choices are 'additive' or
#'                        'multiplicative'. The default is additive.
#' @param read_indicator An argument that identifies what columns are responses
#'                       to be plotted. This could be either a vector of integers
#'                       or a character string that is at the beginning of all of the
#'                       column names of the response data.
#' @param protein_column A character string indicating which column denotes the protein.
#' @param position_column A character string indicating which column corresponds to the position within a protein.
#' @xamples
#' file <- system.file("extdata", "example_counts.csv", package = "PepSeq")
#' plot_pulldown_Shiny(file, read_indicator='Y_')
#'
#' @export
plot_pulldown_Shiny <- function(input_file,
                                standardize.method='additive',
                                read_indicator='X',
                                protein_column = 'protein_ID',
                                position_column = 'position'){

  input_file = "/Library/Frameworks/R.framework/Versions/3.5/Resources/library/PepSeq/extdata/example_counts.csv"

  df <- read.csv(input_file) %>%
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
    drop_na()                                                                     # Get rid of missing values
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
                                         type = standardize.method)) %>%          ## standardize to combine cleaved/uncleaved values
    select(-cleaved, -uncleaved) %>%
    spread(key=Group, value=signal)

  df <- full_join(df1, df2, by=c('protein_ID', 'position','index')) %>%
    gather('Group', 'signal', -protein_ID, -position, -index)

  Proteins <- df %>%
    group_by(protein_ID) %>%
    summarize(position = min(position),
              index = min(index))


  ymin=floor(   min(df$signal))
  ymax=ceiling( max(df$signal))
  n <- max(df$index)

  # Ctrl <- list(window_start = 1,
  #              window_width = 2000,
  #              protein = Proteins %>% pull(protein_ID) %>%.[1])


  # create the UI
  ui <- fluidPage(
    titlePanel("Peptide Sequence Exploration"),
    sidebarLayout(
      # Sidebar panel for inputs ----
      sidebarPanel(width=2,
                   numericInput(inputId = 'window_width', label = 'Window Width:',
                                value = 2000, min = 1, max=n, step = 100),
                   sliderInput(inputId = 'window_start', label = 'Window Start', min=1,
                               max=n, value = 1, step = 2000),
                   actionButton(inputId = 'window_previous', label = 'Previous'),
                   actionButton(inputId = 'window_next', label = 'Next'),
                   selectInput(inputId = 'protein', label  = 'Protein Selected',
                               choices = levels(df$protein_ID),
                               selected = Proteins %>% pull(protein_ID) %>%.[1] ),
                   numericInput(inputId = 'ymin', label = 'Minimum response:',
                                value = round(min(df$signal)),
                                min = round(min(df$signal)), round(max=max(df$signal))),
                   numericInput(inputId = 'ymax', label = 'Maximum response:',
                                value = round(max(df$signal)),
                                min = round(min(df$signal)), round(max=max(df$signal)))
      ),
      mainPanel(width=10,
                plotOutput(outputId = "Plot")
      )
    )
  )


  # define the server function
  server <- function(input, output, session) {

    # Update the controls as other controls get changed by the user
    observeEvent( input$window_start, {
      # isolate(updateSelectInput(session, 'protein',
      #                   selected = df[input$window_start, 'protein_ID']))
    })
    observeEvent( input$window_next, {
      updateSliderInput(session, 'window_start', value = input$window_start + input$window_width)
    })
    observeEvent( input$window_previous, {
      updateSliderInput(session, 'window_start', value = input$window_start - input$window_width)
    })
    observeEvent( input$protein, {
      isolate(updateSliderInput(session, 'window_start',
                        value=Proteins %>% filter(protein_ID == input$protein) %>% pull(index) ))
    })
    observeEvent( input$window_width, {
      updateSliderInput( session, 'window_start',
                         step = input$window_width)
    })

    output$Plot <- renderPlot({
      xmin <- round(input$window_start )
      xmax <- round(input$window_start + input$window_width)

      df %>%
        filter( index > xmin, index < xmax ) %>%
        ggplot(., aes(x=position, y=signal)) +
        geom_point(size=.2) +
        facet_grid( Group ~ protein_ID, scales='free_x', space='free_x') +
        coord_cartesian(ylim = c(input$ymin, input$ymax))
    })

  }

  shinyApp(ui = ui, server = server)

}
