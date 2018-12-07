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
#' @param file The input .csv file. It can be a full path to a file or a
#'                   relative path from the current working directory.
#' @param standardization_method The method by which the cleaved and uncleaved read
#'                        counts are combined. Valid choices are 'additive' or
#'                        'multiplicative'. The default is additive.
#' @param read_indicator An argument that identifies what columns are responses
#'                       to be plotted. This could be either a vector of integers
#'                       or a character string that is at the beginning of all of the
#'                       column names of the response data.
#' @param protein_column A character string indicating which column denotes the protein.
#' @param position_column A character string indicating which column corresponds to the position within a protein.
#' @examples
#' file <- system.file("extdata", "example_counts.csv", package = "PepSeq")
#' plot_pulldown_Shiny(file, read_indicator='Y_')
#'
#' @export
plot_pulldown_Shiny <- function(file,
                                standardization_method='additive',
                                read_indicator='X',
                                protein_column = 'protein_ID',
                                position_column = 'position',
                                peaks=TRUE, peak_method='PoT', peak_param=NA){

  # file = "/Library/Frameworks/R.framework/Versions/3.5/Resources/library/PepSeq/extdata/example_counts.csv"
  # file = '~/Dropbox/NAU/Research/PepSeq/Pulldown Visualization/counts_annotated.csv'
  # peak_method = 'PoT'
  # peak_param  = c(1, NA, NA)

  # Import the data
  df <- PepSeq::import_pulldown(file=file,
                                standardization_method = standardization_method,
                                read_indicator = read_indicator,
                                protein_column=protein_column,
                                position_column=position_column) %>%
    arrange(Group, index)

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
  }else{
    Peaks <- data.frame(Group=NULL, protein_ID=NULL, Peak=NULL, Start=NULL, End=NULL, Start.index=NULL, End.index=NULL)
  }


  Proteins <- df %>%
    group_by(protein_ID) %>%
    summarize(position = min(position),
              index = min(index))


  ymin=floor(   min(df$signal, na.rm=TRUE))
  ymax=ceiling( max(df$signal, na.rm=TRUE))
  n <- max(df$index)


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
                                value = ymin, min = ymin, max=ymax),
                   numericInput(inputId = 'ymax', label = 'Maximum response:',
                                value = ymax, min = ymin, max=ymax)
      ),
      mainPanel(width=10,
                plotOutput(outputId = "Plot")
      )
    )
  )


  # define the server function
  server <- function(input, output, session) {
    # Update the controls as other controls get changed by the user
    # observeEvent( input$window_start, {
    #   # isolate(updateSelectInput(session, 'protein',
    #   #                   selected = df[input$window_start, 'protein_ID']))
    # })
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
    # observeEvent( input$window_width, {
    #   updateSliderInput( session, 'window_start', step = input$window_width)
    # })



    output$Plot <- renderPlot({
      xmin <- round(input$window_start )
      xmax <- round(input$window_start + input$window_width)

      df.small <- df %>% filter( index > xmin, index < xmax )

      P <- ggplot(df.small) +
        geom_point(size=.2, aes(x=position, y=signal)) +
        facet_grid( Group ~ protein_ID, scales='free_x', space='free_x') +
        coord_cartesian(ylim = c(input$ymin, input$ymax))

      if( peaks == TRUE ){
        Peaks.small <- Peaks %>% filter( End.index > xmin, Start.index < xmax )
        P <- P +  geom_rect( data=Peaks.small, alpha=0.4, aes(xmin=Start, xmax=End, ymin=-Inf, ymax=Inf), fill='salmon')
      }

      P
    })

  }

  shinyApp(ui = ui, server = server)

}


