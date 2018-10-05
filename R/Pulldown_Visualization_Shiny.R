#' Create a Shiny application showing a scatterplot representing a pull-down
#'
#' The pulldown script created by John Altin creates a .csv file
#' with several columns that gives the number of reads for various
#' peptide sequences. This function reads the .csv file and creates
#' a Shiny application that allows the user to explore the pulldown.
#'
#' @param input_file The input .csv file. It can be a full path to a file or a
#'                   relative path from the current working directory.
#' @param standardize.method The method by which the cleaved and uncleaved read
#'                        counts are combined. Valid choices are 'additive' or
#'                        'multiplicative'. The default is additive.
#' @param read_indicator What character string indicates the column is data.
#' @export
plot_pulldown_Shiny <- function(input_file,
                                standardize.method='additive',
                                read_indicator='X'){

  df <- read.csv(input_file) %>%
    select( protein_ID, position, starts_with(read_indicator)) %>%                # count column names all start with an X
    mutate(index = 1:n()) %>%                                                     # create a peptide_ID column
    gather(key='Type', value='Value', -protein_ID, -position, -index) %>%         # convert to a _long_ orientation
    drop_na() %>%                                                                 # Get rid of missing values
    separate(Type, c('Group','Cleave','Treatment','Rep'), sep=fixed('_')) %>%     # Break column name into component information
    spread(key=Cleave, value=Value) %>%                                           # Spread CleaveState/Value into Cleaved/Uncleaved
    mutate( signal = PepSeq::standardize(.$Cleaved, .$Uncleaved,                  ##
                                         type = standardize.method)) %>%          ## standardize to combine cleaved/uncleaved values
    mutate( Group = factor(Group),                                                # clean-up in preparation for ggplot facets
            Treatment = factor(Treatment),                                        #
            Rep = factor(Rep) )

  ymin=floor(   min(df$signal))
  ymax=ceiling( max(df$signal))


  # create the UI
  ui <- fluidPage(
    titlePanel("Peptide Sequence Exploration"),
    sidebarLayout(
      # Sidebar panel for inputs ----
      sidebarPanel(
        sliderInput(inputId = "window_width",
                    label = "Window Width:",
                    min = 1,
                    max = n,
                    value = 1000),
        sliderInput(inputId = 'window_center',
                    label = 'Window Center',
                    min=1,
                    max=n,
                    value=500),
        sliderInput(inputId = 'ymin',
                    label = 'y-axis minimum',
                    min=ymin,
                    max=ymax,
                    value=ymin),
        sliderInput(inputId = 'ymax',
                    label = 'y-axis maximum',
                    min=ymin,
                    max=ymax,
                    value=ymax)
      ),
      mainPanel(
        plotOutput(outputId = "Plot")
      )
    )
  )


  # define the server function
  server <- function(input, output) {

    output$Plot <- renderPlot({
      xmin <- round(input$window_center - input$window_width/2)
      xmax <- round(input$window_center + input$window_width/2)

      p <- df %>%
        filter( index > xmin, index < xmax ) %>%
        ggplot(., aes(x=position, y=signal)) +
        geom_point(size=.2) +
        facet_grid( Group*Treatment*Rep ~ protein_ID, scales='free_x', space='free_x') +
        coord_cartesian(ylim = c(input$ymin, input$ymax))
      print(p)

    })

  }

  shinyApp(ui = ui, server = server)

}
