#' file = system.file('extdata', 'example_counts.csv', package='PepSeq')
#' standardization_method='additive'
#' read_indicator='Y'
#' protein_column = 'protein_ID'
#' position_column = 'position'
#' peak_method = 'PeakSeg'
#' peak_param = NA

# df <- import_pulldown( file=system.file('extdata', 'example_counts.csv', package='PepSeq'),
#                        read_indicator = 'Y', protein_column = 'protein_ID',
#                        position_column = 'position' )

plot_pulldown_Shiny(
  file = system.file('extdata', 'PepSeq_vs_NN.csv', package='PepSeq'),
  read_indicator = 3:6,
  protein_column = 'Protein_ID', position_column='Start_Loc',
  peak_param=c(.5, 95, 10))

plot_pulldown_Shiny(
  file = system.file('extdata', 'PepSeq_vs_NN.csv', package='PepSeq'),
  read_indicator = 3:6,
  protein_column = 'Protein_ID', position_column='Start_Loc',
  peaks=FALSE, peak_param=c(.5, 95, 10))

df %>% group_by(index) %>% count() %>% pull(n) %>% hist()

df %>% filter( protein_ID == 'AAV34704.1') %>%
  arrange(index) %>%
  View()

df %>% filter( protein_ID == 'AAV34704.1') %>%
  filter(position == 555)

  arrange(index) %>%
  View()


  pull(index) %>% unique()

unique(df$protein_ID)

out %>% filter( protein_ID == 'AAV34704.1') %>%
  filter(position >550)
