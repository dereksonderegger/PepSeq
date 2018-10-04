#' Find Variance Change Point
#'
#' I was surprised to find that most all of the the work in detecting changes
#' in varance are almost always associated with time-series objects. There are
#' several packages that address change-in-means situations in regression problems.
#' Given this, we will convert the problem to a change-in-mean by calculating
#' a moving window of standard deviation values and then find the change-point
#' on these.
#'
#' Once the sliding-window standard deviations are calculated, we use a "mob"
#' regression tree that fits a linear model on each leaf is run and we just fit a
#' simple intercept only model. Thus the leaves will have an intercept and
#' different variance terms.
#'
#' To calculate a confidence interval, we do a bootstrap procedure where we
#' bootstrap by resampling the peptides, so all repeats of a peptide is removed
#' included as a group.
#'
#' @param data A data frame where each row is a Peptide and column represents
#'             a repeated measurement.
#' @param window.length The proportion of the data to be used for the sliding
#'                      window.
#' @param interval A TRUE/FALSE flag that indicates if a bootstrap confidence
#'   interval should be calculated.
#'
#' @examples
#' df <- data.frame(x=seq(0,1,by=0.005)) %>%
#'   mutate(y1 = x + rnorm(length(x), sd=ifelse(x<0.5, .1, .3)),
#'          y2 = x + rnorm(length(x), sd=ifelse(x<0.5, .1, .3)),
#'          y3 = x + rnorm(length(x), sd=ifelse(x<0.5, .1, .3)))
#' ggplot(df, aes(x=x)) +
#'   geom_point(aes(y=y1)) +
#'   geom_point(aes(y=y2)) +
#'   geom_point(aes(y=y3))
#'
#' changePoint(df[, 2:4], .01)
#' data = df[,2:4]
#' @export
changePoint <- function(data, window.length=0.01, interval=TRUE){
  nreps <- ncol(data)
  colnames(data) <- paste('Y', 1:nreps, sep='')
  groupsize = floor( nrow(data) * window.length * nreps )

  data2 <- data %>%
    mutate( Peptide = 1:n() ) %>%
    gather('Rep','Value', starts_with('Y') )  %>%  # tidy the data
    group_by(Peptide) %>% mutate( yhat = median(Value) ) %>% group_by() %>%
    arrange(yhat) %>%
    mutate( rownum = 1:n(),
            group = rownum %/% groupsize + 1 ) %>%
    group_by(group) %>%
    summarise(yhat = mean(yhat),
              SD = sd(Value))

  # plot(ggplot(data2, aes(x=yhat, y=SD)) + geom_point())

  # model <- lm( SD ~ 1, data=data2 )
  # model2 <- segmented(model, seg.Z = ~yhat)
  # summary(model2)
  # segmented::slope(model2)
  # return(confint(model2)$yhat)

  # Using a Regression Tree and just a single split
  #model <- rpart::rpart(SD ~ yhat, data=data2)
  #return(model$splits[1, 'index'])

  # Using partykit and a model tree
  model <- partykit::lmtree(SD ~ 1 | yhat, data=data2)
  output <- partykit:::.list.rules.party(model, i = 2) %>%
      stringr::str_sub(start=8) %>%
      str_squish() %>% as.numeric()
  return(output)
}

