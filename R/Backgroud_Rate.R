#' Calculate the background rate in a sequence
#'
#' @param x A vector of response values
#' @return A single numeric value
#' @export
background_rate <- function(x){
  probs=seq(.75,.99, by=.005)
  thresholds <- quantile(x, probs)
  # plot(thresholds, probs)
  model <- lm( probs ~ thresholds )
  psi <- segmented(model) %>% summary()
  psi <- psi$psi[1,2]
  return(psi)
}
