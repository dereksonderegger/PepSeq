#' Calculate MLE Estimates for Gamma-Poisson distribution
#'
#' @param x Input vector of observations
#' @export
Est_GammaPoisson <- function(x){

}

foo <- function(a, b){
  b-a
}


#' Perform Fisher's Chi-squared test on each peptide sequence
#'
#' @param cleaved Vector of cleaved observatinos
#' @param uncleaved Vector of uncleaved observations
#' @export
fisher_chi_sq <- function(cleaved, uncleaved){
  N_c = sum(cleaved)
  N_u = sum(uncleaved)

  data.frame(cleaved, uncleaved) %>%
    mutate(Z = row_number() ) %>%
    group_by(Z) %>%
    mutate( p.value = matrix( c(cleaved, N_c - cleaved, uncleaved, N_u - uncleaved), nrow=2 ) %>%
                      chisq.test() %>% broom::glance() %>% pull(p.value) ) %>%
    pull(p.value)
}
# data.frame(x=c(1,2,2.5, 4), y=6:9) %>% mutate( signal = fisher_chi_sq(x, y) )




x <- c(A = 20, B = 15, C = 25)
chisq.test(x)
chisq.test(as.table(x))
