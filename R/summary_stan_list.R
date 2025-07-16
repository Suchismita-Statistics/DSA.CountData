#' summary_stan_list
#'
#' It runs summary_stan functions for all four parameters - beta, gamma, rho, R0 and returns in the form of a matrix.
#' @import rstan
#' @import deSolve
#' @import ggplot2
#' @import mcmcse
#' @param result_list List of \code{stanfit} objects that we want to summarize.
#' @param true_values  Vector of true values of the parameters in the order beta, gamma, rho, R0.
#' @param not_cap Logical. If TRUE, it returns the pairs of 95% CI for all four parameters, which do not contain the true parameter. Otherwise, it returns mean of the means, standard deviations and the proportion of 95% CIs captures the true parameter for all four parameters.
#' @export


summary_stan_list = function(result_list, true_value, notcap = FALSE)
{
  para_beta = summary_stan(result_list, index = 1, true_value = true_value[1], not_cap = notcap)
  para_gamma = summary_stan(result_list, index = 2, true_value = true_value[2], not_cap = notcap)
  para_rho = summary_stan(result_list, index = 3, true_value = true_value[3], not_cap = notcap)
  para_R0 = summary_stan(result_list, index = 4, true_value = true_value[4], not_cap = notcap)
  colnam = c("mean", "sd", "95% coverage")
  if(notcap == FALSE){
    mat = matrix(c(para_beta, para_gamma, para_rho, para_R0), ncol = 7, byrow = TRUE)
    colnames(mat) = colnam
    return(mat)
  }else{
    not_captured = list(beta = para_beta, gamma = para_gamma, rho = para_rho, R0 = para_R0)
    return(not_captured)
  }
}
