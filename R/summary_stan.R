#' summary_stan
#'
#' Accepting a list of \code{stanfit} objects, it extracts summaries for a specified parameter and return the combined summaries by averaging across all the objects in the list.
#' @import rstan
#' @import deSolve
#' @import ggplot2
#' @import mcmcse
#' @param result_list List of \code{stanfit} objects that we want to summarize.
#' @param index Integer. It takes value 1, 2, 3 for beta, gamma, rho respectively.
#' @param true_values The true value of the parameter.
#' @param not_cap Logical. If \code{TRUE}, it returns the pairs of 95% CI, which do not contain the true parameter. Otherwise, it returns mean of the means, standard deviations and the proportion of 95% CIs captures the true parameter.
#' @export


summary_stan = function(result_list, index, true_value, not_cap = FALSE)
{
  result = lapply(result_list, function(x) extract(x))
  temp1 = lapply(result, function(x) quantile(x[[index]], 0.025))
  temp2 = lapply(result, function(x) quantile(x[[index]], 0.975))
  temp3 = as.numeric( temp1 < rep(true_value, length(lengths(result))) & rep(true_value, length(lengths(result))) < temp2 )

  if(not_cap == TRUE)
  {
    ret =  (1:length(result_list))[which(temp3 == 0)]
    mat = matrix(c(temp1[ret], temp2[ret]), ncol = 2)
    return(list(ret, mat))
    }else{
      mn = mean(foo)
      sd = mean(unlist( lapply(result, function(x) sd(x[[index]])) ))
      covg = mean(temp3)
      return(c(mn, sd, covg))
    }
}
