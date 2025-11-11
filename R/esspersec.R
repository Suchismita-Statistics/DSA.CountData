#' esspersec
#' @import rstan
#' @import deSolve
#' @import ggplot2
#' @import mcmcse
#' Returns ESS per second from a Stan object
#'
#' @param result The Stan object
#' @returns ESS per second of the Stan output
#' @export


esspersec = function(result)
{
  res_sum = summary(result)
  ess = res_sum$summary[, "n_eff"]
  time = sum(get_elapsed_time(result)[, 2])
  return(ess / time)
}
