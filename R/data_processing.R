#' sellke
#'
#' Generates a dataset using the Sellke construction and then produce a data frame with two columns: days and incidence counts. If \code{nu} is specified, it simulate trajectories incorporating frailty variable.
#' @import rstan
#' @import deSolve
#' @import ggplot2
#' @import mcmcse
#' @param n Interger. Number of initial susceptible.
#' @param beta Numeric. Infection rate parameter.
#' @param gamma Numeric. Recovery rate parameter.
#' @param rho Numeric. Limiting proportion of initially infected individuals relative to the initially susceptible population. Must satisfy 0 < rho < 1.
#' @param Tmax Numeric. Final observation time of epidemic.
#' @param nu Numeric. A parameter denoting the standard deviation of the frailty variable which is assumed to follow Gamma distribution.
#' @return  A data frame with two columns - days and daily infection counts respectively.
#' @export

data_processing = function(n, rho, beta, gamma, T.max_analysis = 10, nu = NULL)
{
  data_Sellke = sellke(n = n, rho = rho, beta = beta, gamma = gamma, Tmax = T.max_analysis, nu = nu)
  initial_sus = data_Sellke[data_Sellke[, 1] != 0, ]
  N = dim(initial_sus)[1]
  M = dim(data_Sellke)[1] - dim(initial_sus)[1]

  duplicates <- duplicated(initial_sus[which(initial_sus[, 1] < T.max_analysis), 1])
  initial_sus = initial_sus[order(initial_sus[, 1]), ]
  infect_during_ep = subset(initial_sus, initial_sus[, 1] < T.max_analysis)
  t = infect_during_ep[, 1]

  infection_days = as.numeric(names(table(floor(t)))) + 1
  infection_count = as.vector(table(floor(t)))
  inf_time_final = numeric(length = T.max_analysis)


  inf_time_final[infection_days] = infection_count
  # temp = matrix(c(0:(T.max_analysis - 1), inf_time_final), ncol = 2)
  df = data.frame(time = 1:(T.max_analysis), cases = inf_time_final)
  return(df)
}
