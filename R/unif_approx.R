#' unif_approx
#'
#' Taking a matrix or data frame with exact infection times in first column, the function prepares a data set keeping only infection count, then do inference of parameters using Stan assuming the infection is uniformly distributed throughout the day.
#' @import rstan
#' @import deSolve
#' @import ggplot2
#' @import mcmcse
#' @param data Matrix/Data Frame with first column continuous infection times.
#' @param T.max_analysis  Numeric. Final observation time of epidemic
#' @param iteration Integer. The number of iterations for each chain (including warmup). The default is 10000.
#' @param num_chain Integer specifying the number of Markov chains. The default is 4.
#' @param num_cores Same as "cores" in Sampling function of Stan. Default is 1.
#' @param stan_seed Same as "seed" in Sampling function of Stan. Default is 12.
#' @param stan_warmup Same as "warmup" in Sampling function of Stan. Default is floor(iteration/2).
#' @return  An object of class \code{stanfit} returned by \code{rstan::sampling()},
#' containing the posterior samples for the parameters \code{beta}, \code{gamma}, \code{rho}, and \code{R_0}.
#' @export


unif_approx = function(data, T.max_analysis, frailty = FALSE,
                       iteration = 1e4, num_chain = 4, num_cores = 1, stan_seed = 12, stan_warmup = 0)
{
  initial_sus = data[data[, 1] != 0, ]
  N = dim(initial_sus)[1]
  M = dim(data)[1] - dim(initial_sus)[1]

  duplicates <- duplicated(initial_sus[which(initial_sus[, 1] < T.max_analysis), 1])
  initial_sus = initial_sus[order(initial_sus[, 1]), ]
  infect_during_ep = subset(initial_sus, initial_sus[, 1] < T.max_analysis)
  t = infect_during_ep[, 1]

  infection_days = as.numeric(names(table(floor(t))))
  infection_count = as.vector(table(floor(t)))
  inf_time_final = numeric(length = T.max_analysis)

  inf_time_final[infection_days+1] = infection_count

  approx_inf_time = list()
  for(i in 1:T.max_analysis)
  {
    approx_inf_time[[i]] = runif(inf_time_final[i], (i-1), i)
  }

  K= nrow(infect_during_ep)
  stan_file = system.file("stan", "unif_apx.stan", package = "DSA.CountData")
  if(stan_warmup == 0)
  {
    stan_warmup = floor(iteration/2)
  }

  sm = rstan::stan_model(file = stan_file)

  data_stan_full = list(N = N, K = K, t0 = 0, infection_time = c(sort(unlist(approx_inf_time)), T.max_analysis))
  fit_full = rstan::sampling(sm_inf_sus, data = data_stan_full, iter = iteration, chain = num_chain, cores = num_cores, seed = stan_seed, warmup = stan_warmup)
  return(fit_full)
}
