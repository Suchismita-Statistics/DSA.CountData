#' full_lkd
#'
#' Taking a matrix or data frame with exact infection times, recovery times respectively in two columns, the function do inference of parameters using complete likelihood in Stan.
#'
#' @import rstan
#' @import deSolve
#' @import ggplot2
#' @import mcmcse
#' @param data Matrix/Data Frame with two column continuous infection and recovery times respectively.
#' @param T.max_analysis  Numeric. Final observation time of epidemic
#' @param iteration Integer. The number of iterations for each chain (including warmup). The default is 10000.
#' @param num_chain Integer specifying the number of Markov chains. The default is 4.
#' @param num_cores Same as "cores" in Sampling function of Stan. Default is 1.
#' @param stan_seed Same as "seed" in Sampling function of Stan. Default is 12.
#' @param stan_warmup Same as "warmup" in Sampling function of Stan. Default is floor(iteration/2).

#' @return  An object of class \code{stanfit} returned by \code{rstan::sampling()},
#' containing the posterior samples for the parameters \code{beta}, \code{gamma}, \code{rho}, and \code{R_0}.
#' @export


full_lkd = function(data, T.max_analysis, iteration = 1e4, num_chain = 4, num_cores = 1, stan_seed = 12, stan_warmup = 0)
{
  initial_sus = data[data[, 1] != 0, ]
  N = dim(initial_sus)[1]
  M = dim(data)[1] - dim(initial_sus)[1]

  initial_sus = initial_sus[order(initial_sus[, 1]), ]
  infect_during_ep = subset(initial_sus, initial_sus[, 1] < T.max_analysis)


  t = infect_during_ep[, 1]
  K = nrow(infect_during_ep)

  w = apply(matrix(c(rep(10, K), infect_during_ep[, 2]), ncol = 2), 1, min) - infect_during_ep[, 1]

  L = sum(initial_sus[, 1] < T.max_analysis & initial_sus[, 2] < T.max_analysis)

  ## Epsilons
  initial_inf = data[data[, 1] == 0, ]
  epsilon = initial_inf[which(initial_inf[, 2] < T.max_analysis), 2]
  L_tilde = length(epsilon)
  if(stan_warmup == 0)
  {
    stan_warmup = floor(iteration/2)
  }
  stan_path <- system.file("stan", "full.stan", package = "DSA.CountData")
  sm = rstan::stan_model(file = stan_path)

  data_stan_full = list(N = N, K = K, L = L, M = M, T_max = T.max_analysis,
                        t0 = 0.0, infection_times = c(sort(t), T.max_analysis), w = w, L_tilde = L_tilde, epsilon = epsilon)
  fit_HMC_fulllike_twoeqn = rstan::sampling(sm, data = data_stan_full, iter = iteration,
                                     chain = num_chain, cores = num_cores, seed = stan_seed, warmup = stan_warmup)

  return(fit_HMC_fulllike_twoeqn)
}
