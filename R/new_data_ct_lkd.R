#'new_data_ct_lkd
#'
#' The function returns a Stan object containing the posterior samples for the parameters beta, gamma, rho ( and N if not available).
#'
#' @import rstan
#' @import deSolve
#' @import ggplot2
#' @import mcmcse
#' @param time_points A numeric vector giving the right endpoints of the time intervals into which the observation period is divided. For example, if intervals are \((0,1], (1,2], \dots\), then \code{time_points = c(1, 2, ...)}. It is assumed that intervals are of equal length and starts from 0.
#' @param infection_count A vector of integers giving the infection counts. We suggest to keep the time points equally spaced between one to the final time and adjust the infection count vector such that it takes value zero, if we do not observe infections in some of the days.
#' @param Final_time Numeric. Final observation time of epidemic.
#' @param initial_sus If available, give the initial number of susceptible, otherwise, put 0.
#' @param frailty Logical. If True, makes inference using frailty model and standard SIR model otherwise.
#' @param iteration Integer. The number of iterations for each chain (including warmup). The default is 10000.
#' @param num_chain Integer specifying the number of Markov chains. The default is 4.
#' @param num_cores Same as "cores" in Sampling function of Stan. Default is 1.
#' @param stan_seed Same as "seed" in Sampling function of Stan. Default is 12.
#' @param stan_warmup Same as "warmup" in Sampling function of Stan. Default is floor(iteration/2).
#' @return  An object of class \code{stanfit} returned by \code{rstan::sampling()},
#' containing the posterior samples for the parameters \eqn{\beta}{beta}, \eqn{\gamma}{gamma}, \eqn{\rho}{rho}.
#' @export


new_data_ct_lkd = function(time_points,
                           infection_count,
                           Final_time,
                           initial_sus = 0,
                           frailty = FALSE,
                           iteration = 1e4,
                           num_chain = 4,
                           num_cores = 1,
                           stan_seed = 12,
                           stan_warmup = 0)
{
  K = sum(infection_count)

  if (Final_time %% 1 != 0)
  {
    Final_time = ceiling(Final_time)
  }
  if (initial_sus == 0)
  {
    if (frailty == FALSE)
    {
      stan_file = system.file("stan", "count_wo_n.stan", package = "DSA.CountData")
    } else{
      stan_file = system.file("stan", "frailty_wo_n.stan", package = "DSA.CountData")
    }
    data_stan_full = list(
      inf_len = length(infection_count),
      infection_count = infection_count,
      K = K,
      t0 = 0,
      int_time = time_points
    )
  } else{
    if (frailty == FALSE)
    {
      stan_file = system.file("stan", "count.stan", package = "DSA.CountData")
    } else{
      stan_file = system.file("stan", "frailty_count.stan", package = "DSA.CountData")
    }

    data_stan_full = list(
      N = initial_sus,
      K = K,
      T_max_int = Final_time,
      infection_count = infection_count,
      t0 = 0,
      int_time = time_points
    )
  }
  if (stan_warmup == 0)
  {
    stan_warmup = floor(iteration / 2)
  }
  sm = rstan::stan_model(file = stan_file)
  fit_full = rstan::sampling(
    sm,
    data = data_stan_full,
    iter = iteration,
    chain = num_chain,
    cores = num_cores,
    seed = stan_seed,
    warmup = stan_warmup
  )
  return(fit_full)
}
