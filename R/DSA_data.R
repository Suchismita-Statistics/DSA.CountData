#' DSA_data
#'
#' Generates a data set of exact infection and recovery times according to DSA approximation for simulating epidemic trajectories.
#' @import rstan
#' @import deSolve
#' @import ggplot2
#' @param n Interger. Number of initial susceptible.
#' @param beta Numeric. Infection rate parameter.
#' @param gamma Numeric. Recovery rate parameter.
#' @param rho Numeric. Limiting proportion of initially infected individuals relative to the initially susceptible population. Must satisfy 0 < rho < 1.
#' @param Tmax Numeric. Final observation time of epidemic
#' @param dt Numeric. time increment. Default is 0.1.
#' @return  A data set with two columns - exact infection times and recovery times respectively.
#' @export

DSA_data = function(n, beta, gamma, rho, Tmax, dt = 0.1)
{
  m = n * rho
  data_gen = matrix(0, nrow = (n + m), ncol = 2)
  time_pts = seq(0, Tmax, by = dt)
  parameters = c(beta, gamma, rho)
  state  = c(S = 1, I = rho, R = 0)


  Lorenz = function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      ## rate of change
      dS = -beta * S * I
      dI = beta * S * I - gamma * I
      dR = gamma * I
      list(c(dS, dI, dR))
    })
  }
  out <- ode(
    y = state,
    times = time_pts,
    func = Lorenz,
    parms = parameters
  )
  s_T = tail(out, 1)[2]

  u = runif(n)
  not_inf_during_epi = sum(as.numeric(u < s_T))

  u_infect = runif((n - not_inf_during_epi))
  temp = u_infect * (1 - s_T) + s_T
  samples = approx(x = out[, 2], y = out[, 1], xout = temp)

  infect_times = samples$y
  infect_times[is.na(infect_times)] = Tmax
  # sum(as.numeric(infect_times>10))
  data_gen[, 1] = sort(c(rep(0, m), infect_times, rep(Tmax, not_inf_during_epi)))

  data_gen[, 2] = data_gen[, 1] + rexp((n + m), rate = gamma)

  return(data_gen)
}
