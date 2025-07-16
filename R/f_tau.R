#' f_tau
#'
#' Assuming standard SIR model, the function calculates the value of the density of infection times for all values from zero to Final observation time given the parameters.
#'
#' @import rstan
#' @import deSolve
#' @import ggplot2
#' @param beta Numeric. Infection rate parameter.
#' @param gamma Numeric. Recovery rate parameter.
#' @param rho Numeric. Limiting proportion of initially infected individuals relative to the initially susceptible population. Must satisfy 0 < rho < 1.
#' @param Tmax Numeric. Final observation time of epidemic
#' @param dt Numeric. Time increment.
#' @return  A matrix with time points from zero to final time, and its corresponding value of the density  under standard SIR model at every time point.
#' @export

f_tau = function(beta, gamma, rho, Tmax, dt = 0.1)
{
  time_pts = seq(0, Tmax, by = dt)
  parameters = c(beta, gamma, rho)
  state  = c(S = 1, I = rho, R = 0)


  Lorenz = function(t, state, parameters){
    with(as.list(c(state, parameters)),{
      ## rate of change
      dS = -beta*S*I
      dI = beta*S*I - gamma*I
      dR = gamma*I
      list(c(dS, dI, dR))
    })
  }
  out <- deSolve::ode(y = state, times = time_pts, func = Lorenz, parms = parameters)
  sT = tail(out, 1)[2]
  l = (beta*out[, 2]*out[, 3])/(1 - sT)
  return(matrix(c(time_pts, l), ncol = 2))
}
