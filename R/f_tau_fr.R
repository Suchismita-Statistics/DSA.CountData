#' f_tau_fr
#'
#' When we assume a frailty model, the function calculates the value of the density of infection times for all values from zero to Final observation time given the value of the parameters.
#'
#' @import rstan
#' @import deSolve
#' @import ggplot2
#' @param beta Numeric. Infection rate parameter.
#' @param gamma Numeric. Recovery rate parameter.
#' @param rho Numeric. Limiting proportion of initially infected individuals relative to the initially susceptible population. Must satisfy 0 < rho < 1.
#' @param nu Numeric. A parameter denoting the standard deviation of the frailty variable which is assumed to follow Gamma distribution.
#' @param Tmax Numeric. Final observation time of epidemic.
#' @param dt Numeric. time increment
#' @return  A matrix with time points from zero to final time, and its corresponding value of the density under frailty model.
#' @export

f_tau_fr = function(beta, gamma, rho, nu, Tmax, dt = 0.1)
{
  time_pts = seq(0, Tmax, by = dt)
  parameters = c(beta, gamma, rho, nu)
  state  = c(S = 1, I = rho, R = 0)


  Lorenz = function(t, state, parameters){
    with(as.list(c(state, parameters)),{
      ## rate of change
      dS = -beta*(S^(1+nu^2))*I
      dI = beta*(S^(1+nu^2))*I - gamma*I
      dR = gamma*I
      list(c(dS, dI, dR))
    })
  }
  out <- deSolve::ode(y = state, times = time_pts, func = Lorenz, parms = parameters)
  sT = tail(out, 1)[2]
  l = (beta*(out[, 2]^(1+nu^2))*out[, 3])/(1 - sT)
  return(matrix(c(time_pts, l), ncol = 2))
}
