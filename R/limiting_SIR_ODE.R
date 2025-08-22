#' limiting_SIR_ODE
#'
#' Assuming standard SIR model, given the parameters, calculates the value of the limiting ODEs.
#'
#' @import rstan
#' @import deSolve
#' @import ggplot2
#' @import mcmcse
#' @param beta Numeric. Infection rate parameter.
#' @param gamma Numeric. Recovery rate parameter.
#' @param rho Numeric. Limiting proportion of initially infected individuals relative to the initially susceptible population. Must satisfy 0 < \eqn{\rho}{rho} < 1.
#' @param Tmax Numeric. Final observation time of epidemic
#' @param dt Numeric. time increment
#' @return Returns the limiting SIR ODEs.
#' @export


limiting_SIR_ODE = function(beta, gamma, rho, Tmax, dt = 0.01)
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
  return(out)
}
