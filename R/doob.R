#' doob
#'
#' Simulates SIR epidemic trajectories using Gillespie Algorithm.
#' @import rstan
#' @import deSolve
#' @import ggplot2
#' @param N Interger. Number of initial susceptible.
#' @param beta Numeric. Infection rate parameter.
#' @param gamma Numeric. Recovery rate parameter.
#' @param rho Numeric. Limiting proportion of initially infected individuals relative to the initially susceptible population. Must satisfy 0 <  \eqn{\rho}{rho} < 1.
#' @param Tmax Numeric. Final observation time of epidemic.
#' @param dt Numeric. Time increment.
#' @return  Infection times, recovery times, number of Susceptible, Infected, Recovered at every time point.
#' @export


doob = function(N, beta, rho, gamma, Tmax)
{
  S = numeric()
  I = numeric()
  R = numeric()
  infection_time = numeric()
  recovery_time = numeric()
  time_tracker = numeric()

  S[1] = N
  I[1] = round(N*rho)
  R[1] = 0

  lambda_I = numeric()
  lambda_I[1] = 0
  lambda_R = numeric()
  lambda_R[1] = 0

  k = 1
  time_tracker[1] = 0
  h = 1
  while (tail(time_tracker, 1) < Tmax && S[k] > 0 && I[k] > 0) {
    lambda_I[k+1] = beta*S[k]*I[k]/N
    lambda_R[k+1] = gamma*I[k]

    lambda_s = lambda_I[k+1] + lambda_R[k+1]
    if (lambda_s == 0) break

    time_tracker[k + 1] = time_tracker[k] + rexp(1, lambda_s)
    p = lambda_I[k+1]/lambda_s
    inf_or_rec = rbinom(1, size = 1, prob = p)
    #print(frailty)

    if(inf_or_rec == 1)
    {
      S[k+1] = S[k] - 1
      I[k+1] = I[k] + 1
      R[k+1] = R[k]

      #frailty = frailty[-1]
      h = h+1
      infection_time = c(infection_time, time_tracker[k+1])
    }else{
      S[k+1] = S[k]
      I[k+1] = I[k] - 1
      R[k+1] = R[k] + 1
      recovery_time = c(recovery_time, time_tracker[k+1])
    }
    k = k+1
    # print(k)
  }
  return(list(infection_time = infection_time, recovery_time = recovery_time, S = S, I = I, R = R))
}
