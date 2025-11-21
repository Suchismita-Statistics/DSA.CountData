#' sellke
#'
#' Generates a dataset using the Sellke construction for simulating epidemic trajectories. If \code{nu} is specified, it simulate trajectories incorporating frailty variable.
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
#' @return  A data set with two columns - exact infection times and recovery times respectively.
#' @export


sellke = function(n, beta, gamma, rho, Tmax, nu = NULL)
{
  if(is.null(nu)){
    Q = sort(rexp(n, rate = 1))
  }else{
    u = rgamma(n, shape = 1/nu^2, rate = 1/nu^2)
    Q = sort(rexp(n, rate = u))
  }
  m = round(n * rho)
  mattime = matrix(Tmax, ncol = 2, nrow = (n + m))

  mattime[1:m, 1] = rep(0, m)
  mattime[1:m, 2] = rexp(m, rate = gamma)

  S = numeric()
  I = numeric()
  R = numeric()
  X = numeric()

  time_track = numeric()

  S[1] = n
  I[1] = m
  R[1] = 0

  k = 1
  h = 1
  time_track[1] = 0
  Lambda_past = 0


  while (I[k] > 0 && time_track[k] < Tmax && S[k] > 0)
  {
    min_rec = min(subset(mattime[, 2], mattime[, 2] > time_track[k]))
    Lambda_t = Lambda_past + beta * I[k] * (min_rec - time_track[k]) / n

    if (Q[h] < Lambda_t)
    {
      time_track[k + 1] = time_track[k] + (Q[h] - Lambda_past) / (beta * I[k] /
                                                                    n)
      S[k + 1] = S[k] - 1
      I[k + 1] = I[k] + 1
      R[k + 1] = R[k]
      mattime[m + h, 1] = time_track[k + 1]
      mattime[m + h, 2] = mattime[m + h, 1] + rexp(1, gamma)
      min_rec = min(subset(mattime[, 2], mattime[, 2] > time_track[k + 1]))
      Lambda_past = Q[h]
      h = h + 1
      X = c(X, "inf")
    } else{
      time_track[k + 1] = min_rec
      S[k + 1] = S[k]
      I[k + 1] = I[k] - 1
      R[k + 1] = R[k] + 1

      min_rec = min(subset(mattime[, 2], mattime[, 2] > time_track[k + 1]))
      X = c(X, "rec")
      Lambda_past = Lambda_t
    }
    k = k + 1
  }
  return(mattime)
}
