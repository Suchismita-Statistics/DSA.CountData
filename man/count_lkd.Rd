% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_lkd.R
\name{count_lkd}
\alias{count_lkd}
\title{count_lkd}
\usage{
count_lkd(
  data,
  T.max_analysis = 10,
  frailty = FALSE,
  iteration = 10000,
  num_chain = 4,
  num_cores = 1,
  stan_seed = 12,
  stan_warmup = 0
)
}
\arguments{
\item{data}{Matrix/Data Frame with first column continuous infection times.}

\item{T.max_analysis}{Numeric. Final observation time of epidemic}

\item{frailty}{Logical. If \code{TRUE}, inference is performed using the frailty model; otherwise, the standard SIR model is used.}

\item{iteration}{Integer. The number of iterations for each chain (including warmup). The default is 10000.}

\item{num_chain}{Integer specifying the number of Markov chains. The default is 4.}

\item{num_cores}{Same as "cores" in Sampling function of Stan. Default is 1.}

\item{stan_seed}{Same as "seed" in Sampling function of Stan. Default is 12.}

\item{stan_warmup}{Same as "warmup" in Sampling function of Stan. Default is floor(iteration/2).}
}
\value{
An object of class \code{stanfit} returned by \code{rstan::sampling()},
containing the posterior samples for the parameters \code{beta}, \code{gamma}, \code{rho}, and \code{R_0}.#' @export
}
\description{
This function takes a matrix or data frame whose first column contains exact infection times.
It processes the data to retain only the infection counts and then performs Bayesian inference
on the model parameters using Stan.
}
