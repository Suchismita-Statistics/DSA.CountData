% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doob.R
\name{doob}
\alias{doob}
\title{doob}
\usage{
doob(N, beta, rho, gamma, Tmax)
}
\arguments{
\item{N}{Interger. Number of initial susceptible.}

\item{beta}{Numeric. Infection rate parameter.}

\item{rho}{Numeric. Limiting proportion of initially infected individuals relative to the initially susceptible population. Must satisfy 0 < \code{rho} < 1.}

\item{gamma}{Numeric. Recovery rate parameter.}

\item{Tmax}{Numeric. Final observation time of epidemic.}

\item{dt}{Numeric. Time increment.}
}
\value{
Infection times, recovery times, number of Susceptible, Infected, Recovered at every time point.
}
\description{
Simulates SIR epidemic trajectories using Gillespie Algorithm.
}
