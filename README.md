# DSA.CountData
This package is aim to apply DSA models for count data. The SIR Model or frailty model can be used to infer about the parameters -  infection rate (\eqn{\beta}{beta}) and recovery rate (\eqn{\gamma}{gamma}). It further estimates the limiting proportion of susceptible and infected ( \eqn{\rho}{rho}). DSA method does not necessarily requires the information on the number of initial susceptible. It can also be estimated.
This package contains all R codes requires to impelement the paper.
new_data_ct_lkd() is the main function to apply count likelihood for a new daily count data which simulates HMC posterior samples via stan. We used rstan::version 2.21.
Please look at the Example.R file to see the usage of the important functions.
