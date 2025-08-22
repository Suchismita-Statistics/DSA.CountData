## An example of running PDSIR for simulation
## Use the package: https://github.com/rmorsomme/PDSIR

## # install.packages("devtools")
#devtools::install_github("rmorsomme/PDSIR")

library(PDSIR)

## This function make inference using PDSIR method using one chain.
## data_Sellke: a matrix or data frame of two columns with continuous infection and recovery times.
## rho: Proportion of observations we want to update; Please look the paper for further details.
## Tmax: Final observation time of epidemic.
## iter: Number of iteration of the chain.
## b1, g1: Initial choice of parameter beta and gamma.


PDSIR_calc = function(data_Sellke, rho, Tmax = 10, iter = 5e4, b1 = 1, g1 = 1)
{
  initial_sus = data_Sellke[data_Sellke[, 1] != 0, ]
  N = dim(initial_sus)[1]
  M = dim(data_Sellke)[1] - dim(initial_sus)[1]


  initial_sus = initial_sus[order(initial_sus[, 1]), ]
  infect_during_ep = subset(initial_sus, initial_sus[, 1] < Tmax)
  t = infect_during_ep[, 1]

  infection_count = as.vector(table(floor(t)))
  day_num = as.numeric(names(table(floor(t))))
  inf_ct_final = numeric(length = (Tmax) )
  inf_ct_final[1+day_num] = infection_count

  Y = list(I_k = inf_ct_final, I0 = M, S0 = N, ts = 0:Tmax, t_end = Tmax)
  R1 = b1/g1
  theta1 <- list(
    R0 = R1,               # basic reproduction number
    gamma = g1 # parameters of the Weibull distribution for the infection periods
  )
  out = run_DAMCMC(Y, N = iter, iota_dist = "exponential", rho = rho, theta_0 = theta1, par_prior = list(
    a_beta = 0.1, b_beta = 0.1,
    a_gamma = 0.1, b_gamma = 0.1,
    a_R0 = 0.1, b_R0 = 0.1,
    a_lambda = 0.1, b_lambda = 0.1
  ))
  return(out)
}


## "extract_param" function takes the output of "PDSIR_calc" and extract the parameters in the form of a matrix or calculates the ESS.
## result: output of "PDSIR_calc"
## chain_num: number of chains simulated. If you have more than one chain make a list of the outputs.
## iter: Number of iteration of the chain.
## warmup: number of warm-yp samples. 50% by default.
## init_sus: Initial number of susceptible.
## esscalc: Logical. If TRUE, returns the ESS. Otherwise, returns the matrix of parameters.

library(mcmcse)
extract_param = function(result, chain_num, iter, warmup = iter/2, init_sus, esscalc = FALSE)
{
  diff_parameter = list()
  ess1 = matrix(0, ncol = 3, nrow = chain_num)
  if(chain_num == 1)
  {
    result[[1]] = list(theta = result[[1]])
  }

  for(i in 1:chain_num)
  {
    diff_parameter[[i]] = matrix(0, ncol = 3, nrow = iter) ## nrow needs to be changed with iteration
    for(j in 1:iter)
    {
      diff_parameter[[i]][j, 1] =  result[[i]]$theta[[j]]$R0
      diff_parameter[[i]][j, 2] =  result[[i]]$theta[[j]]$gamma
      diff_parameter[[i]][j, 3] =  result[[i]]$theta[[j]]$beta*init_sus
    }
    diff_parameter[[i]] = tail(diff_parameter[[i]], (iter - warmup))
    colnames(diff_parameter[[i]]) = c("R0", "gamma", "beta")
    if(esscalc == TRUE)
    {
      ess1[i, ] = ess(diff_parameter[[i]])
    }
  }
  if(esscalc == TRUE)
  {
    return(ess1)
  }else{
    final = do.call( rbind, diff_parameter)
    return(final)
  }
}

## When we make inference for more than one data sets using PDSIR, the function is used to summarize the results.

## result_list: list of outputs of PDSIR_calc.
## init_sus: Initial number of susceptible.
## true_val: True values of the parameters. A vector of three dimension with true value of R_0, beta, gamma respectively.
## Tmax: Final observation time of epidemic.
## iter: Number of iteration of the chain.
## warmup: number of warm-yp samples. 50% by default.
## chain_num: Number of chains in each output

## returns means, sds, 95% coverage and ESS of the parameters in the order R_0, beta, gamma respectively.



PDSIR_summary = function(result_list, init_sus, true_val, Tmax = 10, iter = 1e5, warmup = iter/2, chain_num = 1)
{
  exPm = lapply(result_list, function(x) extract_param(x, chain_num, iter, warmup, init_sus) )
  mn = apply(matrix(unlist(   lapply(exPm, function(x) apply(x, 2, mean)) ), byrow = TRUE, ncol = 3), 2, mean)

  stddev = apply(matrix(unlist(   lapply(exPm, function(x) apply(x, 2, sd)) ), byrow = TRUE, ncol = 3), 2, mean)
  temp1 = matrix(unlist(   lapply(exPm, function(x) apply(x, 2, function(y) quantile(y, 0.025) )) ), byrow = TRUE, ncol = 3)
  temp2 = matrix(unlist(   lapply(exPm, function(x) apply(x, 2, function(y) quantile(y, 0.975) )) ), byrow = TRUE, ncol = 3)

  R0per = mean(temp1[, 1] < true_val[1] & temp2[, 1] > true_val[1])
  gamma_per = mean(temp1[, 2] < true_val[2] & temp2[, 2] > true_val[2])
  beta_per = mean(temp1[, 3] < true_val[3] & temp2[, 3] > true_val[3])
  esslist = apply(matrix(unlist(   lapply(exPm, ess )), byrow = TRUE, ncol = 3), 2, mean)

  return(list(means = mn, sds = stddev, coverage = c(R0per, gamma_per, beta_per), ESS = esslist))
}

simulated_data = list()
PDSIR_result = list()
for(i in 1:50)
{
  simulated_data[[i]] =  sellke(n = 1000, rho = 0.05, beta = 2, gamma = 0.5, Tmax = 10)
  PDSIR_result[[i]] = PDSIR_calc(data_Sellke =  simulated_data[[i]], rho = 0.08, Tmax = 10, iter = 100, b1 = 0.1, g1 = 0.1) ## rho is calculated so that we get around 20% acceptance probability.
}


temp = extract_param(result = PDSIR_result[[2]] , chain_num = 1, iter = 1e5, warmup = 5e4, init_sus = 1e3, esscalc = FALSE )
PDSIR_summary(result_list = PDSIR_result, init_sus = 1000, true_val = c(4, 0.5, 2))

