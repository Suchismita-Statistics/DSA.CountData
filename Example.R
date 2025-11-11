## Analysis using simulated data
library(DSA.CountData)
simulated_data = sellke(
  n = 1000,
  rho = 0.05,
  beta = 2,
  gamma = 0.5,
  Tmax = 10
)
head(simulated_data)

T.max_analysis = max(simulated_data[, 1])
initial_sus = simulated_data[simulated_data[, 1] != 0, ]
N = dim(initial_sus)[1]
initial_sus = initial_sus[order(initial_sus[, 1]), ]
infect_during_ep = subset(initial_sus, initial_sus[, 1] < T.max_analysis)
t = infect_during_ep[, 1]


#### We make the time point vector all integers between one to the final time and will adjust the vector of infection count so that,
## if we do not observe infections in a day, if will take value zero there.
table_t = table(floor(t))
infection_days = as.numeric(names(table_t))
infection_count = as.vector(table_t)
inf_time_final = numeric(length = T.max_analysis)

inf_time_final[infection_days + 1] = infection_count

out = new_data_ct_lkd(
  time_points = 1:T.max_analysis,
  infection_count = inf_time_final,
  Final_time = T.max_analysis
)


### Likewise new_data_ct_lkd function can be used to analyze a new count data.



## Analysis via count likelihood of the SIR Model

count_result = count_lkd(simulated_data,
                         T.max_analysis = 10,
                         frailty = FALSE)

## Analysis via count likelihood of the Frailty Model

frailty_result = count_lkd(simulated_data,
                           T.max_analysis = 10,
                           frailty = TRUE)


## One example of the simulation study without frailty. Make frailty =.TRUE for using frailty model.
count_result = list()


full_result = list()
uniform_result = list()

data_sellke = list()
for (i in 1:5)
{
  data_sellke[[i]] =  sellke(
    n = 1000,
    rho = 0.05,
    beta = 2,
    gamma = 0.5,
    Tmax = 10
  )

  ## Do inference using count likelihood when only daily infection counts are available
  count_result[[i]] = count_lkd(data_sellke[[i]],
                                T.max_analysis = 10,
                                frailty = FALSE)

  ## Do inference when only daily infection counts are available, assuming infection times are uniformly distributed throughout the day
  uniform_result[[i]] = unif_approx(data_sellke[[i]], T.max_analysis = 10)

  ## Do inference using complete likelihood when exact infection and recovery times are available
  full_result[[i]] = full_lkd(data_sellke[[i]], T.max_analysis = 10)
}



######
summary_stan(result_list = count_result,
             index = 1,
             true_value = 2)
summary_stan_list(result_list = count_result, true_value = c(2, 0.5, 0.05))
