# DSA.CountData
To install:

```
install.packages("devtools")
devtools::install_github("Suchismita-Statistics/DSA.CountData")
```

# About the package

This package aims to apply DSA models for count data. The SIR Model or frailty model can be used to infer about the parameters -  infection rate ($\beta$) and recovery rate ($\gamma$). It further estimates the limiting proportion of susceptible and infected ($\rho$). The DSA method does not necessarily require information on the number of initial susceptible individuals. It can also be estimated.  

This package contains all the R codes required to implement the paper.  
new_data_ct_lkd() is the main function to apply count likelihood for a new daily count data, which simulates HMC posterior samples via stan. We used rstan::version 2.21.8.  
Please look at the Example.R file to see the usage of the important functions.  

```
git clone https://github.com/Suchismita-Statistics/DSA.CountData.git
cd DSA.CountData
Rscript Example.R
```

To estimate the infection and recovery rates of a new dataset, use the following code. The code illustrated using simulated data, but it is easy to change to a new one. Just change the lines mentioned in steps 1 and 2.

```
library(DSA.CountData)
simulated_data = sellke(
  n = 1000,
  rho = 0.05,
  beta = 2,
  gamma = 0.5,
  Tmax = 10
)
head(simulated_data)

T.max_analysis = max(simulated_data[, 1]) ## Step 1: Exchange with final observation time or the maximum time.
initial_sus = simulated_data[simulated_data[, 1] != 0, ]
N = dim(initial_sus)[1]
initial_sus = initial_sus[order(initial_sus[, 1]), ]
infect_during_ep = subset(initial_sus, initial_sus[, 1] < T.max_analysis)
t = infect_during_ep[, 1]


#### Important: We make the time point vector all integers between one to the final time and will adjust the vector of infection count so that,
## if we do not observe infections in a day, if will take value zero there.
table_t = table(floor(t))
infection_days = as.numeric(names(table_t))
infection_count = as.vector(table_t)
inf_time_final = numeric(length = T.max_analysis)

inf_time_final[infection_days + 1] = infection_count  ## Step 2: Give a vector containing the infection counts per day (or any equal interval).


out = new_data_ct_lkd(
  time_points = 1:T.max_analysis,
  infection_count = inf_time_final,
  Final_time = T.max_analysis
)


```

