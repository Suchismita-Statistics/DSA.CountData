functions{
  real[] two_ode(real t, real[] y, real[] parms, real[] x_r, int[] x_i){
      real b = parms[1];
      real g = parms[2];
      real r = parms[3];
      real dydt[2];

      dydt[1] = -b*y[1]*y[2];
      dydt[2] = b*y[1]*y[2] - g*y[2];
      return dydt;
      }
}

data{
  int<lower = 0> inf_len;
  int<lower = 0> infection_count[inf_len];
  int<lower = 0> K;
  real t0;
  real int_time[inf_len];
}

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {
   real<lower=0, upper = 50> b;
   real<lower = 0, upper = b> g; //Assuming it is a big Epidemic.
  real<lower=0, upper = 1 > r;
}

transformed parameters{
  real<lower = 0> R0;
  R0 = b/g;
}

model{
  real parms[3];
  real ic[2];
  real smax;
  real ode_sol[inf_len+1, 2];


  parms[1] = b;
  parms[2] = g;
  parms[3] = r;

  ic[1] = 1.0;
  ic[2] = r;

  ode_sol[1, ] = ic;
  ode_sol[2:(inf_len+1), ] = integrate_ode_rk45(two_ode, ic, t0, int_time, parms, x_r, x_i, 1e-8, 1e-8, 1e6);
  smax = ode_sol[inf_len+1, 1];

  for(i in 1:(inf_len)){
   target += infection_count[i]*(log(ode_sol[i, 1] - ode_sol[i+1, 1]) - log(1-smax));
  }
 target +=  gamma_lpdf(b|0.1, 0.1) + gamma_lpdf(g|0.1, 0.1) + uniform_lpdf(r|0, 1);
}

generated quantities{
  real<lower = 0> N_est;
  real parms[3];
  real ic[2];
  real smax;
  real ode_sol[inf_len+1, 2];


  parms[1] = b;
  parms[2] = g;
  parms[3] = r;


  ic[1] = 1.0;
  ic[2] = r;

  ode_sol[1, ] = ic;
  ode_sol[2:(inf_len+1), ] = integrate_ode_rk45(two_ode, ic, t0, int_time, parms, x_r, x_i, 1e-8, 1e-8, 1e6);
  smax = ode_sol[inf_len+1, 1];

  N_est = K/(1-smax);
}
