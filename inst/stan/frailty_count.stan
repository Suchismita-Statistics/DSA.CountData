functions{
  real[] two_ode(real t, real[] y, real[] parms, real[] x_r, int[] x_i){
      real b = parms[1];
      real g = parms[2];
      real r = parms[3];
      real n = parms[4];
      real dydt[2];

      dydt[1] = -b*y[2]*((y[1])^(1+n^2));
      dydt[2] = b*y[2]*((y[1])^(1+n^2)) - g*y[2];
      return dydt;
      }
}

data{
  int<lower = 0> N;
   int<lower = 0> K; // newly infected
  real<lower = 0> T_max;
  int<lower = 0> T_max_int;
  int<lower = 0> infection_count[T_max_int];
  real t0;
  real int_time[T_max_int];
}

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {
    real<lower = 0> b;
  real<lower = 0> g;
  real<lower = 0, upper = 0.1> r;
   real<lower = 0> n;
}

 transformed parameters{
   real<lower = 0> R0;
   R0 = b/g;
 }

model{
  real parms[4];
  real ic[2];
  real smax;
  real st_at_int[T_max_int + 1];
  real ode_sol[T_max_int, 2];



  parms[1] = b;
  parms[2] = g;
  parms[3] = r;
  parms[4] = n;

  ic[1] = 1.0;
  ic[2] = r;

  st_at_int[1] = 1.0;

  ode_sol = integrate_ode_rk45(two_ode, ic, t0, int_time, parms, x_r, x_i);
  st_at_int[2:(T_max_int + 1)] = ode_sol[1:T_max_int, 1];
  smax = ode_sol[T_max_int, 1];

  for(i in 1:T_max_int){
    target += infection_count[i]*log(st_at_int[i] - st_at_int[i+1]);
  }
  target += (N-K)*log(smax) + gamma_lpdf(b | 1, 1) + gamma_lpdf(g | 1, 1) + beta_lpdf(r | 1, 1) + gamma_lpdf(n | 1, 1);
}

