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
  int<lower = 0> N; // initial susceptible
  int<lower = 0> K; // newly infected
  int<lower = 0> L; // recovered among K
  int<lower = 0> M; // initial infected


  real<lower = 0> T_max;
  real<lower = 0> t0;
  real<lower = 0> infection_times[K+1];
  real<lower = 0> w[K];
  int<lower = 0> L_tilde;
  real<lower = 0> epsilon[L_tilde];
}

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {
  real<lower = 0> b;
  real<lower = 0> g;
  real<lower = 0, upper = 1.0> r;
}
transformed parameters{
  real<lower = 0> R0;
  R0 = b/g;
}

model{
  real parms[3];
  real ic[2];
  real s[K+1, 2];
  real smax;

  parms[1] = b;
  parms[2] = g;
  parms[3] = r;

  ic[1] = 1.0;
  ic[2] = r;

  s = integrate_ode_rk45(two_ode, ic, t0, infection_times, parms, x_r, x_i, 1e-8, 1e-8, 1e6);

  smax = s[K+1, 1];

  for(i in 1:K){
    target += log(s[i, 1]) + log(s[i, 2]) - g*w[i] + log(b);
  }
  for(i in 1:L_tilde)
  {
    target += log(g) -  g*epsilon[i];
  }
  target += +(N - K)*log(smax) + L*log(g) - g*(M - L_tilde)*T_max + gamma_lpdf(b | 0.1, 0.1) + gamma_lpdf(g | 0.1, 0.1) + gamma_lpdf(r | 0.1, 0.1);
}

