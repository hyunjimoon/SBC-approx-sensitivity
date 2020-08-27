data{
  real hp;
}

transformed data {
  int<lower=0> N = 20;               
  int<lower =0> y[N];
  int n_par = 2;
  real<lower=0> theta_ = gamma_rng(hp, 1);
  for (n in 1:N)
    y[n] = neg_binomial_2_rng(theta_, 1);
}

parameters {
  real<lower = 0> theta;
}

model {
  theta ~ gamma(hp, 1);
  y ~ poisson(theta); 
}

generated quantities {
  int y_[N] = y;
  vector[n_par] pars_;
  pars_[1] = theta_;
  pars_[2] = y_[1];
  int ranks_[n_par] = {theta < theta_, y[1] < y_[1]}; 
}  