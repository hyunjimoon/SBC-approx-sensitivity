data{
  real<lower = 0> hp_nu;
}

transformed data {
  real<lower = 0> sigma_ = lognormal_rng(0, 1);
  real<lower = 0> nu_ = lognormal_rng(hp_nu, 1);
  int n_par = 2;                                                                  
  int<lower = 0> N = 10;
  vector[N] y;
  for (n in 1:N)
    y[n] = student_t_rng(nu_, 0, sigma_);
}

parameters {
  real<lower = 0> sigma;
}

model {
  sigma ~ lognormal(0, 1);
  y ~ normal(0, sigma);
}

generated quantities {
  vector[N] y_ = y;
  vector[n_par] pars_;
  pars_[1] = sigma_;
  pars_[2] = y_[1];
  int ranks_[n_par] = {sigma <sigma_, y[1] < y_[1]};  
}  