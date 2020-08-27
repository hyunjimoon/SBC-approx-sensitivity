data{
  real<lower = 0> hp_sd;
}

transformed data {
  real mu_ = normal_rng(0, 1);
  real<lower = 0> sigma_ = lognormal_rng(0, hp_sd);
  int n_par = 2;                                                                      
  int<lower = 0> N = 10;
  vector[N] y;
  for (n in 1:N)
    y[n] = student_t_rng(4, mu_, sigma_);
}

parameters {
  real mu;
  real<lower = 0> sigma;
}

model {
  mu ~ normal(0, 1);
  sigma ~ lognormal(0, hp_sd);
  y ~ normal(mu, sigma);
}

generated quantities {
  vector[N] y_ = y;
  vector[n_par] pars_;
  pars_[1] = mu_;
  pars_[2] = sigma_;
  int ranks_[n_par] = {mu < mu_, sigma <sigma_};  
}  