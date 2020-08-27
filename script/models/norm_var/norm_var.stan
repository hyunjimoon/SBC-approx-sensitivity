data{
  real<lower = 0> hp;
}

transformed data {
  int n_par = 2;                                                        
  int<lower = 0> N = 10;
  real theta_ = normal_rng(1, hp);
  vector[N] y;
  for (n in 1:N)
    y[n] = normal_rng(theta_,2);
}

parameters {
  real theta;
}

model {
  theta ~ normal(1, hp);
  y ~ normal(theta, 1);
}

generated quantities {
  vector[N] y_ = y;
  vector[n_par] pars_;
  pars_[1] = theta;
  pars_[2] = y_[1];
  int ranks_[n_par] = {theta <theta_, y[1] < y_[1]};  
}  