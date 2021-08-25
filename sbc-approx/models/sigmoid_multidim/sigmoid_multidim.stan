data {
  int <lower =0> n_priorpts;
  int  <lower =0> n_datapoints;
  matrix[n_priorpts, n_datapoints] y;
  vector[n_priorpts] x;
  real<lower=0> prior_width;
  real<lower=0> sigma;
}

parameters {
  real w;
  real b;
}

model {
  vector[n_datapoints] y_true = inv_logit(w * x + b);
  for (i in 1:n_priorpts){
    y[i,] ~ normal(y_true[i], sigma);
  }
  w ~ normal(0,prior_width);
  b ~ normal(0,prior_width);
}
