data {
  int N;
  vector[N] y;
  vector[N] x;
  real<lower=0> prior_width;
  real<lower=0> sigma;
}

parameters {
  real w;
  real b;
}

model {
  vector[N] y_true = inv_logit(w * x + b);
  y ~ normal(y_true, sigma);
  w ~ normal(0,prior_width);
  b ~ normal(0,prior_width);
}
