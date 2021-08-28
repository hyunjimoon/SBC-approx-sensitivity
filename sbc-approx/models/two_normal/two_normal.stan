data {
  int N;
  vector[N] y;
  real<lower=0> prior_width;
}

parameters {
  real loc;
  real <lower = 0> scale;
}

model {
  loc ~ normal(0, prior_width);
  scale ~ exponential(prior_width);
  y ~ normal(loc, scale);
}
