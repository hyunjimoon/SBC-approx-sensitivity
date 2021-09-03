data {
  // hyperparams
  int D;
  real theta_loc;
  real <lower = 0> theta_scale;
  // outcome
  vector[D] y;

}
parameters {
  real theta;
}

model {
  theta ~ normal(theta_loc, theta_scale);
  y ~ normal(theta, 1);
}
