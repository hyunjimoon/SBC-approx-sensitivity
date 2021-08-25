data {
  int D;
  vector[D] y;
  real theta_loc;
  real theta_scale;
}
parameters {
  real theta;
}

model {
  theta ~ normal(theta_loc, theta_scale);
  y ~ normal(theta, 1);
}