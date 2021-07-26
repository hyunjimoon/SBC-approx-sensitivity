data {
  int<lower = 1> D;
  vector[D] X;
  vector[D] y;
  //real alpha_loc;
  //real alpha_scale;
  real beta_loc;
  real beta_scale;
}
parameters {
  real alpha;
  real beta;
}
model {
  alpha ~ normal(0, 1);
  beta ~ normal(beta_loc, beta_scale);
  y ~ normal(alpha + beta * X, 1.2);
}
