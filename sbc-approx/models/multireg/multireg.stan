data {
  int<lower = 1> N;
  int<lower = 1> P;
  matrix [N,P]X;
  vector[N] y;
  vector[P] param_loc;
  vector[P] param_scale;
}

parameters {
  vector[P] param;
}
// construct posterior with given prior and data model
model {
  param ~ normal(param_loc, param_scale);
  y ~ normal(X * param, 1.2);
}
// generate data with given prior
generated quantities {
  vector[P] param_;
  vector[N] y_;
  param_ = to_vector(normal_rng(param_loc, param_scale));
  y_ = to_vector(normal_rng(X * param, 1.2));
}
