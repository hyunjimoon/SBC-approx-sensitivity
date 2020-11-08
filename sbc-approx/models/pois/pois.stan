data{
  int N; 
  int y[N];
}
parameters{
  real<lower = 0> lambda;
}
model{
  y ~ poisson(lambda);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = poisson_lpmf(y[n] |lambda);
  }
}