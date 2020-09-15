data{
  int N; 
  int B; // num of Boostrap
  vector[N] w;
}

transformed data {
  real mu_ = normal_rng(0, 1);
  real<lower = 0> sigma_ = lognormal_rng(0, 1);
  int K = 5;
  vector[N] y;
  matrix[K,B] BS;
  int I[K];
  for (n in 1:N)
    y[n] = normal_rng(mu_, sigma_);
  
  for(b in 1:B){
    for(k in 1:K)
      I[k] = categorical_rng(w);
      print(I);
      BS[,b] = y[I];
  }
}

generated quantities {
  matrix[K,B] BS0;
  BS0 = BS;
}
