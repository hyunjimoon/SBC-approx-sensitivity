library(rstan); library(cmdstanr); library(parallel); library(tidyverse); library(dplyr); library(posterior);library(bayesplot)
scriptDir <- getwd()
modelDir <- file.path(scriptDir, "models")
s_fit <- function(modelName, approxName, data, S){
  chains <- 1
  parallel_chains <- min(chains, detectCores())
  scriptDir <- getwd()
  delivDir <- file.path(scriptDir, "deliv", modelName, approxName)
  prefit <- file.path(delivDir, paste0(approxName, '-', S, '.rda'))
  stanfile <- file.path(modelDir, modelName, paste0(modelName, ".stan"))
  dir.create(delivDir)
  if (file.exists(prefit)){
    fit <- readRDS(prefit)
  }else{ 
    mod <- cmdstan_model(stanfile, quiet = FALSE)
    fit <- mod$sample(data, chains = chains, iter_warmup = 300, iter_sampling = 500, parallel_chains = parallel_chains, save_warmup = FALSE)
    dir.create(delivDir)
    fit$save_object(file = prefit)
  }
  fit
}

modelName <- "pois"
N <- 100000 # Number of replications
L <- 25     # Number of posterior samples
B <- 25     # Number of rank bins
shape = 8
rate = 2

#1. correct SBC
ranks <- rep(0, 100)
for (n in 1:N) {
  lambda <- rgamma(1, shape, rate)
  y <- rpois(1, lambda)
  post_shape <- shape + y 
  post_rate <- rate + 1
  post_lambda <- rgamma(L, post_shape, post_rate)
  ranks[n] <- sum(post_lambda < lambda)
}

pdf(file = file.path(scriptDir, "deliv", modelName,"plot",  paste0("correct", ".pdf")), width = 5, height = 5) 
sbc_hist <- hist(ranks,  breaks=seq(0, L + 1, 1) - 0.5,  plot=FALSE) 
plot(sbc_hist, main="SBC (Correct)", xlab="Prior Rank", yaxt='n', ylab="")
dev.off() 

# BS approx, rank compare
approxAlg = "BS1"
N_theta = 10
N_y = 25
N_data = 30 # ok for N_data !=1 ?
approxName  = paste0(approxAlg, "_", N_theta,"_", N_y, "_",N_data)
ranks <- rep(NA, N_theta)
post_lambda <- matrix(, nrow = N_y, ncol = N_theta)
for (n in 1:N_theta) {
  lambda <- rgamma(1, shape, rate)
  y <- rpois(N_data, lambda)
  for (m in 1:N_y){ # N_y = L
    y_BS <- sample(y, replace = T)
    sf_BS <- s_fit(modelName, approxName, data = list(N =N_data, y = as.array(y_BS)), S = n*N_y + m) 
    seq = floor(seq(1,500, length.out  = L))
    post_lambda_s <- as_draws_df(sf_BS$draws(variables ="lambda")) %>% slice(seq) %>% pull(lambda)
    post_lambda[m, n] <- mean(post_lambda_s)
  }
  ranks[n] <- sum(post_lambda[,n] < lambda)
}                                                                                                                                                                                   
sbc_hist <- hist(ranks,breaks=seq(0, L + 1, 1) - 0.5, plot=FALSE) 
title = paste0(approxName, " for theta, y, data sim#")
pdf(file = file.path(scriptDir, "deliv", modelName,"plot",  paste0(approxName, ".pdf")), width = 5, height = 5) 
plot(sbc_hist, main= title, xlab="Prior Rank", yaxt='n', ylab="")
dev.off() 

# BS approx, sample distribution compare
approxAlg = "BS2"
N_theta = 5
N_y = 50
N_data = 30
approxName  = paste0(approxAlg, "_", N_theta,"_", N_y, "_",N_data)
ranks <- rep(0, 100)
post_lambda_lst <- rep(NA, N_y)
for (n in 1:N_theta) {
  rank = 0
  lambda <- rgamma(1, shape, rate)
  y <- rpois(N_data, lambda)
  for (m in 1:N_y){ 
    y_BS <- sample(y, replace = T)
    sf_BS <- s_fit(modelName, approxName, data = list(N =N_data, y = as.array(y_BS)), S = m)  
    post_lambda <- as_draws_df(sf_BS$draws(variables ="lambda"))$lambda 
    post_lambda_lst[m] <- post_lambda
  }
}
sbc_hist <- hist(ranks,  plot=FALSE)
title = paste0(approxName, " for theta, y, data sim#")
pdf(file = file.path(scriptDir, "deliv", modelName,"plot", paste0(approxName, ".pdf")), width = 5, height = 5) 
hist(post_lambda_lst, main= title, xlab="Data avg. prior", yaxt='n', ylab="")
dev.off() 

# BS approx, sample distribution compare2
approxAlg = "BS3"
N_theta = 5
N_y = 50
N_data = 30
approxName  = paste0(approxAlg, "_", N_theta,"_", N_y, "_",N_data)
ranks <- rep(0, 100)
post_lambda_lst <- rep(NA, N_y)
for (n in 1:N_theta) {
  rank = 0
  lambda <- rgamma(1, shape, rate)
  y <- rpois(N_data, lambda)
  for (m in 1:N_y){ 
    y_BS <- sample(y, replace = T)
    sf_BS <- s_fit(modelName, approxName, data = list(N =N_data, y = as.array(y_BS)), S = m)  
    post_lambda <- as_draws_df(sf_BS$draws(variables ="lambda"))$lambda 
    post_lambda_lst[m] <- mean(post_lambda)
  }
}
sbc_hist <- hist(ranks,  plot=FALSE)
title = paste0(approxName, " for theta, y, data sim#")
pdf(file = file.path(scriptDir, "deliv", modelName,"plot", paste0(approxName, ".pdf")), width = 5, height = 5) 
hist(post_lambda_lst, main= title, xlab="Data avg. prior", yaxt='n', ylab="")
dev.off()


# IJ approx N_theta for loop
approxAlg = "IJ1.1"
N_theta = 10
N_data = 30
approxName  = paste0(approxAlg, "_", N_theta, "_",N_data)
ranks <- rep(NA, N_theta)
for (n in 1:N_theta) {
  lambda <- rgamma(1, shape, rate)
  y <- rpois(N_data, lambda) 
  y_BS <- sample(y, replace=TRUE, size=N_data)
  sf_IJ <- s_fit(modelName, approxName, data = list(N =N_data, y = y_BS), S = n) # "IJ var captures BS var better"
  loglik_draws <-  select(as_draws_df(sf_IJ$draws(variables = "log_lik")), contains("log_lik"))
  param_draws <- as_draws_df(sf_IJ$draws(variables ="lambda"))$lambda 
  infl_draws_mat <- N_data * cov(loglik_draws, param_draws) #infl_draws_mat scale is  much smaller comparable to lambda
  post_lambda <- infl_draws_mat[,1]
  ranks[n] <- sum(post_lambda< lambda)
}   
pdf(file = file.path(scriptDir, "deliv", modelName, "plot", paste0(approxName, ".pdf")), width = 5, height = 5) 
sbc_hist <- hist(c(ranks), plot=FALSE) # breaks=seq(0, L + 1, 1) - 0.5,
title = paste0(approxName, " for theta, data sim#")
plot(sbc_hist, main= title, xlab="Prior Rank", yaxt='n', ylab="")
dev.off() 

# IJ approx in N_theta and N_y for loop
approxAlg = "IJ1_2"
N_theta = 5
N_y = 25
N_data = 30
approxName  = paste0(approxAlg, "_", N_theta,"_", N_y, "_",N_data)
ranks <- rep(NA, N_theta)
for (n in 1:N_theta) {
  lambda <- rgamma(1, shape, rate)
  y <- rpois(N_data, lambda) 
  for (m in 1:N_y){
    y_BS <- sample(y, replace=TRUE, size=N_data)
    sf_IJ <- s_fit(modelName, approxName, data = list(N =N_data, y = y_BS), S = n*N_y + m) # "IJ var captures BS var"
    loglik_draws <-  select(as_draws_df(sf_IJ$draws(variables = "log_lik")), contains("log_lik"))
    param_draws <- as_draws_df(sf_IJ$draws(variables ="lambda"))$lambda
    infl_draws_mat <- N_data * cov(loglik_draws, param_draws) #infl_draws_mat scale is  much smaller comparable to lambda
    post_lambda <- infl_draws_mat[,1]
  }
}   
pdf(file = file.path(scriptDir, "deliv", modelName, "plot", paste0(approxName, ".pdf")), width = 5, height = 5) 
sbc_hist <- hist(c(ranks), plot=FALSE) # breaks=seq(0, L + 1, 1) - 0.5,
title = paste0(approxName, " for theta, y, data sim#")
plot(sbc_hist, main= title, xlab="Prior Rank", yaxt='n', ylab="")
dev.off() 

# IJ assume known mean, estimate var with ij_cov 
approxAlg = "IJ2"
N_theta = 30
N_y = 20
N_data = 30 # ok for N_data !=1 ?
approxName  = paste0(approxAlg, "_", N_theta,"_", N_y, "_",N_data)
ranks <- matrix(, nrow = N_y, ncol = N_theta)

for (n in 1:N_theta) {
  lambda <- rgamma(1, shape, rate)
  y <- rpois(N_data, lambda) # original was single sampling
  for (m in 1:N_y){
    y_BS <- sample(y, replace=TRUE, size=N_data)
    sf_IJ <- s_fit(modelName, approxName, data = list(N =N_data, y = y_BS), S = n*N_y + m) # IJ var captures BS var well.
    loglik_draws <- as_draws_df(sf_IJ$draws(variables = "log_lik"))%>%select(1:N_data)
    param_draws <- as_draws_df(sf_IJ$draws(variables ="lambda"))$lambda
    infl_draws_mat <- N_data * cov(loglik_draws, param_draws)
    ij_cov <- cov(infl_draws_mat, infl_draws_mat)
    #assume known mean 
    post_mean = post_shape/post_rate
    post_var = ij_cov
    post_shape_IJ = post_mean^2/post_var
    post_scale_IJ = post_mean/post_var
    post_lambda <- rgamma(L, post_shape_IJ, post_scale_IJ)
    ranks[m, n] <- sum(post_lambda < lambda)
  }
}
pdf(file = file.path(scriptDir, "deliv", modelName, "plot", paste0(approxName, ".pdf")), width = 5, height = 5) 
sbc_hist <- hist(c(ranks), plot=FALSE) # breaks=seq(0, L + 1, 1) - 0.5,
title = paste0(approxName, " for theta, y, data sim#")
plot(sbc_hist, main= title, xlab="Prior Rank", yaxt='n', ylab="")
dev.off() 

# incorrect SBC
ranks <- rep(0, 100)
for (n in 1:N) {
  lambda <- rgamma(1, shape, rate)
  y <- rpois(L, lambda)
  post_shape <- shape + y 
  post_rate <- rate + 1
  post_lambda <- rgamma(L, post_shape, post_rate)
  ranks[n] <- sum(post_lambda < lambda)
}
pdf(file = file.path(scriptDir, "deliv", modelName,"plot",  paste0("Incorrect", ".pdf")), width = 5, height = 5) 
sbc_hist <- hist(ranks,  breaks=seq(0, L + 1, 1) - 0.5,  plot=FALSE) 
title = paste0(approxName, " for theta, y, data sim#")
plot(sbc_hist, main="SBC (InCorrect)", xlab="Prior Rank", yaxt='n', ylab="")
dev.off() 