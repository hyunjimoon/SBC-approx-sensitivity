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


# IJ assume known mean, estimate var with ij_cov 
approxAlg = "IJ"

approxName  = paste0(approxAlg, "_", N,"_", M, "_",N_data)
ranks <- matrix(, nrow = M, ncol = N)

for (n in 1:N) {
  lambda <- rgamma(1, shape, rate)
  y <- rpois(N_data, lambda) 
  for (m in 1:M){
    y_BS <- sample(y, replace=TRUE, size=N_data)
    sf_IJ <- s_fit(modelName, approxName, data = list(N =N_data, y = y_BS), S = n*M + m) # IJ var captures BS var well.
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
