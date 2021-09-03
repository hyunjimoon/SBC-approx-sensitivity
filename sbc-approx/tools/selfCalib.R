source("tools/functions.r")
##' Auto calibrate the initial prior samples using SBC iteration
##'
##' @param param rvars<S>[N] prior values i.e. parameter values to be tested "draws_rvars"
##' @param generator function that generates datasets given each value in `param`
##' @param predictor rvars<N>[S] predictor values
##' @param stan_mod stan_model that samples posterior with `datasets` simulated from the former two
##' @param target_vars function of parameters with which SBC iteration convergence are judged
##' @param N the number of simulated data points from each prior value (default 4000)
##' @param n_sample the number of simulated posterior from each prior value (default 4000)
##' @param cnt function of parameters with which SBC iteration convergence are judged
##' @param evolve_df datafraame holding `median, mad (or sd)` of samples from every iteration
##' @param delivDir save location of output result
##' @return  next param summarized from `S` * `n_sample` posterior samples
##' @export

selfCalib <- function(generator, hyperparam, param, predictor, backend, target_vars, cnt, thin, evolve_df, delivDir){
  result <- compute_results(generator(hyperparam, param, predictor), backend, thin = thin)
  for (tv in target_vars){
    summ <- summarise_draws(param, median, sd) %>% filter(variable == tv)
    evolve_df[[tv]]$median[cnt] <- as.numeric(summ["median"])
    evolve_df[[tv]]$sd[cnt] <- as.numeric(summ["sd"])
  }
  next_param <- update_param(param, result, target_vars, cnt, delivDir, thin = thin)
  param_next_param_close <- iter_stop(param, next_param, result, target_vars = NULL)
  # iter_stop much stabilized when S vs S compared to S vs S * 4000 (= n_sample)
  if (param_next_param_close || is.na(param_next_param_close)){ #NA if the two are the same `draws_rvars`
    for (tv in target_vars){
      summ <- summarise_draws(param, median, sd) %>% filter(variable == tv)
      evolve_df[[tv]]$median[cnt] <- as.numeric(summ["median"])
      evolve_df[[tv]]$sd[cnt] <- as.numeric(summ["sd"])
      csv_save(param, delivDir, cnt, type = "each")
      csv_save(evolve_df, delivDir, cnt, type = "each")
      #intv_plot_save(evolve_df[[v]])
    }
    return (param) # calibrated only for the target
  }else{
    cnt = cnt + 1
    return (selfCalib(generator, hyperparam, next_param, predictor, backend, target_vars, cnt, evolve_df, delivDir))
  }
}
##' Judge whether the SBC iteration have converged
##'
##' @param param numeric vector of prior values i.e. parameter values to be tested
##' @param posteriors numeric vector of posterior values i.e. sampled parameter values
##' @param result computed result with `generator_truepoints(param), backend)
##' @param target_vars function of parameters with which SBC iteration convergence are judged
##' @param n_sample the number of posterior samples for each prior value (default 4000)
##' @param cnt function of parameters with which SBC iteration convergence are judged
##' @param sumtype types are classified by object and size of compared distribution
##'                object: generic parameters vs targeted functions of parameter
##'                size: `niterations` are same (`prior2prior`) vs different (`prior2post`)
##' @param bins the number of bins to discretize samples
##' @return distance between prior and posterior samples
##' @export
iter_stop <- function(param, next_param, target_vars, bins = 20){
  if (is.null(target_vars)){
    post_r_loc <- lapply(next_param, mean)
    post_r_scale <- lapply(next_param, sd)
    r_loc <- list()
    r_scale <- list()
    for (par in names(param)){
      r_loc <- append(r_loc, E(param[[par]]) / post_r_loc[[par]])
      r_scale <- append(r_scale, sd(param[[par]]) / post_r_scale[[par]])
    }
    return (all(r_loc > 0.9 && r_loc < 1.1 && r_scale > 0.9 && r_scale < 1.1 ))
  }else{
    return(all(unlist(lapply(target_vars, FUN = function(tv) {cjs_dist(draws_of(param[[tv]]), draws_of(next_param[[tv]])) < 0.1}))))
  }
}

# Resample for valid compare between `param`, `next_param` and to control `S`
# S * n_sample posterior  Ss as comparison threshold is possible for the same number of samples
##'
##' @param param numeric vector of prior values i.e. parameter values to be tested
##' @param result computed result with `generator_truepoints(param), backend)
##' @param target_vars function of parameters with which SBC iteration convergence are judged
##' @param sumtype types are classified by object and size of compared distribution
##'                object: generic parameters vs targeted functions of parameter
##'                size: `niterations` are same (`prior2prior`) vs different (`prior2post`)
##' @param cnt needed to keep track of iteration counts
##' @return resampled posterior with prior information
##' @export
update_param <-function(param, result, target_vars, cnt = 0, delivDir, thin = 10){
  S <- ndraws(param[[1]])
  M <- dim(SBC_fit_to_draws_matrix(result$fits[[1]]))[1] / thin
  next_param <- param # template
  #next_param_mtr <- matrix(NA, nrow = S, ncol = S) #summarize n_sample to S
  g <- list()
  for (tv in target_vars){
    param_ord <- sort(c(as_draws_df(param)[[tv]]))
    post_mtr <- matrix(NA, nrow = S, ncol = M)
    for (i in 1:S){
      fit_mtr <- subset_draws(SBC_fit_to_draws_matrix(result$fits[[i]]), variable = tv) # change for tv, i order?
      fit_thinned <- posterior::thin_draws(fit_mtr, thin)
      post_mtr[i,] <- c(fit_thinned)
    }
    next_param_vec <- draws_of(resample_draws(as_draws(rvar(param_ord)), tabulate(ecdf(param_ord)(c(post_mtr)) * S, nbins = S))[[1]])
    g[[tv]] <- ppc_hist(param_ord, matrix(next_param_vec, ncol = length(param_ord)))
    next_param[[tv]] <- rvar(next_param_vec)
  }
  ggarrange(g[[target_vars[1]]], g[[target_vars[2]]], nrow = 2)
  ggsave(file =  file.path(delivDir, paste0(paste0(paste0(paste(target_vars[1], target_vars[2]), cnt, "_"), "pp.png", sep = ""))),  bg = "white")
  return (next_param)
}


