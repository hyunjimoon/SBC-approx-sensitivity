source("tools/functions.r")
##' Auto calibrate the initial prior samples using SBC iteration
##'
##' @param priors rvars<n_priorval>[n_dataval] prior values i.e. parameter values to be tested "draws_rvars"
##' @param predictor rvars<n_dataval>[n_priorval] predictor values
##' @param generator_priorvals function that generates datasets given each value in `priors`
##' @param stan_mod stan_model that samples posterior with `datasets` simulated from the former two
##' @param target_vars function of parameters with which SBC iteration convergence are judged
##' @param n_dataval the number of simulated data points from each prior value (default 4000)
##' @param n_sample the number of simulated posterior from each prior value (default 4000)
##' @param cnt function of parameters with which SBC iteration convergence are judged
##' @param evolve_df datafraame holding `median, mad (or sd)` of samples from every iteration
##' @param delivDir save location of output result
##' @return  next priors summarized from `n_priorval` * `n_sample` posterior samples
##' @export

selfCalib <- function(priors, predictor, generator_priorvals, stan_mod, target_vars, n_dataval, n_sample, cnt, evolve_df, delivDir){
  n_priorval <- niterations(priors)
  backend <- SBC_backend_cmdstan_sample(stan_mod, iter_sampling = n_sample, chains = 1) # M/4, chains = 4)
  results <- compute_results(generator_priorvals(priors, n_dataval, predictor), backend)
  for (tv in target_vars){
    summ <- summarise_draws(priors, median, sd) %>% filter(variable == tv)
    evolve_df[[tv]]$median[cnt] <- as.numeric(summ["median"])
    evolve_df[[tv]]$sd[cnt] <- as.numeric(summ["sd"])
  }

  next_priors <- adj_post(priors, results, target_vars, sumtype = "adj_postor_reweight", cnt)

  priors_next_priors_close <- iter_stop(priors, next_priors, results, target_vars = NULL, sumtype = "adj_postor_reweight")
  # iter_stop much stabilized when n_priorval vs n_priorval compared to n_priorval vs n_priorval * 4000 (= n_sample)
  if (priors_next_priors_close || is.na(priors_next_priors_close)){ #NA if the two are the same `draws_rvars`
    for (tv in target_vars){
      summ <- summarise_draws(priors, median, sd) %>% filter(variable == tv)
      evolve_df[[tv]]$median[cnt] <- as.numeric(summ["median"])
      evolve_df[[tv]]$sd[cnt] <- as.numeric(summ["sd"])
      csv_save(priors, delivDir, cnt, type = "each")
      csv_save(evolve_df, delivDir, cnt, type = "each")
      #intv_plot_save(evolve_df[[v]])
    }

    return (priors) # calibrated only for the target
  }else{
    cnt = cnt + 1
    return (selfCalib(next_priors, predictor, generator_priorvals, stan_mod, target_vars, n_dataval, n_sample, cnt, evolve_df, delivDir))
  }
}
##' Judge whether the SBC iteration have converged
##'
##' @param priors numeric vector of prior values i.e. parameter values to be tested
##' @param posteriors numeric vector of posterior values i.e. sampled parameter values
##' @param results computed results with `generator_truepoints(priors), backend)
##' @param target_vars function of parameters with which SBC iteration convergence are judged
##' @param n_sample the number of posterior samples for each prior value (default 4000)
##' @param cnt function of parameters with which SBC iteration convergence are judged
##' @param sumtype types are classified by object and size of compared distribution
##'                object: generic parameters vs targeted functions of parameter
##'                size: `niterations` are same (`prior2prior`) vs different (`prior2post`)
##' @param bins the number of bins to discretize samples
##' @return distance between prior and posterior samples
##' @export
iter_stop <- function(priors, next_priors, results, target_vars, sumtype, bins = 20){
  if (is.null(target_vars)){
    post_r_loc <- lapply(next_priors, mean)
    post_r_scale <- lapply(next_priors, sd)
    r_loc <- list()
    r_scale <- list()
    for (par in names(priors)){
      r_loc <- append(r_loc, E(priors[[par]]) / post_r_loc[[par]])
      r_scale <- append(r_scale, sd(priors[[par]]) / post_r_scale[[par]])
    }
    return (all(r_loc > 0.9 && r_loc < 1.1 && r_scale > 0.9 && r_scale < 1.1 ))
  }else{
    return(all(unlist(lapply(target_vars, FUN = function(tv) {cjs_dist(draws_of(priors[[tv]]), draws_of(next_priors[[tv]])) < 0.1}))))
  }
}

# Resample for valid compare between `priors`, `next_priors` and to control `n_priorval`
# n_priorval * n_sample posterior  n_priorvals as comparison threshold is possible for the same number of samples
##'
##' @param priors numeric vector of prior values i.e. parameter values to be tested
##' @param results computed results with `generator_truepoints(priors), backend)
##' @param target_vars function of parameters with which SBC iteration convergence are judged
##' @param sumtype types are classified by object and size of compared distribution
##'                object: generic parameters vs targeted functions of parameter
##'                size: `niterations` are same (`prior2prior`) vs different (`prior2post`)
##' @param cnt needed to keep track of iteration counts
##' @return resampled posterior with prior information
##' @export
adj_post <-function(priors, results, target_vars, sumtype, cnt = 0){
  n_priorval <- niterations(priors)
  post_mtr <- SBC_fit_to_draws_matrix(results$fits[[1]])
  n_post <- nchains(post_mtr) * niterations(post_mtr) # default nchains = 4
  priors_tpl <- priors # template
  for (tv in target_vars){
    post_mtr <- matrix(NA, nrow = n_post, ncol = n_priorval)
    for (i in 1:n_priorval){
      post_mtr[,i] <- c(subset_draws(SBC_fit_to_draws_matrix(results$fits[[i]]), variable = tv))
    }
    priors_v <- c(as_draws_df(priors)[[tv]])
    pp_overlay_save(priors_v, post_mtr, tv, cnt)

    if(n_priorval < 10 || sumtype == "MtoM"){
      priors_tpl[[tv]] <- rvar(c(post_mtr))
    }else if(grepl("Mto1", sumtype, fixed = TRUE)){
      if (sumtype == "Mto1_randpick"){
        for (i in 1:n_priorval){
          post_mtr_i <- subset_draws(SBC_fit_to_draws_matrix(results$fits[[i]]), varaible = tv)
          priors_tpl[[tv]] <- rvar(post_mtr_i)[sample(1:n_post, n_priorval)]
        }
      }else if (sumtype == "Mto1_reweight"){
        prior_sort  <- sort(results$stats %>% filter(parameter == tv) %>% pull(simulated_value))
        ar <- resample_draws(as_draws_rvars(sort(priors[[tv]])), tabulate(ecdf(prior_sort)(c(post_mtr)) * n_priorval, nbins = n_priorval))
        priors_tpl[[tv]] <- ar[[1]]
      }
    }
  }
  return (priors_tpl)
}
