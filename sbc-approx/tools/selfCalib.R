source("tools/functions.r")
selfCalib <- function(stan_model, prior, pars, data, N, M, cnt, evolve_df, delivDir, is_param = NULL){
  data$theta_loc = mean(prior[[pars]])  #mean(prior) #summarise_draws(workflow$prior_samples)$mean
  data$theta_scale = sd(prior[[pars]]) #summarise_draws(workflow$prior_samples)$sd
  if(is_param){
    generator <- function(){
      function(theta_hp){
        theta <- rnorm(1, theta_hp)
        list(
          generated = rnorm(8, theta, 1),
          parameters = list(
            theta = theta
          )
        )
      }
    }
    theta_hp <- unlist(summarise_draws()[c("mean", "sd")])
    workflow <- SBCWorkflow$new(stan_model, generator())
    workflow$simulate(n_sbc_iterations = N, theta_hp = theta_hp)
    workflow$fit_model(sample_iterations = M, warmup_iterations = M, data)
    next_prior <- workflow$posterior_samples[pars]
  }else{
    # nonparameteric generator but parameteric stan prior
    generator_np <- function(){
      function(theta){
        list(
          generated = rnorm(8, theta, 1),
          parameters = list(
            theta = theta
          )
        )
      }
    }
    workflow <- SBCWorkflow$new(stan_model, generator_np())
    workflow$simulate(n_sbc_iterations = N, param = is_param, as_draws_df(prior)[[pars]]) #custome_prior
    workflow$fit_model(sample_iterations = M, warmup_iterations = M, next_data)
    next_prior <- post_summ(workflow, pars, sumtype = "filtering")
  }
  t_prior <- workflow$prior_samples[pars]
  t_post <- workflow$posterior_samples[pars]
  overlay_rvar_pri_post(t_prior, t_post, pars, cnt)
  if(cnt == 1){
    evolve_df[row(evolve_df)==cnt] <-summarise_draws(t_prior, median, mad)[2:3]
  }else{
    d <- summarise_draws(t_prior, median, mad)[2:3]
    evolve_df <-rbind(evolve_df, as.numeric(d))
  }
  if (iter_stop(t_prior, t_post)){
    csv_store(t_prior, delivDir, cnt)
    csv_store(evolve_df, delivDir, cnt,  type = "evolve")
    intv_plot_save(evolve_df)
    return (prior) # calibrated only for the target
  }
  else{
    csv_store(t_prior, delivDir, cnt)
    cnt = cnt + 1
    return (selfCalib(stan_model, next_prior, pars, data, N, M, cnt, evolve_df, delivDir, is_param = is_param))
  }
}

iter_stop <- function(prior, post, pars){
  post_r_loc <- lapply(post, mean)
  post_r_scale <- lapply(post, sd)
  r_loc <- list()
  r_scale <- list()
  for (par in names(prior)){
    r_loc <- append(r_loc, E(prior[[par]]) / post_r_loc[[par]])
    r_scale <- append(r_scale, sd(prior[[par]]) / post_r_scale[[par]])
  }
  #print(paste0(paste0("r_loc ", r_loc, paste0(" r_scale ", r_scale))))
  # NMP_G  %>% group_by(par, prior) %>% summarise("post_median" = median(value)) %>%
  #   ggplot(aes(x = prior, y = post_median, colour = par))  + geom_point() + xlim(min(prior_theta), max(prior_theta)) + ylim(min(prior_theta), max(prior_theta)) + coord_equal()
  # ggsave(file = file.path(delivDir, paste0(paste0(paste0(modelName, "_"), ".png"))), width = 5, height = 5)
  return (all(r_loc > 0.9 && r_loc < 1.1 && r_scale > 0.9 && r_scale < 1.1 ))
}

# summarize NM posterior samples to N for each parameter
post_summ <- function(workflow, pars, sumtype){
  prior_rv <- workflow$prior_samples
  prior <- as_draws_df(workflow$prior_samples[[pars]])$x
  post <- as_draws_df(workflow$posterior_samples[[pars]])$x
  if (sumtype == "filtering"){
    if(length(prior) < 10){
      return (workflow$posterior_samples[names(prior_rv)])
    }else{
      q_prior <-as_draws_df(workflow$prior_samples[[pars]]) # to borrow the frame
      q_prior$x  <- sort(prior)
      post <- as_draws_df(workflow$posterior_samples[[pars]])$x
      prior_rv[pars] <- resample_draws(as_draws_df(workflow$prior_samples), tabulate(ecdf(prior)(post) * N, nbins = N))[[par]]
      return (prior_rv)
    }
  } else if(sumtype == "sample"){ # only choose the first
    return (subset_draws(post,variable = names(prior), iteration = 1))
  } else if(sumtype == "median"){ # median of M from each N prior
  }
}

initDf <-function(L, summary, pars = NA){
  if(summary == "ms"){
    df <- data.frame(
                     mean = rep(NA,L),
                     sd = rep(NA,L)
    )
  }else if (summary == "q"){
    df <- data.frame(
                     q1 = rep(NA,L),
                     q3 = rep(NA,L)
    )
  }else if (summary == "pars"){
    df <- data.frame(
      median = rep(NA,L),
      mad = rep(NA,L)
      #tau_mean = rep(NA,L),
      #tau_sd = rep(NA,L)
    )
  }
  df
}
