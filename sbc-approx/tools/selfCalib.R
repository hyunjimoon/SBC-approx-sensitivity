selfCalib <- function(workflow, prior, data, pars, N, M, cnt, evolve_df, delivDir){ # replace data with its placeholder
  # targeting params 
  t_prior <- prior[pars]
  workflow$simulate(n_sbc_iterations = N, custom_prior = prior)
  #clamped stan could be applied for inference
  D <- dim(as_draws_array(workflow$simulated_y))[1]
  y <- workflow$simulated_y
  prior_summ <- summarise_draws(workflow$prior_samples)
  theta_hp <- prior_summ$mean
  data = list("D"=D, "y"=y, "theta_hp"=theta_hp)
    
  workflow$fit_model(sample_iterations = M, warmup_iterations = M, data) #better way: y_sim <-  workflow$simulated_y 
  t_post <- workflow$posterior_samples[pars]
  if(cnt == 1){
    evolve_df[row(evolve_df)==cnt] <-c(cnt, summarise_draws(t_prior, median, mad)[2:3])
  }else{
    d <- summarise_draws(t_prior, median, mad)
    d$variable <- cnt
    evolve_df <-rbind(evolve_df, as.numeric(d))
  }
  if (iter_stop(t_prior, t_post)){
    write.csv(as_draws_df(t_prior), file =  file.path(delivDir, paste0(paste0(cnt, "_"), "final_scp.csv", sep = "")))
    write.csv(apply(evolve_df,2,as.character), file =  file.path(delivDir, paste0(paste0(cnt, "_"), "evolve.csv", sep = "")))
    intv <- subset_draws(mutate_variables(as_draws_df(lapply(evolve_df, as.numeric)), low1sd = (median + mad), up1sd = (median - mad)) , c("iter", "low1sd", "up1sd"))
    intv <- reshape2::melt(intv, id.vars = "iter") 
    intv <- filter(intv, variable == "low1sd" | variable == "up1sd")
    ggplot(intv, aes(x = as.numeric(iter), y = value,  color = variable) ) +
      geom_line() + ggtitle(sprintf("target par: %s, N: %s, M: %s ", pars, N, M))
    ggsave(file = file.path(delivDir, paste0(paste0(paste0(modelName, cnt), "_"), "evolove.png")), width = 5, height = 5)
    return (prior, evolve_df) # calibrated only for the target
  }
  else{
    write.csv(as_draws_df(t_prior), file =  file.path(delivDir, paste0(paste0(cnt, "_"), "scp.csv", sep = "")))
    cnt = cnt + 1
    #return (selfCalib(workflow, post_summ(workflow$prior_samples, workflow$posterior_samples, pars, sumtype = "sample"), data, pars, N, M, cnt, evolve_df, delivDir))
    return (selfCalib(workflow, post_summ_wf(workflow, pars, sumtype = "sample"), data, pars, N, M, cnt, evolve_df, delivDir))
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
  return (all(r_loc > 0.7 && r_loc < 13 && r_scale > 0.7 && r_scale < 1.3 ))
}

# summarize NM posterior samples to N for each parameter
post_summ_tmp <- function(prior, post, pars, sumtype){
  if (sumtype == "filtering"){ 
    prior <- as_draws_df(workflow$prior_samples[[pars]])$x
    q_prior <-as_draws_df(workflow$prior_samples[[pars]]) # to borrow the frame
    q_prior$x  <- sort(prior)
    post <- as_draws_df(workflow$posterior_samples[[pars]])$x
    q_post <- resample_draws(q_prior, tabulate(ecdf(prior)(post) * N))
    return (q_post) 
  } else if(sumtype == "sample"){ 
    return (subset_draws(post,variable = names(prior), iteration = 1))
  } else if(sumtype == "median"){ # median of M from each N prior 
  } else if(sumtype == "reweightsample"){
    post_df <- as_draws_df(posterior::subset_draws(post))
    w <- runif(ndraws(prior), 0, 10)
    return (resample_draws(post_df, w))
  }
}
post_summ_wf <- function(workflow, pars, sumtype){
  prior_rv <- workflow$prior_samples
  prior <- as_draws_df(workflow$prior_samples[[pars]])$x
  post <- as_draws_df(workflow$posterior_samples[[pars]])$x
  if (sumtype == "filtering"){ 
    if(length(prior) < 200){
      return (workflow$posterior_samples[names(prior_rv)])
    }else{
      q_prior <-as_draws_df(workflow$prior_samples[[pars]]) # to borrow the frame
      q_prior$x  <- sort(prior)
      post <- as_draws_df(workflow$posterior_samples[[pars]])$x
      q_post <- resample_draws(q_prior, tabulate(ecdf(prior)(post) * N))
      return (resample_draws(as_draws_df(workflow$prior_samples), tabulate(ecdf(prior)(post) * N))) 
    }
  } else if(sumtype == "sample"){ # only choose the first
    return (subset_draws(post,variable = names(prior), iteration = 1))
  } else if(sumtype == "median"){ # median of M from each N prior 
  } 
}

initDf <-function(L, summary, pars = NA){
  if(summary == "ms"){
    df <- data.frame(iter = 1:L,
                     mean = rep(NA,L),
                     sd = rep(NA,L)
    )
  }else if (summary == "q"){
    df <- data.frame(iter = 1:L,
                     q1 = rep(NA,L),
                     q3 = rep(NA,L)
    )
  }else if (summary == "pars"){
    df <- data.frame(iter = 1:L,
      median = rep(NA,L),
      mad = rep(NA,L)
      #tau_mean = rep(NA,L),
      #tau_sd = rep(NA,L)
    )
  }
  df
}


ws <- function(pri, post,  minBin = -1, maxBin = 1, numBins = 40){
  #min_s = min(c(pri,post)) reltarive range does not give globale entropy
  pribin = discretize(pri, numBins, r=c(minBin, maxBin)) / length(pri)
  postbin = discretize(post, numBins, r=c(minBin, maxBin))/ length(post)
  tb <- tibble(pribin, postbin)
  df <- melt(df, id.vars = 'iter')
  ggplot(df, aes(x = iter, y = value,  color = variable) ) +
    geom_point() + ggtitle(sprintf("priorSd: %s, N: %s, M: %s ", s0, N, M))
  bin_diff_abs <- Vectorize(function(i)  abs(pribin[i] - postbin[i]))
  wasserstein <- integrate(bin_diff_abs,1,numBins, rel.tol=.Machine$double.eps^.05)$value
  return (wasserstein - entropy::entropy(postbin))
}