loc_scale_post <- function(post_theta){
  return (list("loc" = t(apply(post_theta, c(2), median)), "scale" = t(apply(post_theta, c(2), sd))))
}

post_summ <- function(NMP_post_theta, sumtype){
  if (sumtype == "identity"){
      return (NMP_post_theta)
  } else if(sumtype == "resample"){
    NP_samp_post_theta <- sample_n(NMP_post_theta %>% group_by(prior, par), 1)  %>% select(par, value)
    return (NP_samp_post_theta) 
  } else if(sumtype == "median"){            
    return(NMP_tb %>%  gather(variable) %>% group_by(par, variable) %>% summarise('median'  = round(median(value),3)) %>% select(par, median))
  }
}

selfCalib <- function(modelName, sbc_obj, prior_theta, par_names, N, M, data, cnt){
  delivDir <- set_get_Dir(modelName)$delivDir
  prior_out = loc_scale_post(prior_theta) # sample summary for stan input
  data$param_loc = prior_out$loc
  data$param_scale = prior_out$scale
  sampled_y = sbc_obj$sample_y_tilde(prior_theta, data=data)
  post_theta = sbc_obj$sample_theta_bar_y(sampled_y, data=data, pars=par_names, fit_iter=M)
  NMP_post_theta <- as_tibble(post_theta)  %>%  gather(variable)
  NMP_post_theta$par <- gsub("\\..*","", NMP_post_theta$variable)
  NMP_post_theta$prior <- sapply(NMP_post_theta$variable, function(x){prior_theta[as.integer(gsub(".*\\.","",x)), gsub("\\..*","",x)]})
  postVar <- NMP_post_theta %>% group_by(par) %>% summarise('postVar' = sd(value)^2)

  priVar = NMP_post_theta %>% select(prior, par) %>% group_by(par) %>% unique() %>% summarise('priVar' = sd(prior)^2)
  ratioVar = inner_join(postVar, priVar) %>% summarise("ratioVar" =postVar/priVar)

  #summarize NM posterior samples to N for each parameter
  s_post_theta = NMP_post_theta %>% group_by(par, variable) %>% summarise('median'  = round(median(value),3)) %>% select(par, median)
  #as.matrix(s_post_theta %>% pivot_wider(names_from = par, values_from =mean)) #%>% summarise('mean'  = round(mean(value),3))
  
  #plot post by prior
  NMP_post_theta  %>% group_by(par, prior) %>% summarise("post_median" = median(value)) %>%
    ggplot(aes(x = prior, y = post_median, colour = par))  + geom_point() + xlim(min(prior_theta), max(prior_theta)) + ylim(min(prior_theta), max(prior_theta)) + coord_equal()

  medVarRatio <- round(median(ratioVar$ratioVar), 2)
  if(any(ratioVar$ratioVar<0.9) || any(ratioVar$ratioVar > 1.1)){
    cnt = cnt +1
    return (selfCalib(modelName, sbc_obj, post_summ(NMP_tb, sumtype = "resample"), par_names, N, M, data, cnt))
  }
  else{
    return (post_theta)
  }
}
