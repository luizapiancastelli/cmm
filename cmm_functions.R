dgo_h = function(pi,z){
  return(sum(z[pi] != sort(z)))
}

dgo_h_mat = function(Pi,z){
  d = apply(Pi, 1, dgo_h, z)
  return(d)
}

dgo_k = function(pi, z){
  return(dk_ties(group_from_pi(pi, z), group_from_z(z)))
}

dgo_k_mat = function(Pi,z){
  d = apply(Pi, 1, dk_cpp, z)
  return(d)
}

pair_matrix = function(group1, group2){
  r = length(group1)
  pairs = data.frame(expand.grid("i"=1:r, j = 1:r))
  for(index in 1:nrow(pairs)){
    pairs$value[index] = pair_count(index, pairs, group1, group2)
  }
  pairs = pairs %>% arrange(i)
  return(pairs)
}


dk_ties = function(group1, group2){
  pairs = pair_matrix(group1, group2)
  L = length(group1)
  indexes = expand.grid('i' =1:L, 'j' =1:L)
  indexes = indexes[indexes$i < L,  ]
  
  d =sum(sapply(1:nrow(indexes), term_, pairs = pairs))
  return(d)
}

rcmm = function(q, theta, z, dist, gibbs_iter = 1000){

if(dist == 'k'){
Pi = rcmm_k_gibbs_cpp(q, theta, z, gibbs_iter)
} else if(dist == 'h'){
Pi = rcmm_h_gibbs_cpp(q, theta, z, gibbs_iter)  
}

return(Pi)

}

#MLE 

MLE_fit_cmm = function(CT, Pi, theta_anneal, eps, dist, gibbs_iter = 200){
  
  if(dist == 'h'){
    run = MLE_fit_h(CT, Pi, theta_anneal, gibbs_iter, eps)
  } else if(dist== "k"){
    run = MLE_fit_k(CT, Pi, theta_anneal, gibbs_iter, eps)
  }
  
  return(run)
  
}


MLE_fit_h = function(nz, Pi, theta_anneal, gibbs_iter, eps){
  
  L = length(nz); z_init = c(); n = sum(nz)
  for(l in 1:L){
    z_init = c(z_init, rep(l,nz[l]))
  }
  z_init = sample(z_init, replace = F, size = n)
  
  opt_z =  z_annealing(z_init, Pi, 2, theta_anneal, "h")
  
  Obs_dist = mean(apply(Pi, 1, function(x, z){dgo_h(x, z)}, z= opt_z))
  
  theta_opt = opt_theta_h(Obs_dist, gibbs_iter, opt_z, eps)
  
  return(list("z_mle" = opt_z, 'theta_mle' = theta_opt, "Pi" = Pi, "dist" = "h"))
  
}


MLE_fit_k = function(nz, Pi, theta_anneal, gibbs_iter, eps){
  
  L = length(nz); z_init = c(); n = sum(nz)
  for(l in 1:L){
    z_init = c(z_init, rep(l,nz[l]))
  }
  z_init = sample(z_init, replace = F, size = n)
  
  opt_z =  z_annealing(z_init, Pi, 2, theta_anneal, "k")
  #opt_z =  z_annealing(z_init, Pi, 2, theta_anneal, "h")
  
  Obs_dist = mean(apply(Pi, 1, function(x, z){dk_cpp(x, z)}, z= opt_z))
  
  theta_opt = opt_theta_k(Obs_dist, gibbs_iter, opt_z, eps)
  
  return(list("z_mle" = opt_z, 'theta_mle' = theta_opt, "Pi" = Pi, "dist" = "k"))
  
}


z_annealing = function(z, Pi, m, theta_grid,dist){
  #z: initial value
  
  i =1
  while(i <= length(theta_grid)){
    
    MH = z_switch(z, theta_grid[i], Pi, m, dist)
    z = MH$z
    
    i = i +1
  }
  
  return(z)
  
}


opt_theta_k = function(Obs_dist, gibbs_iter, z, eps){
  
  theta_mle= 1; converge = FALSE
  while(converge ==FALSE){
    
    Expected_dist = Exp_appMC_k( 2000, theta_mle, z, gibbs_iter)
    
    #Expected distance is smaller: decrease theta
    new_theta = theta_mle*(Expected_dist/Obs_dist)
    
    converge =  abs(theta_mle -new_theta) < eps 
    
    theta_mle = new_theta

  }
  
  return(new_theta)
  
}


opt_theta_h = function(Obs_dist, gibbs_iter, z, eps){
  
  theta_mle= 1; converge = FALSE
  while(converge ==FALSE){
    
    Expected_dist = Exp_appMC_h( 1000, theta_mle, z, gibbs_iter)
    
    #Expected distance is smaller: decrease theta
    new_theta = theta_mle*(Expected_dist/Obs_dist)
    
    converge =  abs(theta_mle -new_theta) < eps 
    
    
    theta_mle = new_theta
    
  }
  
  return(new_theta)
  
}

Exp_appMC_h = function(MC, theta, z, gibbs_iter){
  sample = rcmm_h_gibbs_cpp(MC,   theta,  z, gibbs_iter)
  val =mean(apply(sample, 1, function(x, z){dgo_h(x, z)}, z= z))
  return(val)
}


Exp_appMC_k = function(MC, theta, z, gibbs_iter){
  sample = rcmm_k_gibbs_cpp(MC, theta, z, gibbs_iter)
  val =mean(apply(sample, 1, function(x, z){dk_cpp(x, z)}, z= z))
  return(val)
}



plot_rank_probs = function(Pi, z){
  probs_df = data.frame(rank_type_probs_data(Pi,z))
  
  probs_df$rank = factor(1:nrow(probs_df))
  df = melt(probs_df, id.vars = 'rank')
  names(df)[2:3] = c('cluster', 'p')
  df$cluster = gsub("X", "", df$cluster)
  
  plot = ggplot(df, aes(x = rank, y = p, fill = cluster))+geom_bar(stat = 'identity', position = 'stack', alpha = 0.7)+
    theme_bw()+theme(text = element_text(size = 16))+labs(x = 'Rank', y = 'Probability', fill = 'Cluster')
  return(plot)
  
}


rank_type_probs_data = function(Pi, z){
  z_Pi =t(apply(Pi, 1, function(x,z){z[x]},z))
  n = ncol(Pi)
  L = length(unique(z))
  
  mat = matrix(NA, ncol = L, nrow = n)
  
  i =1
  for(i in 1:n){
    l =1
    while(l <= L){
      mat[i, l] = mean( z_Pi[,i]==l )
      l= l +1
    }
    
  }
  
  return(mat)
}

Psi_IS_conditionals = function(N, z, theta, dist){
  
  reps = replicate(N,IS_conditionals_aux(z, theta, dist))
  return(mean(reps))
}

IS_conditionals_aux = function(z, theta, dist){
  sample = pseudo_sampler(z, theta)
  
  if(dist == 'h'){
    num = exp(-theta*dgo_h(sample$pi, z))
  } else if(dist == 'k'){
    num = exp(-theta* dk_cpp(sample$pi, z) )
  }
  is_ratio = num/prod(sample$probs)
  return(is_ratio)
}


Info_fit_h = function(fit,N_Psi){
  
  dists = apply(fit$Pi, 1, function(x, z){dgo_h(x, z)}, z= fit$z_mle)
  
  try =TRUE
  while(try){
    psi_est = Psi_IS_conditionals(N_Psi, fit$z_mle, fit$theta_mle, 'h')
    try = ifelse(is.na(psi_est), TRUE, FALSE)
  }
  
  loglik = -fit$theta_mle*sum(dists) - nrow(fit$Pi)*log(psi_est)
  #bic  = n_params(table(fit$z_mle))*log(nrow(fit$Pi)) - 2*loglik
  
  #table_penalisation = sum(lfactorial(table(fit$z_mle))) 
  
  n_perms = factorial(length(fit$z_mle))/prod(factorial(n_z(fit$z_mle, max(fit$z_mle))))
  
  log_n_perms = lfactorial(length(fit$z_mle)) - sum(lfactorial(table(fit$z_mle)))
  log_p_z = - log_n_perms
  
  val = loglik + log_p_z
  return(val)
}

Info_fit_k = function(fit,N_Psi){
  
  dists = apply(fit$Pi, 1, function(x, z){dk_cpp(x, z)}, z= fit$z_mle)
  
  try =TRUE
  while(try){
    psi_est = Psi_IS_conditionals(N_Psi, fit$z_mle, fit$theta_mle, 'k')
    try = ifelse(is.na(psi_est), TRUE, FALSE)
  }
  
  loglik = -fit$theta_mle*sum(dists) - nrow(fit$Pi)*log(psi_est)
  
  n_perms = factorial(length(fit$z_mle))/prod(factorial(n_z(fit$z_mle, max(fit$z_mle))))
  
  log_n_perms = lfactorial(length(fit$z_mle)) - sum(lfactorial(table(fit$z_mle)))
  log_p_z = - log_n_perms
  
  val = loglik + log_p_z
  return(c(val))
}

I_criterion = function(fit, N_Psi){
  
  if(fit$dist == 'k'){
    i =Info_fit_k(fit,N_Psi)
  } else if(fit$dist == 'h'){
    i = Info_fit_h(fit,N_Psi)
  }
  
  return(i)
}


CMM_MCMC_h = function(nz, Pi, niter, gibbs_iter, prior_theta){
  
  if(sum(is.na(Pi))>0){
    stop('NA values in Pi')
  }
  
  n = sum(nz)
  
  #Hamming case-------------
  #Initialise from MLEs
  mle_fit = MLE_fit_h(nz, Pi, seq(0.01, 3, 0.01), gibbs_iter, 10^(-3))
  theta_init = mle_fit$theta_mle
  z_init =mle_fit$z_mle
  
  z_chain = matrix(NA, ncol =n, nrow = niter)
  theta_chain = numeric(niter)
  accepted_z = 0; accepted_theta = 0;
  sigma_theta = 0.5 
  s=prior_theta$shape; r=prior_theta$rate;
  iter =1; theta = theta_init; z = z_init
  while(iter <= niter){
    
    #Z move (tractable):
    z_prop = switch_cpp(z)$prime
    num = -theta*sum(dgo_h_mat(Pi, z_prop))
    den = -theta*sum(dgo_h_mat(Pi, z))
    
    if(log(runif(1))< (num-den)){
      z = z_prop
      accepted_z = accepted_z + 1
    }
    
    #Theta move (AEA):
    step = theta_exchange_h(z, theta, Pi, sigma_theta, gibbs_iter,s,r)
    
    theta = step$theta
    accepted_theta = accepted_theta + step$accepted
    sigma_theta = adjust_proposal(sigma_theta, accepted_theta/iter, iter, 0.44)
    
    #Store values
    theta_chain[iter] = theta
    z_chain[iter, ] = z
    
    iter = iter +1
    if(iter %%100 == 0){
      print(iter)
    }
    
  }
  output = list("theta"=theta_chain, "z" = z_chain, 
                "accepted_theta" =accepted_theta/niter, dist = 'h', "Pi" =Pi )
  
}



theta_exchange_h = function(z, theta, Pi, sigma_theta, gibbs_iter,s,r){
  
  q=nrow(Pi); 
  
  theta_prime = rlnorm(1, meanlog = log(theta), sigma_theta)
  
  #Proposal ratio:
  q_ratio = dlnorm(theta, log(theta_prime), sigma_theta,log= T) -
    dlnorm(theta_prime, log(theta), sigma_theta,log= T)
  prior_ratio = dgamma(theta_prime, s, r, log = TRUE)-dgamma(theta, s, r, log = TRUE)
  
  #Exchange
  Pi_prime = rcmm_h_gibbs_cpp(q, theta_prime, z, gibbs_iter)
  
  num = -theta*sum(dgo_h_mat(Pi_prime, z))  -theta_prime*sum(dgo_h_mat(Pi, z))
  den = -theta_prime*sum(dgo_h_mat(Pi_prime, z)) - theta*sum(dgo_h_mat(Pi, z))
  log_alpha = num-den  + q_ratio+ prior_ratio
  
  if(log(runif(1))< log_alpha){
    theta = theta_prime
    accepted_theta = 1
    #if(theta>10){
    #  print( num-den )
    #  print( q_ratio+ prior_ratio )
    # }
  } else{
    theta = theta
    accepted_theta = 0
  }
  
  return(list('theta' = theta, 'accepted' = accepted_theta))
}

CMM_MCMC_k = function(nz, Pi, niter, gibbs_iter, prior_theta){
  
  if(sum(is.na(Pi))>0){
    stop('NA values in Pi')
  }
  
  n = sum(nz)
  
  #Hamming case-------------
  #Initialise from MLEs
  mle_fit = MLE_fit_h(nz, Pi, seq(0.01, 2, 0.01), gibbs_iter, 10^(-2))
  theta_init =  0.5
  z_init =mle_fit$z_mle
  
  z_chain = matrix(NA, ncol =n, nrow = niter)
  theta_chain = numeric(niter)
  accepted_z = 0; accepted_theta = 0;
  sigma_theta = 0.2
  s=prior_theta$shape; r=prior_theta$rate;
  iter =1; theta = theta_init; z = z_init
  while(iter <= niter){
    
    #Z move (tractable):
    z_prop = switch_cpp(z)$prime
    num = -theta*sum(dgo_k_mat(Pi, z_prop))
    den = -theta*sum(dgo_k_mat(Pi, z))
    
    if(log(runif(1))< (num-den)){
      z = z_prop
      accepted_z = accepted_z + 1
    }
    
    #Theta move (AEA):
    step = theta_exchange_k(z, theta, Pi, sigma_theta, gibbs_iter,s,r)
    
    theta = step$theta
    accepted_theta = accepted_theta + step$accepted
    sigma_theta = adjust_proposal(sigma_theta, accepted_theta/iter, iter, 0.3)
    
    #Store values
    theta_chain[iter] = theta
    z_chain[iter, ] = z
    
    iter = iter +1
    print(iter)
    
    
  }
  output = list("theta"=theta_chain, "z" = z_chain, 
                "accepted_theta" =accepted_theta/niter , dist = 'k', "Pi" = Pi)
  
}


theta_exchange_k = function(z, theta, Pi, sigma_theta, gibbs_iter,s,r){
  
  q=nrow(Pi); 
  theta_l_prime = rlnorm(1, meanlog = log(theta), sigma_theta)
  
  #Proposal ratio:
  q_ratio = dlnorm(theta, log(theta_l_prime), sigma_theta,log= T) -
    dlnorm(theta_l_prime, log(theta), sigma_theta,log= T)
  prior_ratio = dgamma(theta_l_prime, s, r, log = TRUE)-dgamma(theta, s, r, log = TRUE)
  
  #Exchange
  Pi_prime = rcmm_k_gibbs_cpp(q, theta_l_prime, z, gibbs_iter)
  
  num = -theta*sum(dgo_k_mat(Pi_prime, z))  -theta_l_prime*sum(dgo_k_mat(Pi, z))
  den = -theta_l_prime*sum(dgo_k_mat(Pi_prime, z)) - theta*sum(dgo_k_mat(Pi, z))
  log_alpha = num-den  + q_ratio+ prior_ratio
  
  if(log(runif(1))< log_alpha){
    theta = theta_l_prime
    accepted_theta = 1
    #if(theta>10){
    #  print( num-den )
    #  print( q_ratio+ prior_ratio )
    #}
  } else{
    theta = theta
    accepted_theta = 0
  }
  
  return(list('theta' = theta, 'accepted' = accepted_theta))
}

cmm_MCMC = function(CT, Pi, niter, dist, prior_theta, gibbs_iter = 100){
  
  if(dist == 'k'){
    run = CMM_MCMC_k(CT, Pi, niter, gibbs_iter, prior_theta)
  }else if(dist == 'h'){
    run = CMM_MCMC_h(CT, Pi, niter, gibbs_iter, prior_theta)
  }
  return(run)
  
}

z_switch = function(z, theta, Pi, switches, dist){
  
  if(sum(is.na(Pi))>0){
    stop('NA in Pi')
  }
  z_prime = switch_some(z, switches)
  
  if(dist == 'h'){
    alpha = -theta*(sum(apply(Pi,1, dgo_h_cpp, z_prime)) -sum(apply(Pi,1, dgo_h_cpp, z)))
  } else if(dist == 'k'){
    alpha = -theta*(sum(apply(Pi,1, dk_cpp, z_prime)) -sum(apply(Pi,1, dk_cpp, z)))
  }
  
  if(log(runif(1)) < alpha){ 
    z_out = z_prime
    accepted = ifelse(all(z==z_out),0, 1)
  } else {
    z_out = z
    accepted = 0
  }
  
  return(list('z' = z_out,'accepted' = accepted))
}

switch_some = function(x,m){
  x_switch = x
  indexes = sample(1:length(x), m)
  
  x_switch[indexes] =  sample(x[indexes])
  
  return(x_switch)
}

pseudo_sampler= function(z, theta){
  
  n = length(z)
  available = rep(TRUE, n) #Items: 1...n
  pi = numeric(n)
  probs = numeric(n)
  L = length(unique(z))
  
  tab = NULL
  step = 1
  z_sorted = sort(z)
  while(step <n){
    
    #tab = cbind(sort(unique(z[available])))
    tab = cbind(1:L, n_z(z[available], L))
    tab = tab[sort(unique(z[available])),]
    tab = matrix(tab, nrow = length(unique(z[available])))
    
    tab[which(tab[,1] ==  min(z[available])),2] = 0
    
    #tab = cbind(tab, c(0, rep(1,nrow(tab)-1)))
    tab = cbind(tab, rep(theta, nrow(tab)))
    
    p_un = exp(-tab[,3]*tab[,2])
    p = p_un/sum(p_un)
    tab = cbind(tab, p)
    
    if(nrow(tab) > 1){
      group_sampled = sample(tab[,1], size =1, prob = p)
      group_prob = tab[tab[,1] == group_sampled,4]
    }else{
      group_sampled = tab[1,1]
      group_prob = tab[1, 4]
    }
    
    group_filter = z == group_sampled
    av_filter = available
    filter = group_filter&av_filter
    item_set = (1:n)[filter]
    
    if(length(item_set)>1){
      pi[step] = sample(item_set, size =1)
      item_prob = 1/length(item_set)
    } else {
      pi[step] = item_set
      item_prob=1
    }
    
    
    probs[step] = group_prob*item_prob
    
    available[pi[step]]=FALSE
    step = step + 1
    
  }
  pi[n] = (1:n)[available]
  probs[n]=1
  
  return(list('pi'= pi, 'probs' = probs))
  
}

partition_from_z = function(z, item_names = 1:length(z)){
  L =length(unique(z))
  groups = list();
  for(l in 1:L){
    groups = append(groups, list(item_names[which(z==l)]  ))
  }
  return(groups)
}



