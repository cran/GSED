library(rootSolve)
library(memoise) 
library(survival)

catch_entries_commun = function(K_stages, N_subsets, f, ratio_Delta_star_d1, ordering, increasing_theta, seed){
  if(as.integer(K_stages) != K_stages){
    K_stages = as.integer(K_stages)
    print(paste("The number of stages was transformed into integer: K_stages=", K_stages, sep=""))
  }
  if(K_stages<1){
    stop("The number of stages must be > 0")
  }
  if(as.integer(N_subsets) != N_subsets){
    N_subsets = as.integer(N_subsets)
    print(paste("The number of distinct subgroups was transformed into integer: N_subsets=", N_subsets, sep=""))
  }
  if(N_subsets<1){
    stop("The number of subgroups must be > 0")
  }
  if(length(f) != N_subsets){
    stop(paste("The prevalence rates, f, must be a vector of length: ", N_subsets, sep=""))
  }
  if(any(f>1) || any(f<0)){
    stop("Prevalance rates must be between 0 and 1")
  }
  if(length(ratio_Delta_star_d1) != K_stages-1){
    stop(paste("The ratio between information at stage k (>1) with information at stage 1 must be a vectors of length: ", K_stages-1, sep=""))
  }
  if(is.logical(ordering)==FALSE){
    stop("ordering indicating if subgroup are assumed ordered must be a boolean")
  } 
  if(is.logical(increasing_theta)==FALSE){
    stop("increasing_theta indicating if response is assumed better with increasing parameter must be a boolean")
  } 
  if(as.integer(seed) != seed){
    seed = as.integer(seed)
    print(paste("The seed was transformed into integer: seed=", seed, sep=""))
  }
  return(list(K_stages, N_subsets, f, ratio_Delta_star_d1, ordering, increasing_theta, seed))
}


catch_entries_boundaries = function(K_stages, N_subsets, f, ratio_Delta_star_d1, ordering, increasing_theta, seed, n_trials, alpha_spending, one_minus_alpha_spending){
  catch = catch_entries_commun(K_stages, N_subsets, f, ratio_Delta_star_d1, ordering, increasing_theta, seed)
  K_stages = catch[[1]]
  N_subsets = catch[[2]]
  f = catch[[3]]
  ratio_Delta_star_d1 = catch[[4]]
  ordering = catch[[5]]
  increasing_theta = catch[[6]]
  seed = catch[[7]]
  if(as.integer(n_trials) != n_trials){
    n_trials = as.integer(n_trials)
    print(paste("The number of simulated trials was transformed into integer: n_trials=", n_trials, sep=""))
  }
  if(length(alpha_spending) != K_stages+1 || length(one_minus_alpha_spending) != K_stages+1){
    stop(paste("alpha_spending and one_minus_alpha_spending must be vectors of length: ", K_stages+1, " representing alpha spending and 1-alpha spending respectively, with values 0 in first position and values alpha and 1-alpha respectively at position K_stages+1", sep=""))
  }
  if(any(alpha_spending>1) || any(alpha_spending<0) || any(one_minus_alpha_spending>1) || any(one_minus_alpha_spending<0)){
    stop("Values for alpha_spending and one_minus_alpha_spending must be between 0 and 1")
  }
  return(list(K_stages, N_subsets, f, ratio_Delta_star_d1, ordering, increasing_theta, seed, n_trials, alpha_spending, one_minus_alpha_spending))
}


catch_entries_FI = function(K_stages, N_subsets, f, ratio_Delta_star_d1, l, u, type_outcome, param_theta, pow, ordering, increasing_theta, seed, n_trials, rule){
  catch = catch_entries_commun(K_stages, N_subsets, f, ratio_Delta_star_d1, ordering, increasing_theta, seed)
  K_stages = catch[[1]]
  N_subsets = catch[[2]]
  f = catch[[3]]
  ratio_Delta_star_d1 = catch[[4]]
  ordering = catch[[5]]
  increasing_theta = catch[[6]]
  seed = catch[[7]]
  if(length(l) != K_stages || length(u) != K_stages){
    stop(paste("The lower and upper boundaries for decision stages must be two vectors of length: ", K_stages, sep=""))
  }
  if(as.integer(n_trials) != n_trials){
    n_trials = as.integer(n_trials)
    print(paste("The number of simulated trials was transformed into integer: n_trials=", n_trials, sep=""))
  }
  if(pow<0 || pow>1){
    stop("Power, pow, must be between 0 and 1")
  }
  if(as.integer(rule) != rule){
    rule = as.integer(rule)
    print(paste("The rule number for power calculation was transformed into integer: rule=", rule, sep=""))
  }
  if(rule != 1 && rule != 2){
    stop("Only rule number 1 or 2 are implemented")
  }
  if(type_outcome != "binary" && type_outcome != "continuous" && type_outcome != "survival"){
    stop("Type of outcome not supported")
  }
  return(list(K_stages, N_subsets, f, ratio_Delta_star_d1, l, u, param_theta, pow, ordering, increasing_theta, seed, n_trials, rule, type_outcome))
}



catch_entries_MT = function(K_stages, N_subsets, f, l, u, ratio_Delta_star_d1, type_outcome, param_outcome=NA,
                         n_max=NA, incl_rate=NA, mean_cur_c=NA, HR=NA, 
                         nb_required=NA, nmax_wait=+Inf, ordering, increasing_theta=FALSE, nsim=1000, seed=42){
  catch = catch_entries_commun(K_stages, N_subsets, f, ratio_Delta_star_d1, ordering, increasing_theta, seed)
  K_stages = catch[[1]]
  N_subsets = catch[[2]]
  f = catch[[3]]
  ratio_Delta_star_d1 = catch[[4]]
  ordering = catch[[5]]
  increasing_theta = catch[[6]]
  seed = catch[[7]]
  if(length(l) != K_stages || length(u) != K_stages){
    stop(paste("The lower and upper boundaries for decision stages must be two vectors of length: ", K_stages, sep=""))
  }
  if(type_outcome != "survival" && type_outcome != "binary" && type_outcome != "continuous"){
    stop("The type of outcome must be either: survival, binary or continuous")
  }
  if(is.na(param_outcome)==FALSE){
    if(type_outcome=="binary"){
      if(length(param_outcome) != 1){
        stop("The parameters supplied for the binary outcome, param_outcome, must be a list of one element containing a matrix of size 2xN_subsets. The matrix should contain the probabilities of response for in row control or treatment, and in column the subgroup number")
      }
      else{
        if(N_subsets==1){
          if(is.matrix(param_outcome[[1]])==FALSE){
            param_outcome[[1]] = t(as.matrix(param_outcome[[1]]))
          }
          if(any(param_outcome[[1]]>1) || any(param_outcome[[1]]<0)){
            stop("Response probabilty in param_outcome must be between 0 and 1")
          }
        }
        else{
          if(all(dim(param_outcome[[1]]) == c(2,N_subsets))==FALSE){
            stop("In param_outcome, the matrix inside the list containing the probabilities of response for each treatment x subgroup is of wrong size")
          }
          if(any(param_outcome[[1]]>1) || any(param_outcome[[1]]<0)){
            stop("Response probabilties in param_outcome must be between 0 and 1")
          }
        }
      }
    }
    else if(type_outcome=="continuous"){
      if(length(param_outcome) != 2){
        stop("The parameters supplied for the continuous outcome, param_outcome, must be a list of two elements containing two matrices of size 2xN_subsets. The two matrices should contain the means and variances respectively, for in row control or treatment, and in column the subgroup number")
      }
      else{
        if(N_subsets==1){
          if(is.matrix(param_outcome[[1]])==FALSE){
            param_outcome[[1]] = t(as.matrix(param_outcome[[1]]))
          }
          if(is.matrix(param_outcome[[2]])==FALSE){
            param_outcome[[2]] = t(as.matrix(param_outcome[[2]]))
          }
        }
        else{
          if(all(dim(param_outcome[[1]]) == c(2,N_subsets))==FALSE || all(dim(param_outcome[[2]]) == c(2,N_subsets))==FALSE){
            stop("In param_outcome, the matrices inside the lists containing means or variances for each treatment x subgroup are of wrong size")
          }
        }
      }
    }
    else{
      stop("With the provided type of outcome, param_outcome should not be provided but set as NA")
    }
  }
  else{
    if(type_outcome=="binary" || type_outcome=="continuous"){
      stop("param_outcome should be provided. For binary outcome, it must be a list of size one containing a matrix of size 2xN_subsets. The matrix should contain the probabilities of response for in row control or treatment, and in column the subgroup number. For continuous outcome, it must be a list of size two containing two matrices of size 2xN_subsets. The two matrices should contain the means and variances respectively, for in row control or treatment, and in column the subgroup number")
    }
  }
  if(is.na(n_max) == FALSE && as.integer(n_max) != n_max){
    n_max = as.integer(K_stages)
    print(paste("The maximum number of patients was transformed into integer: n_max=", n_max, sep=""))
  }
  if((type_outcome=="binary" || type_outcome=="continuous") && is.na(n_max)){
    stop("The maximum number of patients to enroll, n_max, must be supplied")
  }
  if(type_outcome=="survival" && is.na(n_max)==FALSE){
    print("As type of outcome is survival, the maximum number of patients n_max should be NA and will be ignored")
  }
  if(type_outcome=="survival" && is.na(incl_rate)){
    stop("The inclusion rate, incl_rate, must be provided")
  }
  if(type_outcome=="survival" && is.na(mean_cur_c)){
    stop("The median survival for control group, mean_cur_c, must be provided")
  }
  if(type_outcome=="survival" && anyNA(HR)){
    stop("The expected hazard ratios for each subgroup, HR, must be provided")
  }
  if(type_outcome=="survival" && is.na(nb_required)){
    stop("The maximum number of events required, nb_required, must be provided")
  }
  if(type_outcome=="survival" && is.na(nmax_wait)== FALSE && (nmax_wait < nb_required)){
    stop("The maximum number of patients to include in the trial must be superior or equal to the number of events required")
  }
  if(type_outcome=="survival" && length(HR) != N_subsets){
    stop(paste("The expected hazard ratios for each subgroup, HR, must be a vector of length: ", N_subsets, sep=""))
  }
  if(type_outcome=="survival" && as.integer(nb_required) != nb_required){
    nb_required = as.integer(nb_required)
    print(paste("The maximum number of events required was transformed into integer: nb_required=", nb_required, sep=""))
  }
  if((type_outcome=="binary" || type_outcome=="continuous") && 
     (is.na(incl_rate)==FALSE || is.na(mean_cur_c)==FALSE || 
      anyNA(HR)==FALSE || is.na(nb_required)==FALSE)){
    print("Arguments required only for survival outcome are supplied and will be ignored")
  }
  if(as.integer(nsim) != nsim){
    nsim = as.integer(nsim)
    print(paste("The number of simulations to perform was transformed into integer: nsim=", nsim, sep=""))
  }
  return(list(K_stages, N_subsets, f, l, u, ratio_Delta_star_d1, type_outcome, param_outcome, n_max, incl_rate, 
              mean_cur_c, HR, nb_required, ordering, increasing_theta, nsim, seed, nmax_wait))
}


#########
### STAGE 1 - Subpopulation selection
stage_1_selection = function(N_subsets, Z_1j, l, ordering, increasing_theta=FALSE){
  keep = c()
  if(ordering == FALSE){
    for(j in 1:N_subsets){
      if(Z_1j[j] > l[1]){
        keep = c(keep, j)
      }
    }
  }
  else{
    r = NA
    if(increasing_theta==FALSE){
      j = N_subsets
      while(is.na(r) && j >= 1){
        if(Z_1j[j] > l[1]){
          r = j
        }
        j = j-1
      }
      if(!is.na(r)){
        keep = 1:r
      }
    }
    else{
      j = 1
      while(is.na(r) && j <= N_subsets){
        if(Z_1j[j] > l[1]){
          r = j
        }
        j = j+1
      }
      if(!is.na(r)){
        keep = r:N_subsets
      }
    }
  } 
  return(keep)
}


#########
### STAGE 1 - Evaluation of subpopulation selected
stage_1_evaluation = function(keep, Z_1j, f, u){
  if(length(keep)==0){
    return(list("stage"=1,"S"=c()))
  }
  else{
    f_keep = f[keep] 
    f_S = sum(f_keep)
    Z_1S = 1/sqrt(f_S) * sum( Z_1j[keep] * sqrt(f_keep) )
    if(Z_1S > u[1]){
      return(list("stage"=1,"S"=keep))
    }
    else{
      return(list("stage"=-1,"S"=keep))
    }
  }
}


#########
### SUBSEQUENT STAGES for global Z stat (shortcut in simulations for time saving)
subsequent_stages_sim = function(stage_k, f, keep, tZ, u, l, part_Z_prev, ratio_Delta_star_d1){
  f_S = sum(f[keep])
  part_Z_prev = part_Z_prev + tZ
  sum_ratios = sum(ratio_Delta_star_d1[1:(stage_k-1)])
  Z_kS = 1/sqrt(f_S+sum_ratios) * part_Z_prev
  if(Z_kS > u[stage_k]){
    return(list("stage"=stage_k,"S"=keep)) # stop at stage k with rejection
  }
  else if(Z_kS <= l[stage_k]){
    return(list("stage"=stage_k,"S"=c())) # stop at stage k with acceptation
  }
  else{
    return(list("stage"=-1,"S"=keep,"part_Z_prev"=part_Z_prev))
  }
}


#########
### Differences in spending functions and observed fisher informations
diff_spending = function(K_stages, alpha_spending, one_minus_alpha_spending){
  diff_alpha_U = numeric(K_stages)
  diff_alpha_L = numeric(K_stages)
  for(k in 1:K_stages){
    diff_alpha_U[k] = alpha_spending[k+1] - alpha_spending[k]
    diff_alpha_L[k] = one_minus_alpha_spending[k+1] - one_minus_alpha_spending[k]
  } 
  return(list(diff_alpha_U, diff_alpha_L))
}


#########
### Simulation of one trial for boundaries determination
# 0 stop for futility, 1 stop for efficacy
sim_one_trial_boundaries = function(K_stages, N_subsets, f, ratio_Delta_star_d1, u, l, ordering, increasing_theta=FALSE){
  ### STAGE 1
  Z_1j = rnorm(N_subsets, 0, 1)
  subpop = stage_1_selection(N_subsets, Z_1j, l, ordering, increasing_theta)
  eval_s1 = stage_1_evaluation(subpop, Z_1j, f, u)
  keep = eval_s1$S
  if(eval_s1$stage==1){
    return(c(1,as.numeric(!is.null(keep))))
  }
  else{
    if(K_stages == 1){
      return(c(1,NA))
    }
    ### SUBSEQUENT STAGES
    else{
      f_keep = f[keep] 
      f_S = sum(f_keep)
      part_Z_prev = sum( Z_1j[keep] * sqrt(f_keep) )
      for(k in 2:K_stages){
        tZ = rnorm(1, 0, sqrt(ratio_Delta_star_d1[k-1]))
        eval_sk = subsequent_stages_sim(k, f, keep, tZ, u, l, part_Z_prev, ratio_Delta_star_d1)    
        if(eval_sk$stage==k){
          return(c(k,as.numeric(!is.null(eval_sk$S))))
        }
        else{
          part_Z_prev = eval_sk$part_Z_prev
          if(k==K_stages){
            return(c(k,NA))
          }
        }
      }
    }
  }
}


#########
### Simulations of n_trials for boundaries determination
sim_trials_boundaries = function(K_stages, N_subsets, f, ratio_Delta_star_d1, u, l, ordering, increasing_theta=FALSE, seed=42, n_trials){
  set.seed(seed)
  prop_eff_k = numeric(K_stages)
  prop_fut_k = numeric(K_stages)
  for(i in 1:n_trials){
    sim_i = sim_one_trial_boundaries(K_stages, N_subsets, f, ratio_Delta_star_d1, u, l, ordering, increasing_theta)     
    stage_i = sim_i[1]
    rule_stop_i = sim_i[2]
    if(!is.na(rule_stop_i)){
      if(rule_stop_i == 1){
        prop_eff_k[stage_i] = prop_eff_k[stage_i]+1
      }
      else{
        prop_fut_k[stage_i] = prop_fut_k[stage_i]+1
      }
    }
  }
  prop_eff_k = prop_eff_k/n_trials
  prop_fut_k = prop_fut_k/n_trials
  return(list(prop_eff_k,prop_fut_k))
}


#########
### Determine stopping BOUNDARIES
boundaries_sim = function(K_stages, N_subsets, f, ratio_Delta_star_d1, ordering, increasing_theta=FALSE, seed=42, n_trials,
                          alpha_spending, one_minus_alpha_spending){
  catch = catch_entries_boundaries(K_stages, N_subsets, f, ratio_Delta_star_d1, ordering, increasing_theta, seed, n_trials, alpha_spending, one_minus_alpha_spending)
  K_stages = catch[[1]]
  N_subsets = catch[[2]]
  f = catch[[3]]
  ratio_Delta_star_d1 = catch[[4]]
  ordering = catch[[5]]
  increasing_theta = catch[[6]]
  seed = catch[[7]]
  n_trials = catch[[8]]
  alpha_spending = catch[[9]]
  one_minus_alpha_spending = catch[[10]]
  
  bound_U = c()
  bound_L = c()
  diff_spend = diff_spending(K_stages, alpha_spending, one_minus_alpha_spending)   
  for(k in 1:K_stages){
    fun_l = function(x) {
      sim = sim_trials_boundaries(k, N_subsets, f, ratio_Delta_star_d1, c(bound_U,+Inf), c(bound_L,x), ordering, increasing_theta, seed, n_trials)
      return(sim[[2]][k] - diff_spend[[2]][k])
    }
    sol_l = uniroot(fun_l, c(-50,50))$root
    bound_L = c(bound_L,sol_l)
    fun_u = function(x) {
      sim = sim_trials_boundaries(k, N_subsets, f, ratio_Delta_star_d1, c(bound_U,x), c(bound_L), ordering, increasing_theta, seed, n_trials)
      return(sim[[1]][k] - diff_spend[[1]][k])
      
    }
    sol_u = uniroot(fun_u, c(-50,50))$root
    bound_U = c(bound_U,sol_u)
    if(k == K_stages){
      m = mean(c(bound_U[k],bound_L[k]))
      bound_U[k] = m
      bound_L[k] = m
    }    
  }
  return(list("l"=round(bound_L,4), "u"=round(bound_U,4))) 
}


#########
### Simulation of one trial for maximum Fisher Information determination
sim_one_trial_max_FI = function(K_stages, N_subsets, f, ratio_Delta_star_d1, l, u, type_outcome, param_theta, Imax, ordering, increasing_theta=FALSE){
  if(type_outcome == "binary"){
    #p_mean = (param_theta[[1]][2,]+param_theta[[1]][1,])/2
    theta = (param_theta[[1]][2,]-param_theta[[1]][1,]) #/ sqrt(p_mean*(1-p_mean))
  }
  else if(type_outcome == "continuous"){
    theta = (param_theta[[1]][2,]-param_theta[[1]][1,]) #/ sqrt( (param_theta[[2]][2,]+param_theta[[2]][1,])/2 )
  }
  else if(type_outcome == "survival"){
    theta = -log(param_theta)
  }
  
  ### STAGE 1
  Z_1j = rep(NA,N_subsets) 
  for(j in 1:N_subsets){
    Z_1j[j] = rnorm(1, theta[j]*sqrt(f[j]*Imax/(1+sum(ratio_Delta_star_d1))), 1)
  }
  subpop = stage_1_selection(N_subsets, Z_1j, l, ordering, increasing_theta)
  eval_s1 = stage_1_evaluation(subpop, Z_1j, f, u)
  keep = eval_s1$S
  if(eval_s1$stage==1){
    return(list(1,as.numeric(!is.null(keep)),keep))
  }
  else{
    if(K_stages == 1){
      return(list(1,NA,keep))
    }
    ### SUBSEQUENT STAGES
    else{
      f_keep = f[keep] 
      f_S = sum(f_keep)
      part_Z_prev = sum( Z_1j[keep] * sqrt(f_keep) )
      for(k in 2:K_stages){
        tZ = rnorm(1, ratio_Delta_star_d1[k-1] * sqrt(Imax/(1+sum(ratio_Delta_star_d1))) * sum(theta[keep]*f_keep), sqrt(ratio_Delta_star_d1[k-1]))
        eval_sk = subsequent_stages_sim(k, f, keep, tZ, u, l, part_Z_prev, ratio_Delta_star_d1)    
        if(eval_sk$stage==k){
          return(list(k,as.numeric(!is.null(eval_sk$S)),keep))
        }
        else{
          part_Z_prev = eval_sk$part_Z_prev
          if(k==K_stages){
            return(list(k,NA,keep))
          }
        }
      }
    }
  }
}


#########
### Simulations of n_trials for maximum Fisher Information determination
sim_trials_max_FI = function(K_stages, N_subsets, f, ratio_Delta_star_d1, l, u, type_outcome, param_theta, Imax, ordering, increasing_theta=FALSE, seed=42, n_trials){
  set.seed(seed)
  prop_eff_k = numeric(K_stages)
  prop_fut_k = numeric(K_stages)
  prop_eff_all_k = numeric(K_stages)
  for(i in 1:n_trials){
    sim_i = sim_one_trial_max_FI(K_stages, N_subsets, f, ratio_Delta_star_d1, l, u, type_outcome, param_theta, Imax, ordering, increasing_theta)     
    stage_i = sim_i[[1]]
    rule_stop_i = sim_i[[2]]
    keep = sim_i[[3]]
    if(!is.na(rule_stop_i)){
      if(rule_stop_i == 1){
        prop_eff_k[stage_i] = prop_eff_k[stage_i]+1
        if(all.equal(sort(keep),1:N_subsets)==TRUE){
          prop_eff_all_k[stage_i] = prop_eff_all_k[stage_i]+1
        }
      }
      else{
        prop_fut_k[stage_i] = prop_fut_k[stage_i]+1
      }
    }
  }
  prop_eff_k = prop_eff_k/n_trials
  prop_fut_k = prop_fut_k/n_trials
  prop_eff_all_k = prop_eff_all_k/n_trials
  return(list(prop_eff_k,prop_fut_k,prop_eff_all_k))
}


#########
### Maximum Fisher Information
max_FI = function(K_stages, N_subsets, f, ratio_Delta_star_d1, l, u, type_outcome, param_theta, pow, ordering, increasing_theta=FALSE, seed=42, n_trials, rule){  
  catch = catch_entries_FI(K_stages, N_subsets, f, ratio_Delta_star_d1, l, u, type_outcome, param_theta, pow, ordering, increasing_theta, seed, n_trials, rule)
  K_stages = catch[[1]]
  N_subsets = catch[[2]]
  f = catch[[3]]
  ratio_Delta_star_d1 = catch[[4]]
  l = catch[[5]]
  u = catch[[6]]
  param_theta = catch[[7]]
  pow = catch[[8]]
  ordering = catch[[9]]
  increasing_theta = catch[[10]]
  seed = catch[[11]]
  n_trials = catch[[12]]
  rule = catch[[13]]
  type_outcome = catch[[14]]
  
  if(rule==1){
    fun_FI = function(x){
      sim = sim_trials_max_FI(K_stages, N_subsets, f, ratio_Delta_star_d1, l, u, type_outcome, param_theta, x, ordering, increasing_theta, seed, n_trials) 
      return(sum(sim[[3]])-pow)
    }
  }
  else if(rule==2){
    fun_FI = function(x){
      sim = sim_trials_max_FI(K_stages, N_subsets, f, ratio_Delta_star_d1, l, u, type_outcome, param_theta, x, ordering, increasing_theta, seed, n_trials) 
      return(sum(sim[[1]])-pow)
    }
  }
  FI_max = uniroot(fun_FI, c(0,1e08))$root
  return(FI_max) 
}


#########
### Application of Magnusson and Turnbull with data
magnusson_turnbull = function(stage_cur, keep=NA, N_subsets, Y, I, l, u, ordering, increasing_theta=FALSE){
  Z = Y/sqrt(I)
  if(stage_cur==0){
    keep = stage_1_selection(N_subsets, Z, l, ordering, increasing_theta)
    if(length(keep) == 0){
      return(list("Rejection"=0, "Acceptation"=1, "Keep"="No subgroup"))
    }
    else{
      return(list("Rejection"=0, "Acceptation"=0, "Keep"=keep))
    }
  }
  else{
    if(stage_cur==1){ 
      if(Z > u[1]){
        return(list("Rejection"=1, "Acceptation"=0, "Keep"=keep))
      }
      else{
        return(list("Rejection"=0, "Acceptation"=0, "Keep"=keep))
      }
    }
    else if(stage_cur > 1){
      if(Z > u[stage_cur]){
        return(list("Rejection"=1, "Acceptation"=0, "Keep"=keep))
      }
      else if(Z <= l[stage_cur]){
        return(list("Rejection"=0, "Acceptation"=1, "Keep"="No subgroup"))
      }
      else{
        return(list("Rejection"=0, "Acceptation"=0, "Keep"=keep))
      }
    }
  } 
}


#########
### Is a vector an element in a list
in_list = function(vec,liste){
  res = FALSE
  i=1
  while(i <= length(liste) && res == FALSE){
    if(setequal(vec,liste[[i]])){
      res = TRUE
    }
    i=i+1
  }
  return(list(res,i-1))
}


## Functions for sampling some random distributions by batches for
## better performance
rexp1 = function(p) {
  b = rexp(100, p)
  i = 0
  return(function() {
    if(i >= 100) {
      b <<- rexp(100, p)
      i <<- 0
    }
    i <<- i + 1
    return(b[i])
  })
}

rbern1 = function(p) {
  b = rbinom(100, 1, p)
  i = 0
  return(function() {
    if(i >= 100) {
      b <<- rbinom(100, 1, p)
      i <<- 0
    }
    i <<- i + 1
    return(b[i])
  })
}

runif1 = function() {
  b = runif(100, 0, 1)
  i = 0
  return(function() {
    if(i >= 100) {
      b <<- runif(100, 0, 1)
      i <<- 0
    }
    i <<- i + 1
    return(b[i])
  })
}

sample1 = function(pop, prob) {
  b = sample(pop, 100, TRUE, prob)
  i = 0
  return(function() {
    if(i >= 100) {
      b <<- sample(pop, 100, TRUE, prob)
      i <<- 0
    }
    i <<- i + 1
    return(b[i])
  })
}


### Simulation of patients' inclusions for survival data until number of events required
incl_until_ev = function(incl_rate, stage, nb_required, nmax_wait, data_1, f, keep, mean_cur_c, HR){
  if(stage == 1){
    data = list("nb_ev" = 0, "duration"=0, "dur_incl"=0, "t_incl"=c(), "t_event"=c(),
                "trt"=c(), "biom_group"=c(), "n_pat"=0, "n_pat_keep"=0)
  }
  else{
    data = data_1
  }
  
  file = data$t_event[data$t_event > data$duration]
  n_pat_keep = data$n_pat_keep
  rincl = rexp1(incl_rate)
  rsurv = rexp1(1/mean_cur_c)
  rtrt = rbern1(0.5)
  rbg = sample1(keep, f)
  runi = runif1()
  next_incl = data$duration + rincl()
  
  while(data$nb_ev < nb_required){
    next_ev = if(length(file) == 0) Inf else min(file)
    if(data$n_pat < nmax_wait && next_incl < next_ev){
      data$duration = next_incl
      bg = rbg()
      if(is.element(bg,keep)){
        data$dur_incl = next_incl
        data$n_pat = data$n_pat + 1
        n_pat_keep = n_pat_keep + 1
        data$t_incl[n_pat_keep] = next_incl
        data$biom_group[n_pat_keep] = bg
        data$trt[n_pat_keep] = rtrt()
        if(data$trt[n_pat_keep] == 0){
          data$t_event[n_pat_keep] = next_incl + rsurv()
        }
        else{
          data$t_event[n_pat_keep] = next_incl + rsurv()/HR[bg]
        }
        file[length(file)+1] = data$t_event[n_pat_keep]
      }
      next_incl = next_incl + rincl()
    }
    else{
      i = which.min(file)
      data$duration = file[i]
      file = file[-i]
      data$nb_ev = data$nb_ev+1
    }
  }
  data$n_pat_keep = n_pat_keep
  return(data)
}


### Calculation of survival quantities
up_survival = function(data){
  follow = data$duration - data$t_incl
  time_event = data$t_event - data$t_incl
  ind_event = as.numeric(time_event <= follow)
  time_min = pmin(time_event, follow, na.rm=TRUE)
  return(list("time_min"=time_min, "ind_event"=ind_event))
}


remove_NS = function(data, keep, nb_ev_S){
  data_next = list()
  data_next$nb_ev = nb_ev_S
  data_next$duration = data$duration
  data_next$dur_incl = data$dur_incl
  data_next$n_pat = data$n_pat
  ind_keep = which(is.element(data$biom_group,keep))
  data_next$n_pat_keep = length(ind_keep)
  data_next$t_incl = data$t_incl[ind_keep]
  data_next$t_event = data$t_event[ind_keep]
  data_next$trt = data$trt[ind_keep]
  data_next$biom_group = data$biom_group[ind_keep]
  return(data_next)
}



#########
### Simulation of one design with Magnusson and Turnbull for survival data 
sim_one_OS_MT = function(K_stages, N_subsets, f, l, u, ratio_Delta_star_d1, incl_rate, 
                         mean_cur_c, HR, nb_required, nmax_wait=+Inf, ordering, increasing_theta=FALSE){
  reject = NA
  
  n_1 = ceiling(nb_required / (1+sum(ratio_Delta_star_d1)))
  n_req_step = c(n_1, nb_required-n_1)
  
  #Step 1
  data = incl_until_ev(incl_rate=incl_rate, stage=1, nb_required=n_req_step[1], nmax_wait=nmax_wait, data_1=NA, 
                       f=f, keep=1:N_subsets, mean_cur_c=mean_cur_c, HR=HR)
  data_k = data
  update = up_survival(data)
  time_min = update$time_min
  ind_event = update$ind_event
    
  Y_1j = numeric(N_subsets)
  I_1j = numeric(N_subsets)
  for(j in 1:N_subsets){
    ind_group_j = which(data$biom_group==j)
    nb_event_j = sum(ind_event[ind_group_j])
    if(length(which(data$trt[ind_group_j]==0))>0 && length(which(data$trt[ind_group_j]==1))>0 && nb_event_j>0){
      temp = survdiff(Surv(time_min[ind_group_j], ind_event[ind_group_j]) ~ data$trt[ind_group_j])
      Z_1j = (temp$obs[1]-temp$exp[1])/sqrt(temp$var[1,1])
      I_1j[j] = nb_event_j/4
      Y_1j[j] = Z_1j*sqrt(I_1j[j])
    }
    else{
      Y_1j[j] = NA
      I_1j[j] = 1
    }
  }
  
  k=1
  selection = magnusson_turnbull(0, NA, N_subsets, Y_1j, I_1j, l, u, ordering, increasing_theta)
  if(selection$Acceptation == 1){
    reject = 0
    keep = c()
  }
  else{
    keep = selection$Keep
    ind_group_1S = c()
    for(j in keep){
      ind_group_1S = c(ind_group_1S, which(data$biom_group==j))
    }
    nb_event_1S = sum(ind_event[ind_group_1S])
    if(length(which(data$trt[ind_group_1S]==0))>1 && length(which(data$trt[ind_group_1S]==1))>1 && nb_event_1S>0){
      temp = survdiff(Surv(time_min[ind_group_1S], ind_event[ind_group_1S]) ~ data$trt[ind_group_1S])
      Zp_1S = (temp$obs[1]-temp$exp[1])/sqrt(temp$var[1,1])
      I_1S = nb_event_1S/4
      Y_1S = Zp_1S*sqrt(I_1S)
    }
    else{
      Y_1S = NA
      I_1S = 1
    }
    step1 = magnusson_turnbull(1, keep, N_subsets, Y_1S, I_1S, l, u, ordering, increasing_theta) 
    if(step1$Rejection == 1){
      reject = 1
    }  
    else{  
      npat_1S = length(which(is.element(data$biom_group,keep)==TRUE))
      npat_1NS = length(which(is.element(data$biom_group,keep)==FALSE))
      data = remove_NS(data, keep, nb_event_1S)
      k = 2
      
      while(k <= K_stages && is.na(reject)){
        data_k = incl_until_ev(incl_rate=incl_rate, stage=k, nb_required=cumsum(n_req_step)[k], 
                               nmax_wait=max(nmax_wait, npat_1NS+cumsum(n_req_step)[k]), data_1=data, 
                               f=f, keep=keep, mean_cur_c=mean_cur_c, HR=HR)
        update = up_survival(data_k)
        time_min = update$time_min
        ind_event = update$ind_event
        nb_event_kS = sum(ind_event)
        if(length(which(data_k$trt==0))>1 && length(which(data_k$trt==1))>1 && nb_event_kS>0){
          temp = survdiff(Surv(time_min, ind_event) ~ data_k$trt)
          Zp_kS = (temp$obs[1]-temp$exp[1])/sqrt(temp$var[1,1])
          I_kS = nb_event_kS/4
          Y_kS = Zp_kS*sqrt(I_kS)
        }
        else{
          Y_kS = NA
          I_kS = 1
        }
        stepk = magnusson_turnbull(k, keep, N_subsets, Y_kS, I_kS, l, u, ordering, increasing_theta)
        if(stepk$Rejection == 1){
          reject = 1  
        }
        else if(stepk$Acceptation == 1){
          reject = 0
        }
        else{
          k = k+1  
        }
      }
    }
  }
  return(list("reject"=reject, "keep"=keep, "stage"=k, "duration"=data_k$duration, 
              "nb_patients"=data_k$n_pat,"dur_incl"=data_k$dur_incl))
}


#########
### Simulation of n design with Magnusson and Turnbull for survival data
sim_trials_OS_MT = function(K_stages, N_subsets, f, l, u, ratio_Delta_star_d1, incl_rate, mean_cur_c, 
                            HR, nb_required, nmax_wait=+Inf, ordering, increasing_theta=FALSE, nsim=1000, seed=42){
  set.seed(seed) 
  prob_rejec = 0
  prob_accep = 0
  list_keep = list()
  pct_keep = c()
  pct_rejec_keep = c()
  pct_rejec_keep1 = c()
  pct_rejec_keep2 = c()
  rejec_stage = numeric(K_stages)
  accep_stage = numeric(K_stages)
  mean_duration = 0
  mean_pat = 0
  mean_dur_incl = 0
  dist_pat = c()
  dist_duration = c()
  dist_dur_incl = c()
  quant_pat = c()
  quant_duration = c()
  quant_dur_incl = c()
  
  for(isim in 1:nsim){
if(isim %% 1000 == 0){
print(isim)
}
    sMT = sim_one_OS_MT(K_stages, N_subsets, f, l, u, ratio_Delta_star_d1, incl_rate,  
                        mean_cur_c, HR, nb_required, nmax_wait, ordering, increasing_theta) 
    reject = sMT$reject
    st = sMT$stage

    if(reject == 1){
      prob_rejec = prob_rejec+1
      rejec_stage[st] = rejec_stage[st]+1
    }
    else{
      prob_accep = prob_accep+1 
      accep_stage[st] = accep_stage[st]+1
    }
    keep = sMT$keep
    inl = in_list(keep,list_keep)
    if(inl[[1]]==TRUE){
      pct_keep[inl[[2]]] = pct_keep[inl[[2]]]+1
      if(reject==1){
        pct_rejec_keep[inl[[2]]] = pct_rejec_keep[inl[[2]]]+1
        if(st==1){
          pct_rejec_keep1[inl[[2]]] = pct_rejec_keep1[inl[[2]]]+1
        }
        else if(st==2){
          pct_rejec_keep2[inl[[2]]] = pct_rejec_keep2[inl[[2]]]+1
        }
      }
    }
    else{
      list_keep = c(list_keep, list(keep))
      pct_keep = c(pct_keep, 1)
      if(reject==1){
        pct_rejec_keep = c(pct_rejec_keep, 1)
        if(st==1){
          pct_rejec_keep1 = c(pct_rejec_keep1, 1)
          pct_rejec_keep2 = c(pct_rejec_keep2, 0)
        }
        else if(st==2){
          pct_rejec_keep1 = c(pct_rejec_keep1, 0)
          pct_rejec_keep2 = c(pct_rejec_keep2, 1)
        }
      }
      else{
        pct_rejec_keep = c(pct_rejec_keep, 0)
        if(st==1){
          pct_rejec_keep1 = c(pct_rejec_keep1, 0)
          pct_rejec_keep2 = c(pct_rejec_keep2, 0)
        }
        else if(st==2){
          pct_rejec_keep1 = c(pct_rejec_keep1, 0)
          pct_rejec_keep2 = c(pct_rejec_keep2, 0)
        }
      }
    }
    mean_duration = mean_duration + sMT$duration
    mean_pat = mean_pat + sMT$nb_patients
    if(is.infinite(nmax_wait)==FALSE){
      mean_dur_incl = mean_dur_incl + sMT$dur_incl
      dist_pat = c(dist_pat, sMT$nb_patients)
      dist_duration = c(dist_duration, sMT$duration)
      dist_dur_incl = c(dist_dur_incl, sMT$dur_incl)
    }
  }
  prob_rejec = prob_rejec/nsim
  prob_accep = prob_accep/nsim
  pct_keep = pct_keep/nsim*100
  pct_rejec_keep = pct_rejec_keep/nsim*100
  pct_rejec_keep1 = pct_rejec_keep1/nsim*100
  pct_rejec_keep2 = pct_rejec_keep2/nsim*100
  rejec_stage = rejec_stage/nsim*100
  accep_stage = accep_stage/nsim*100
  mean_duration = mean_duration/nsim
  mean_pat = mean_pat/nsim
  mean_dur_incl = mean_dur_incl/nsim
  quant_pat = quantile(dist_pat, probs =c(0.5,0.75,1))
  quant_duration = quantile(dist_duration, probs =c(0.5,0.75,1))
  quant_dur_incl = quantile(dist_dur_incl, probs =c(0.5,0.75,1))
  
  return(list("prob_rejec"=prob_rejec, "prob_accep"=prob_accep, "list_keep"=list_keep, "pct_keep"=pct_keep,
              "pct_rejec_keep"=pct_rejec_keep, "pct_rejec_keep1"=pct_rejec_keep1, "pct_rejec_keep2"=pct_rejec_keep2,
              "rejec_stage"=rejec_stage, "accep_stage"=accep_stage, "mean_duration"=mean_duration, "mean_pat"=mean_pat,
              "mean_dur_incl"=mean_dur_incl, "dist_pat"=dist_pat, "dist_duration"=dist_duration, "dist_dur_incl"=dist_dur_incl,
              "quant_pat"=quant_pat, "quant_duration"=quant_duration, "quant_dur_incl"=quant_dur_incl))
}


#########
### Test for continuous or binary outcome 
test_BC = function(ind_trt_group_j, ind_con_group_j, outcome, type_outcome){
    n1 = length(ind_trt_group_j)
    n2 = length(ind_con_group_j)
    outcome_trt1 = outcome[ind_trt_group_j]
    outcome_trt2 = outcome[ind_con_group_j]
    mean1 = mean(outcome_trt1)
    mean2 = mean(outcome_trt2)
    if(n1 > 1 && n2 > 1){
      if(type_outcome=="binary"){
        if( mean1 == mean2 && (mean1 == 0 || mean1 == 1) ){
          Z_1j = NA
          I_1j = 1
        }
        else{
          prop_commune = (mean1*n1+mean2*n2)/(n1+n2)
          var_pool = prop_commune*(1-prop_commune)
          Z_1j = (mean1-mean2) / sqrt( var_pool * (1/n1 + 1/n2) )
          I_1j = (n1+n2) / (4*mean(c(outcome_trt1,outcome_trt2))*(1-mean(c(outcome_trt1,outcome_trt2))))
        }
      }
      else if(type_outcome=="continuous"){ 
        var1 = var(outcome_trt1)
        var2 = var(outcome_trt2)
        var_pool = ( (n1-1)*var1+(n2-1)*var2 ) / (n1+n2-2)
        Z_1j = (mean1-mean2) / sqrt(var_pool*(1/n1+1/n2))
        I_1j = (n1+n2) / 4
      }
    }
    else{
      Z_1j = NA
      I_1j = 1
    }
    return(c(Z_1j,I_1j))
}


#########
### Simulation of one design with Magnusson and Turnbull for binary or continuous data
sim_one_BC_MT = function(K_stages, N_subsets, f, l, u, ratio_Delta_star_d1, n_max, type_outcome, param_outcome, ordering, increasing_theta=FALSE){
  reject = NA
  outcome = c()
  trt = c()
  biom_group = c()
  n_pat = 0
  
  n_1 = ceiling(n_max / (1+sum(ratio_Delta_star_d1)))
  #n_step = c(n_1, n_1*ratio_Delta_star_d1)
  n_step = c(n_1, n_max-n_1)
  
  #Step 1
  while(n_pat < n_step[1]){
    n_pat = n_pat+1
    biom_group_cur = sample(1:length(f), 1, prob=f)
    biom_group = c(biom_group, biom_group_cur)
    trt_cur = rbinom(1,1,0.5)
    trt = c(trt, trt_cur)
    if(type_outcome=="binary"){
        outcome = c(outcome, rbinom(1,1,param_outcome[[1]][trt_cur+1,biom_group_cur]))
    }
    else if(type_outcome=="continuous"){
        outcome = c(outcome, rnorm(1, param_outcome[[1]][trt_cur+1,biom_group_cur], sqrt(param_outcome[[2]][trt_cur+1,biom_group_cur])))
    } 
  }
  Y_1j = numeric(N_subsets)
  I_1j = numeric(N_subsets)
  for(j in 1:N_subsets){
    ind_trt_group_j = which(trt==1 & biom_group==j)
    ind_con_group_j = which(trt==0 & biom_group==j)
    t_stat = test_BC(ind_trt_group_j, ind_con_group_j, outcome, type_outcome)
    I_1j[j] = t_stat[2]
    Y_1j[j] = t_stat[1] * sqrt(I_1j[j])
    if(is.na(Y_1j[j])){
      Y_1j[j] = -Inf
      I_1j[j] = 1
    }
  }
  k=1
  selection = magnusson_turnbull(0, NA, N_subsets, Y_1j, I_1j, l, u, ordering, increasing_theta)
  if(selection$Acceptation == 1){
    reject = 0
    keep = c()
  }
  else{
    keep = selection$Keep
    ind_group_1S = c()
    for(j in keep){
      ind_group_1S = c(ind_group_1S, which(biom_group==j))
    }
    t_statS = test_BC(intersect(ind_group_1S,which(trt==1)), intersect(ind_group_1S,which(trt==0)), outcome, type_outcome)
    I_1S = t_statS[2]
    Y_1S = t_statS[1] * sqrt(I_1S)
    if(is.na(Y_1S)){
      Y_1S = -Inf
      I_1S = 1
    }
    step1 = magnusson_turnbull(1, keep, N_subsets, Y_1S, I_1S, l, u, ordering, increasing_theta) 
    if(step1$Rejection == 1){
      reject = 1
    }  
    else{  
      k = 2
      f_keep = f[keep]
      f_S = sum(f_keep)
      while(k <= K_stages && is.na(reject)){ 
        while(n_pat < cumsum(n_step)[k]){
          n_pat = n_pat+1
          if(length(keep)>1){
            biom_group_cur = sample(keep, 1, prob=f_keep/f_S)
            biom_group = c(biom_group, biom_group_cur)
          }
          else{
            biom_group_cur = keep
            biom_group = c(biom_group, biom_group_cur)
          }
          trt_cur = rbinom(1,1,0.5)
          trt = c(trt, trt_cur)
          if(type_outcome=="binary"){
            outcome = c(outcome, rbinom(1,1,param_outcome[[1]][trt_cur+1,biom_group_cur]))
            
          }
          else if(type_outcome=="continuous"){
            outcome = c(outcome, rnorm(1, param_outcome[[1]][trt_cur+1,biom_group_cur], sqrt(param_outcome[[2]][trt_cur+1,biom_group_cur])))
          }
        }
        ind_group_S = c()
        for(j in keep){
          ind_group_S = c(ind_group_S, which(biom_group==j))
        }
        t_statkS = test_BC(intersect(ind_group_S,which(trt==1)), intersect(ind_group_S,which(trt==0)), outcome, type_outcome)
        I_kS = t_statkS[2]
        Y_kS = t_statkS[1] * sqrt(I_kS)
        stepk = magnusson_turnbull(k, keep, N_subsets, Y_kS, I_kS, l, u, ordering, increasing_theta)
        if(stepk$Rejection == 1){
          reject = 1  
        }
        else if(stepk$Acceptation == 1){
          reject = 0
        }
        else{
          k = k+1  
        }
      }
    }
  }
  return(list("reject"=reject, "keep"=keep, "stage"=k, "nb_patients"=length(trt)))
}


#########
### Simulation of n design with Magnusson and Turnbull for binary or continuous data
sim_trials_BC_MT = function(K_stages, N_subsets, f, l, u, ratio_Delta_star_d1, n_max, type_outcome, param_outcome, ordering, 
                            increasing_theta=FALSE, nsim=1000, seed=42){
  set.seed(seed) 
  prob_rejec = 0
  prob_accep = 0
  list_keep = list()
  pct_keep = c()
  pct_rejec_keep = c()
  pct_rejec_keep1 = c()
  pct_rejec_keep2 = c()
  rejec_stage = numeric(K_stages)
  accep_stage = numeric(K_stages)
  mean_pat = 0
  dist_pat = c()

  for(isim in 1:nsim){
if(isim %% 1000 == 0){
print(isim)
}
    sMT = sim_one_BC_MT(K_stages, N_subsets, f, l, u, ratio_Delta_star_d1, n_max, type_outcome, param_outcome, ordering, increasing_theta)
    reject = sMT$reject
    st = sMT$stage
    if(reject == 1){
      prob_rejec = prob_rejec+1
      rejec_stage[st] = rejec_stage[st]+1
    }
    else{
      prob_accep = prob_accep+1 
      accep_stage[st] = accep_stage[st]+1
    }
    keep = sMT$keep
    inl = in_list(keep,list_keep)
    if(inl[[1]]==TRUE){
      pct_keep[inl[[2]]] = pct_keep[inl[[2]]]+1
      if(reject==1){
        pct_rejec_keep[inl[[2]]] = pct_rejec_keep[inl[[2]]]+1
        if(st==1){
          pct_rejec_keep1[inl[[2]]] = pct_rejec_keep1[inl[[2]]]+1
        }
        else if(st==2){
          pct_rejec_keep2[inl[[2]]] = pct_rejec_keep2[inl[[2]]]+1
        }
      }
    }
    else{
      list_keep = c(list_keep, list(keep))
      pct_keep = c(pct_keep, 1)
      if(reject==1){
        pct_rejec_keep = c(pct_rejec_keep, 1)
        if(st==1){
          pct_rejec_keep1 = c(pct_rejec_keep1, 1)
          pct_rejec_keep2 = c(pct_rejec_keep2, 0)
        }
        else if(st==2){
          pct_rejec_keep1 = c(pct_rejec_keep1, 0)
          pct_rejec_keep2 = c(pct_rejec_keep2, 1)
        }
      }
      else{
        pct_rejec_keep = c(pct_rejec_keep, 0)
        if(st==1){
          pct_rejec_keep1 = c(pct_rejec_keep1, 0)
          pct_rejec_keep2 = c(pct_rejec_keep2, 0)
        }
        else if(st==2){
          pct_rejec_keep1 = c(pct_rejec_keep1, 0)
          pct_rejec_keep2 = c(pct_rejec_keep2, 0)
        }
      }
    }
    mean_pat = mean_pat + sMT$nb_patients
    dist_pat = c(dist_pat, sMT$nb_patients)
  }
  prob_rejec = prob_rejec/nsim
  prob_accep = prob_accep/nsim
  pct_keep = pct_keep/nsim*100
  pct_rejec_keep1 = pct_rejec_keep1/nsim*100
  pct_rejec_keep2 = pct_rejec_keep2/nsim*100
  pct_rejec_keep = pct_rejec_keep/nsim*100
  rejec_stage = rejec_stage/nsim*100
  accep_stage = accep_stage/nsim*100
  mean_pat = mean_pat/nsim
  
  return(list("prob_rejec"=prob_rejec, "prob_accep"=prob_accep, "list_keep"=list_keep, "pct_keep"=pct_keep,
              "pct_rejec_keep"=pct_rejec_keep, "pct_rejec_keep1"=pct_rejec_keep1, "pct_rejec_keep2"=pct_rejec_keep2,
              "rejec_stage"=rejec_stage, "accep_stage"=accep_stage, "mean_pat"=mean_pat, "dist_pat"=dist_pat))
}





#########
### General function to simulate n design with Magnusson and Turnbull for a given type of outcome
sim_magnusson_turnbull = function(K_stages, N_subsets, f, l, u, ratio_Delta_star_d1, type_outcome, param_outcome=NA,
                                  n_max=NA, incl_rate=NA, mean_cur_c=NA, HR=NA, 
                                  nb_required=NA, nmax_wait=+Inf, ordering, increasing_theta=FALSE, nsim=1000, seed=42){
  catch = catch_entries_MT(K_stages, N_subsets, f, l, u, ratio_Delta_star_d1, type_outcome, param_outcome, n_max, incl_rate, 
                        mean_cur_c, HR, nb_required, nmax_wait, ordering, increasing_theta, nsim, seed)
  K_stages = catch[[1]]
  N_subsets = catch[[2]]
  f = catch[[3]]
  l = catch[[4]]
  u = catch[[5]]
  ratio_Delta_star_d1 = catch[[6]]
  type_outcome = catch[[7]]
  param_outcome = catch[[8]]
  n_max = catch[[9]]
  incl_rate = catch[[10]]
  mean_cur_c = catch[[11]]
  HR = catch[[12]]
  nb_required = catch[[13]]
  ordering = catch[[14]]
  increasing_theta = catch[[15]]
  nsim = catch[[16]]
  seed = catch[[17]]
  nmax_wait = catch[[18]]
    
  if(type_outcome=="survival"){
    return(sim_trials_OS_MT(K_stages, N_subsets, f, l, u, ratio_Delta_star_d1, incl_rate, mean_cur_c, HR, 
                            nb_required, nmax_wait, ordering, increasing_theta, nsim, seed))  
  }
  else if(type_outcome=="binary" || type_outcome=="continuous"){
    return(sim_trials_BC_MT(K_stages, N_subsets, f, l, u, ratio_Delta_star_d1, n_max, type_outcome, param_outcome, 
                            ordering, increasing_theta, nsim, seed))
  }
}


