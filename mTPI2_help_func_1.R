library('Iso')



######################################################################
#helper functions:

#########################GET MAX UPM########################

#get max UPM in HI interval
max_upm_HI <- function(HI_seq, Equal_length, alpha, beta) {
  integrand <- function(x) {
    dbeta(x, shape1 = alpha, shape2 = beta)
  }
  max_UPM = 0
  len = length(HI_seq)
  for (i in 1:len) {
    lower = HI_seq[i]
    if (i == len) {
      upper = 1
    } else{
      upper = HI_seq[i + 1]
    }
    upm_i = integrate(integrand, lower = lower, upper = upper)$value / Equal_length
    if (upm_i > max_UPM) {
      max_UPM = upm_i
    }
  }
  return(max_UPM)
}


#get max UPM in LI interval
max_upm_LI <- function(LI_seq, Equal_length, alpha, beta) {
  integrand <- function(x) {
    dbeta(x, shape1 = alpha, shape2 = beta)
  }
  max_UPM = 0
  len = length(LI_seq)
  for (i in 1:len) {
    upper = LI_seq[i]
    if (i == len) {
      lower = 0
    } else{
      lower = LI_seq[i + 1]
    }
    upm_i = integrate(integrand, lower = lower, upper = upper)$value / Equal_length
    if (upm_i > max_UPM) {
      max_UPM = upm_i
    }
  }
  return(max_UPM)
}


##############GET DECISION###############################
#dose-assignment decision
get_decision <-
  function(target = 0.3,
           epsilon1 = 0.05,
           epsilon2 = 0.05,
           alpha = 3,
           beta = 3) {
    # calculate three interval based on target, epsilon1, epsilon2
    LI_upper = target - epsilon1
    Equal_length = epsilon1 + epsilon2
    HI_lower = target + epsilon2
    
    #subintervals split points in LI_interval and HI_interval
    LI_seq = seq(LI_upper, 0, by = -Equal_length)
    HI_seq = seq(HI_lower, 1, by = Equal_length)
    
    
    #calculate posterior probability of each interval
    integrand <- function(x) {
      dbeta(x, shape1 = alpha, shape2 = beta)
    }
    
    #find th max UPM in LI, EI, and HI
    maxUPM_LI = max_upm_LI(LI_seq, Equal_length, alpha, beta)
    
    upm_EI = integrate(integrand,
                       lower = (target - epsilon1),
                       upper = (target + epsilon2))$value / Equal_length
 
    maxUPM_HI = max_upm_HI(HI_seq, Equal_length, alpha, beta)

    #find max upm among all subintervals
    max = max(maxUPM_LI, upm_EI, maxUPM_HI)
    if (max == maxUPM_LI) {
      return("E")
    } else if (max == upm_EI) {
      return("S")
    } else{
      return("D")
    }
  }



##########################################################################
#check if current dose d need to be excluded (check DU)
check_if_exclude_dose_d <-
  function(
    target = 0.3,
    exclusive_certainty = 0.95,
    DLT_d,
    pts_d,
    prior_alpha = 1,
    prior_beta = 1) {
    
    #check if the probability of pi> pt is greater than 95%(exclusive_certainty)
    if (1 - pbeta(target, DLT_d + prior_alpha, pts_d - DLT_d + prior_beta) >
        exclusive_certainty) {
      return(1)
    }
    
    return(0)
  }

#########################################################
#select mtd among all doses
select_mtd <-
  function(num_dose = 5,
           dlt_obs,
           pts_obs,
           elimi,
           target = 0.3,
           prior_alpha = 0.005,
           prior_beta = 0.005) {
    mtd = 0
    beta_posterior_dst_mean = c()
    phat_var = c()
    nadmis = 0
    #calculate mean of each beta posterior distribution
    for (i in 1:num_dose) {
      #check if dose one is too toxic (or Ti == 0)
      if(i ==1 & elimi[1] == 1){
        return(mtd)
      }
      #check if dose 2 and above is too toxic (or Ti == 0)
      if (elimi[i] == 0) {
        alpha_i = prior_alpha + dlt_obs[i]
        beta_i = prior_beta + pts_obs[i] - dlt_obs[i]
        phat = alpha_i / (alpha_i + beta_i)
        beta_posterior_dst_mean = c(beta_posterior_dst_mean, phat)
        var_i = (alpha_i * beta_i) / ((alpha_i + beta_i) ^ 2 * (alpha_i + beta_i + 1))
        phat_var = c(phat_var, var_i)
      }
      
    }

    #using PAVA to transform phat s.t. phat_star increase with the dose level
    trans_phat = pava(beta_posterior_dst_mean, w = 1 / phat_var)
    Pt = rep(target, length(beta_posterior_dst_mean))
    
    #check the difference between trans_hat and target
    diff = abs(trans_phat - Pt)
    min_diff = min(diff)
    index_of_min = which(diff == min_diff)
    
    #find the index of smallest diff
    len = length(index_of_min)

    #if only one smallest number
    if (len == 1) {
      mtd = index_of_min[1]
      #if tow or more doses tie for the smallest differences
    } else if (min_diff < target) {
      mtd = index_of_min[len]
    } else if (min_diff >= target) {
      mtd = index_of_min[1]
    }
    return(mtd)
    
  }

######################CALCULATE COHORT DURATION#############################

cohort_duration <- function(p.true,
                            cohort_size = 3,
                            dlt_obs= 0,
                            dose_level = 0,
                            follow_up_time = 30, 
                            monthly_pts_enroll = 3){
  #generate start date for new cohort pts
  start_date = c()
  
  #if number of patients enrolled per month smaller than cohort size
  if(monthly_pts_enroll < cohort_size){
    num_follow_up_times_needed = ceiling(cohort_size/monthly_pts_enroll)
    time_period_lower_boundry = 1
    time_period_higher_boundry = follow_up_time
    
    #allocate all cohort patients to different month
    for(i in 1: num_follow_up_times_needed){
      if(i == num_follow_up_times_needed){
        num_pts_left = cohort_size - (i-1)*monthly_pts_enroll
        #generate start date for pts in last follow_up_time
        startdate = floor(runif(num_pts_left, 
                                min=time_period_lower_boundry, max=time_period_higher_boundry))
      }
      else{
        #generate start date for pts in each follow_up_time
        startdate = floor(runif(monthly_pts_enroll, 
                                min=time_period_lower_boundry, max=time_period_higher_boundry))
      }
      start_date = c(start_date, startdate)
      time_period_lower_boundry = time_period_higher_boundry+1
      time_period_higher_boundry = time_period_higher_boundry + follow_up_time
      
    }
    #if number of patients enrolled per month larger or equal cohort size
  }else{
    # we can only add cohort_size pts for current dose
    start_date = floor(runif(cohort_size, min=1, max=follow_up_time))
  }
  
  #generate xd(number of toxicity events) for dose d
  dlt_pts = rbinom(1, cohort_size, p.true[dose_level])
  dlt_obs[dose_level] =  dlt_obs[dose_level] + dlt_pts
  
  #generate pts id for dlt
  dlt_pts_id = sample(1:cohort_size ,dlt_pts , replace=F)
  
  #generate end_date and calculate current cohort duration
  end_date = rep(0, cohort_size)
  new_cohort_duration = 0
  
  #if no DLT in current cohort
  if(length((dlt_pts_id))==0){
    end_date = start_date + follow_up_time
    
  }
  #if >=1 DLT in current cohort
  else{
    for (i in 1: cohort_size){
      if(any(i == dlt_pts_id)){
        #generate patients follow up time within follow_up_time
        DLT_follow_up_time = floor(runif(1, min=1, max=follow_up_time))
        end_date[i] = start_date[i] + DLT_follow_up_time
      }else
        end_date[i] = start_date[i] + follow_up_time
    }
    
  }
  
  #calculate the trail duration for current cohort
  new_cohort_duration = max(end_date) - start_date[which.min(start_date)]
  res = list(new_cohort_duration, dlt_obs)
  return(res)
}

########################ONE SIMULATION############################
# one simulation
one_sim_trail <-
  function(p.true,
           startdose = 1,
           num_cohort = 10,
           cohort_size = 3,
           target = 0.3,
           epsilon1 = 0.05,
           epsilon2 = 0.05,
           prior_alpha = 1,
           prior_beta = 1,
           exclusive_certainty = 0.95, 
           max_pts_per_dose = NA, 
           follow_up_time = 30, monthly_pts_enroll = 3) {
    
    ndose = length(p.true)
    dlt_obs = rep(0, ndose)
    pts_obs = rep(0, ndose)
    d = startdose
    elimi = rep(0, ndose)
    trial_duration = 0
    # initialize d_next
    d_next = d = startdose
    
    # for each cohort we make a decision of "D","S" or "E"
    for (i in 1:num_cohort) {
      
      #assign d_next to d so current dose is d
      d = d_next

      #check if current dose reaches max_pts and the decision is stay
      if(!is.na(max_pts_per_dose) & pts_obs[d] >= max_pts_per_dose){
          break
        }
      
      
      #get nd for dose d
      pts_obs[d] = pts_obs[d] + cohort_size
      
      #calculate cohort trial duration
      res = cohort_duration(p.true = p.true,
                                  cohort_size = cohort_size,
                                  dlt_obs=dlt_obs,
                                  dose_level = d,
                                  follow_up_time = follow_up_time, 
                                 monthly_pts_enroll = monthly_pts_enroll)
      
      new_cohort_duration = as.numeric(res[1])
      dlt_obs = as.numeric(unlist(res[2]))
      trial_duration = trial_duration + new_cohort_duration


      #update alpha and beta to get posterior distribution of Pd
      alpha_d = prior_alpha + dlt_obs[d]
      beta_d = prior_beta + pts_obs[d] - dlt_obs[d]
     
      #based on posterior distribution of Pd make decision
      decision_d = get_decision(
        target = target,
        epsilon1 = epsilon1,
        epsilon2 = epsilon2,
        alpha = alpha_d,
        beta = beta_d
      )
      
      #check each decision
      #if we need escalate the dose
      if (decision_d == "E") {
        if(d < ndose){
          #check if d+1 is too toxic
          if (elimi[d + 1] == 1) {
            d_next = d
          }
          else{
            d_next = d + 1
          }
        }else if( d== ndose){
          d_next = d
        }
      }
      #if we need to de-escalate the dose
      else if (decision_d == "D") {
        #check du for current dose when >=3 pts at does d
        if (d > 1 & pts_obs[d] >= 3) {
          #check if current d is too toxic (DU)
          check_du = check_if_exclude_dose_d(
            target = target,
            exclusive_certainty = exclusive_certainty,
            dlt_obs[d],
            pts_obs[d],
            prior_alpha = prior_alpha,
            prior_beta = prior_alpha
          )
          #if it's too toxic, change all the value after/include this dose to 1
          if (check_du == 1) {
            elimi[d:ndose] = 1
          }
          d_next = d - 1
          #check if dose_1 is D or DU
        } else if (d == 1) {
          #check du for cohort_size >1
          if (cohort_size > 1) {
            check_dose1 = check_if_exclude_dose_d(
              target = target,
              exclusive_certainty = exclusive_certainty,
              dlt_obs[1],
              pts_obs[1],
              prior_alpha = prior_alpha,
              prior_beta = prior_beta
            )
            #if dose one is too toxic, end the trail
            if (check_dose1 == 1) {
              elimi[d:ndose] = 1
              break
            } else{
              #if not DU, keep dose level at dose one
              d_next = d
            }
            
          }
          else if (cohort_size == 1) {
            #if cohort_size ==1, dose one patients num <3
            #keep dose level at dose one
            if (pts_obs[1] < 3) {
              d_next = d
            } else if (pts_obs[1] >= 3) {
              #check du if patients at dose one >=3
              if (cohort_size > 1) {
                check_dose1 = check_if_exclude_dose_d(
                  target = target,
                  exclusive_certainty = exclusive_certainty,
                  dlt_obs[1],
                  pts_obs[1],
                  prior_alpha = prior_alpha,
                  prior_beta = prior_beta
                )
                if (check_dose1 == 1) {
                  break
                } else{
                  #continue add patients to dose one
                  d_next = d
                }
              }
            }
          }
        }
      } else{
        d_next = d
      }
      
    }
    r = list(dlt_obs, pts_obs, elimi, trial_duration)
    return(r)
    
  }

