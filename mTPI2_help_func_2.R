setwd("C:/Users/katie.wang/OneDrive")
source("mTPI2_help_func_1.R")
library(pillar)


####################OUTPUT DECISION BOUNDARY AND DECISION TABLE###########################
get_boundary <- function(file_dir, target , epsilon1, epsilon2 , 
                         ncohort , cohortSize , 
                         exclusive_certainty , 
                         prior_alpha , prior_beta,
                         max_pts_per_dose ){
  
  
  #get decision boundary 
  decision_boundary = get_decision_boundary(target = target, epsilon1 = epsilon1, epsilon2 = epsilon2, 
                                            ncohort = ncohort, cohortSize = cohortSize, 
                                            exclusive_certainty = exclusive_certainty, 
                                            prior_alpha = prior_alpha, prior_beta=prior_beta,
                                            max_pts_per_dose = max_pts_per_dose )
  
  
  
  decision_table = get_decision_table(target = target, epsilon1 = epsilon1, epsilon2 = epsilon2, 
                                      ncohort = ncohort, cohortSize = cohortSize, 
                                      exclusive_certainty = exclusive_certainty, 
                                      prior_alpha = prior_alpha, prior_beta=prior_beta,
                                      max_pts_per_dose = max_pts_per_dose )
  
  parameters = c( "ncohort", "cohortsize",
                 "target", "epsilon1", "epsilon2","exclusive_certainty",
                 "max_pts_per_dose", "prior_alpha", "prior_beta")
  meaning = c("the total number of times to add patients",
              "number of patients in each cohort",
              "the target toxicity rate",
              "determine the lower bound of equivalence interval",
              "determine the higher bound of equivalence interval",
              "The cutoff to eliminate the overly toxic dose for safety",
              "the max number of patients assigned to each dose",
              "alpha of prior beta distribution for dose-assignment",
              "beta of prior beta distribution for dose-assignment")
  values = c(ncohort, cohortSize,target, epsilon1, epsilon2,
             exclusive_certainty,max_pts_per_dose, prior_alpha, prior_beta )
  df = data.frame(parameters, meaning, values)
  
  #write out the decision table and decision boundary
  write.table(df, file_dir, na = "NA", append = TRUE, 
              col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
  write.table('\n', file_dir, append = TRUE, 
              col.names = FALSE, row.names = FALSE, sep = ",", quote = FALSE)
  
  
  write.table(decision_boundary,file_dir , na = "NA", append = TRUE, 
              col.names = NA, row.names = TRUE, sep = ",", quote = FALSE)
  
  write.table('\n', file_dir, append = TRUE, 
              col.names = FALSE, row.names = FALSE, sep = ",", quote = FALSE)
  
  write.table(decision_table,file_dir , na = "NA", append = TRUE, 
              col.names = NA, row.names = TRUE, sep = ",", quote = FALSE)
  
  write.table('\n', file_dir, append = TRUE, 
              col.names = FALSE, row.names = FALSE, sep = ",", quote = FALSE)

}



#############################GET DECISION BOUNDARY#########################

#generate escalation/de-escalation boundary
get_decision_boundary <-function(target = 0.3, epsilon1 = 0.05, epsilon2 = 0.05, 
                              ncohort = 4, cohortSize = 3, exclusive_certainty = 0.95,
                              prior_alpha = 1, 
                              prior_beta=1, max_pts_per_dose = NA ){
  npts = ncohort * cohortSize
  result = c()
  
  if(!is.na(max_pts_per_dose)){
    npts = max_pts_per_dose
  }
  
  for(i in 1: npts){
    result_i = rep(0, npts+1)
    boundary_i = c()
    for(ndlt in 0: i ){
      decision = get_decision(target = target, epsilon1 = epsilon1, epsilon2 = epsilon2,
                              alpha =prior_beta + ndlt , beta = prior_beta + i - ndlt)
      
      if(decision == 'D' & i >=3){
        check_du =check_if_exclude_dose_d(
          target = target,
          exclusive_certainty = exclusive_certainty,
          ndlt,
          i,
          prior_alpha = prior_alpha,
          prior_beta = prior_beta) 
        if (check_du == 1){
          decision = "DU"
        }
      }
      result_i[ndlt+1] = decision
    }

    #get the index of the last occurrence of "E"
    res_e = tapply(seq_along(result_i),result_i,max,simplify = FALSE)
    e = res_e$E -1
    boundary_i = c(boundary_i,e)
    
    #get the index of the first occurence of "D" and "DE"
    res_d_du= tapply(seq_along(result_i),result_i,min,simplify = FALSE)
    d = res_d_du$D -1
    du = res_d_du$DU
    if(is.null(du)){
      du = "NA"
    }else{
      du = du -1
    }
    boundary_i = c(boundary_i,d)
    boundary_i = c(boundary_i,du)
    
    #combine all boundary together
    result = cbind(result, i = boundary_i)

  }
  
  colnames(result)= 1:npts
  #convert vector to data frame
  result = as.data.frame(result)
  colnames =1:npts
  result=rbind(colnames, result)
  row.names(result) = c("Number of patients treated","Escalate if # of DLT <=", "Deescalate if # of DLT >=", "Eliminate if # of DLT >=")
  return(result)
}


####################GET DECISION TABLE#########################

#generate dose assignment table
get_decision_table <-function(target = 0.3, epsilon1 = 0.05, epsilon2 = 0.05, 
                          ncohort = 4, cohortSize = 3, exclusive_certainty = 0.95,
                          prior_alpha = 1, 
                          prior_beta=1, max_pts_per_dose = NA ){
  npts = ncohort * cohortSize
  result = c()
  
  if(!is.na(max_pts_per_dose)){
    npts = max_pts_per_dose
  }
  
  for(i in 1: npts){
    result_i = rep(0, npts+1)
    
    for(ndlt in 0: i ){
      decision = get_decision(target = target, epsilon1 = epsilon1, epsilon2 = epsilon2,
                              alpha =prior_beta + ndlt , beta = prior_beta + i - ndlt)
      #print(c("decision: ", decision))
      
      if(decision == 'D' & i >=3){
        check_du =check_if_exclude_dose_d(
          target = target,
          exclusive_certainty = exclusive_certainty,
          ndlt,
          i,
          prior_alpha = prior_alpha,
          prior_beta = prior_beta) 
        if (check_du == 1){
          decision = "DU"
        }
      }
      result_i[ndlt+1] = decision
    }
    #print(c("i: ", i))
    #print(c("dec under i pts: ",result_i))
    result = cbind(result, i = result_i)
    
  }
  #convert vector to data frame
  result = as.data.frame(result)
  result[result == 0] = ""
  colnames(result) = 1:npts
  rownames(result) = 0:npts
  return(result)
}



######################RUN SIMULATION#######################

#run simulations
mTPI2_sim <-function(file_dir, ntrial=1000, true_toxic, startdose = 1, ncohor = 12, cohortsize = 3, 
                     target = 0.3, epsilon1 = 0.05, epsilon2 = 0.05, decison_table_prior_alpha = 1, 
                     decison_table_prior_beta = 1, mtd_prior_alpha = 0.005, mtd_prior_beta = 0.005, 
                     exclusive_certainty = 0.95, max_pts_per_dose = NA, follow_up_time = 30, 
                     monthly_pts_enroll = 3){
  set.seed(6)
  ndose = length(true_toxic)
  selected_dose_level = rep(0, ntrial)
  Y=matrix(rep(0, ndose*ntrial), ncol=ndose); # store toxicity outcome
  N=matrix(rep(0, ndose*ntrial), ncol=ndose); # store the number of patients
  trial_duration = rep(0, ntrial)
  check_toxicity = 0
  pat_in_mtd = rep(0, ntrial)
 
   for (trail in 1:ntrial ){
    # run one trail
    one_trail = one_sim_trail(p.true = true_toxic,
               startdose = startdose,
               num_cohort = ncohor,
               cohort_size = cohortsize,
               target = target,
               epsilon1 = epsilon1,
               epsilon2 = epsilon2,
               prior_alpha = decison_table_prior_alpha,
               prior_beta = decison_table_prior_beta,
               exclusive_certainty = exclusive_certainty,
               max_pts_per_dose=max_pts_per_dose,
               follow_up_time = follow_up_time,
               monthly_pts_enroll = monthly_pts_enroll)
   
    dlt_obs = as.numeric(unlist(one_trail[1]))
    pts_obs = as.numeric(unlist(one_trail[2]))
    elimi = as.numeric(unlist(one_trail[3]))
    trial_duration[trail] = as.numeric(one_trail[4])
    

    #select mtd for this trail
    dselect = select_mtd(num_dose = ndose,
               dlt_obs = dlt_obs,
               pts_obs=pts_obs,
               elimi = elimi,
               target = target,
               prior_alpha = mtd_prior_alpha,# different than prior alpha in decision making
               prior_beta = mtd_prior_beta) 
    #save trail info into matrix
    selected_dose_level[trail] = dselect
    Y[trail, ] = dlt_obs
    N[trail, ] = pts_obs
    

    if(dselect != 0){
        pat_in_mtd[trail] = pts_obs[dselect]
    }
    
    
    if(sum(dlt_obs) > ncohor*cohortsize*target){
          check_toxicity = check_toxicity+1
    }
    

  }
  pats_non_seq = ncohor*cohortsize/ndose
  risk_poor_allocation = sum(pat_in_mtd < pats_non_seq, na.rm=T )/ntrial
  risk_high_toxicity = check_toxicity/ntrial
  #generate output
  selpercent=rep(0, ndose)
  for(i in 1:ndose) {
    selpercent[i]=sum(selected_dose_level==i)/ntrial*100
  }
  perct_pts_per_dose = t(apply(N,1, function(x) x/sum(x)))
  cat("selection percentage at each dose level (%):\n");
  cat(formatC(selpercent, digits=1, format="f"), sep="  ", "\n");
  cat("average number of patients treated at each dose level:\n");
  cat(formatC(apply(N,2,mean), digits=1, format="f"), sep ="  ", "\n");
  cat("true DLT rate:", formatC(true_toxic, digits=2, format="f"), "\n");
  cat("average DLT percentage per dose:", formatC(apply(Y/N,2,mean,na.rm=TRUE), digits=2, format="f"), "\n");
  cat("average percentage of patients per dose:", formatC(apply(perct_pts_per_dose,2,mean), digits=2, format="f"), "\n");
  cat("average number of DLT:", formatC(sum(Y)/ntrial, digits=1, format="f"), "\n");
  cat("average number of patients:", formatC(sum(N)/ntrial, digits=1, format="f"), "\n");
  cat("average trial_duration(month):", formatC(sum(trial_duration)/ntrial/30.4375, digits=1, format="f"), "\n");
  cat("25% quantile trial_duration(month):", formatC(quantile(trial_duration, prob = 0.25)/30.4375, digits=2, format="f"), "\n");
  cat("75% quantile trial_duration(month):", formatC(quantile(trial_duration, prob = 0.75)/30.4375, digits=2, format="f"), "\n");
  cat("risk of poor allocation:", formatC(risk_poor_allocation, digits=2, format="f"), "\n");
  cat("risk of high toxicity:", formatC(risk_high_toxicity, digits=2, format="f"), "\n");
  
  #create a dataframe
  df_title = data.frame(information_name = c("true DLT rate: ",
                                          "selection percentage at each dose level (%): ", 
                                         "average number of patients treated at each dose level: ",
                                         
                                         "average DLT percentage per dose: ",
                                         "average percentage of patients per dose: ",
                                         "average number of DLT: ",
                                         "average number of patients: ",
                                         "average trial_duration(month): ",
                                         "25% quantile trial_duration(month): ",
                                         "75% quantile trial_duration(month): ",
                                         "risk of poor allocation: ",
                                         "risk of high toxicity: "
                                         
                                         ))
  df_content = data.frame(content =
                          c(paste(formatC(true_toxic, digits=2, format="f"),collapse = ","),
                        paste(formatC(selpercent, digits=1, format="f"),collapse = ","),
                         paste(formatC(apply(N,2,mean), digits=1, format="f"),collapse = ","),
                          
                          paste(formatC(apply(Y/N,2,mean,na.rm=TRUE), digits=2, format="f"),collapse = ","),
                          paste(formatC(apply(perct_pts_per_dose,2,mean), digits=2, format="f"),collapse = ","),
                          paste(formatC(sum(Y)/ntrial, digits=1, format="f"),collapse = ","),
                          paste(formatC(sum(N)/ntrial, digits=1, format="f"),collapse = ","),
                          paste(formatC(sum(trial_duration)/ntrial/30.4375, digits=1, format="f"),collapse = ","),
                          paste(formatC(quantile(trial_duration, prob = 0.25)/30.4375, digits=2, format="f"),collapse = ","),
                          paste(formatC(quantile(trial_duration, prob = 0.75)/30.4375, digits=2, format="f"),collapse = ","),
                          paste(formatC(risk_poor_allocation, digits=2, format="f"),collapse = ","),
                          paste(formatC(risk_high_toxicity, digits=2, format="f"),collapse = ",")
                        ))
                          
  

  result = cbind(df_title, df_content)

  write.table(result, file_dir, na = "NA", append = TRUE, 
              col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
  write.table('\n', file_dir, append = TRUE, 
              col.names = FALSE, row.names = FALSE, sep = ",", quote = FALSE)
 
}

########### Run Simulation with different scenario######################

#scenarios: df contains scenario1, scenario2...
mTPI2_summary <- function(scenarios, file_dir, 
                          ntrial=1000, startdose = 1, ncohor = 12, 
                          cohortsize = 3, target = 0.3, epsilon1 = 0.05, 
                          epsilon2 = 0.05, decison_table_prior_alpha = 1, 
                          decison_table_prior_beta = 1, mtd_prior_alpha = 0.005, 
                          mtd_prior_beta = 0.005, 
                          exclusive_certainty = 0.95, 
                          max_pts_per_dose = NA, follow_up_time = 30, 
                          monthly_pts_enroll = 3){
      
      parameters = c("ntrial", "startdose", "ncohor", "cohortsize",
                     "target", "epsilon1", "epsilon2","exclusive_certainty",
                     "max_pts_per_dose", "follow_up_time", "monthly_pts_enroll")
      
      meaning = c("number of simulated trials",
                  "the starting dose of the trial",
                  "the total number of times to add patients",
                  "number of patients in each cohort",
                  "the target toxicity rate",
                  "determine the lower bound of equivalence interval",
                  "determine the higher bound of equivalence interval",
                  "The cutoff to eliminate the overly toxic dose for safety",
                  "the max number of patients assigned to each dose",
                  "follow up time for non DLT patients",
                  "number of patients enrolled per month")
      values = c(ntrial, startdose, ncohor, cohortsize,target, epsilon1, epsilon2,
                 exclusive_certainty,max_pts_per_dose, follow_up_time, monthly_pts_enroll )
      df = data.frame(parameters, meaning, values)
      
      write.table(df, file_dir, na = "NA", append = TRUE, 
                  col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
      write.table('\n', file_dir, append = TRUE, 
                  col.names = FALSE, row.names = FALSE, sep = ",", quote = FALSE)
      
      #scenarios: 
      for(i in 1: ncol(scenarios)){
            true_toxic = scenarios[, i]
            mTPI2_sim(file_dir = file_dir,ntrial= ntrial, true_toxic= true_toxic, startdose = startdose, ncohor = ncohor, 
                      cohortsize = cohortsize, target = target, epsilon1 = epsilon1, epsilon2 = epsilon2, 
                      decison_table_prior_alpha = decison_table_prior_alpha,
                      decison_table_prior_beta = decison_table_prior_beta, mtd_prior_alpha = mtd_prior_alpha, 
                      mtd_prior_beta = mtd_prior_beta,exclusive_certainty = exclusive_certainty, 
                      max_pts_per_dose=max_pts_per_dose,
                      follow_up_time = follow_up_time, 
                      monthly_pts_enroll=monthly_pts_enroll)
      }
      
}

