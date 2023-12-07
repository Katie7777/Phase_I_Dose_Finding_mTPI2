source("mTPI2_help_func_1.R")
source("mTPI2_help_func_2.R")


####################Getting Decision Boundary & Decision Table###############################
###########################################################################################
# The get_boundary() function is to generate 
# dose escalation/de-escalation rule including decision boundary and decision table                                                      #
#                                                                                         #
# The decision boundary will help physicians to assign different dose                     #
# based on the number of DLT observed                                                     #
#                                                                                         #
###########################################################################################

#                                        #############
#                                        # Arguments:#
#                                        #############
###########################################################################################
#     file_dir: location where the decision boundary & table  will be saved into 
#     
#     target:  the target toxicity rate
#
#     epsilon1: determine the lower bound of equivalence interval
#
#     epsilon2: determine the higher bound of equivalence interval
#
#     ncohort:   the total number of times to add patients
#
#     cohortsize:   number of patients in each cohort
#
#			exclusive_certainty:  cutoff to eliminate the overly toxic dose for safety monitoring
#
#     exclusive_certainty: The cutoff to eliminate the overly toxic dose for safety
#
#     prior_alpha: alpha of prior beta distribution for dose-assignment
#
#     prior_beta: beta of prior beta distribution for dose-assignment
#
#     max_pts_per_dose: maximum patients assigned to each dose 
############################################################################################

#parameters:
file_dir = "decision_rule.csv" 
target = 0.3
epsilon1 = 0.05
epsilon2 = 0.05
ncohort = 12
cohortSize = 3
exclusive_certainty = 0.95
prior_alpha = 1
prior_beta=1
max_pts_per_dose = 12

get_boundary(file_dir=file_dir, target = target, epsilon1 = epsilon1, epsilon2 = epsilon2, 
                         ncohort = ncohort, cohortSize = cohortSize, 
                         exclusive_certainty = exclusive_certainty, 
                         prior_alpha = prior_alpha, prior_beta=prior_beta,
                         max_pts_per_dose = max_pts_per_dose)




######################################## Run Simulation #########################################

###########################################################################################
# The mTPI2_summary() function is to run the simulation "ntrial" times                        #
#                                                                                         #
#  Simulation includes the dose assignment based on mTPI dose assignment rule             #
#  and the MTD finding based on dose-finding rule.                                        #
#                                                                                         #
#                                                                                         #
###########################################################################################

#
#                                        #############
#                                        # Arguments:#
#                                        #############
###########################################################################################
#     scenarios: data frame that contains all different scenarios of true DLT vector 
#
#     file_dir: location where the simulation result will be saved into 
#
#     ntrial: number of simulated trials
#
#     true_toxic:  the true toxicity rate for each dose
#
#     startdose:  the starting dose of the trial
#
#     ncohort:   the total number of cohorts
#
#     cohortsize:   cohort size
#
#     target:  the target toxicity rate
#
#     epsilon1: determine the lower bound of equivalence interval
#
#     epsilon2: determine the higher bound of equivalence interval
#
#     decison_table_prior_alpha: prior alpha used to make dose-assignment
#
#     decison_table_prior_beta: prior beta used to make dose-assignment
#
#     mtd_prior_alpha: prior alpha used in dose-finding
#
#     mtd_prior_beta: prior beta used in dose-finding
#
#	exclusive_certainty:  cutoff to eliminate the overly toxic dose for safety monitoring
#
#     max_pts_per_dose: the max number of patients assigned to each dose 
#
#     follow_up_time : follow up time for non DLT patients
#
#     monthly_pts_enroll: number of patients enrolled per month
############################################################################################

#parameters:

scenarios = data.frame(s1 = c(0.3,0.46,0.50,0.54,0.58),
                       s2 = c(0.16,0.3,0.47,0.54,0.60),
                       s3 = c(0.04,0.15,0.3,0.48,0.68),
                       s4 = c(0.02,0.07,0.12,0.3,0.45),
                       s5 = c(0.02,0.06,0.1,0.13,0.3))
file_dir="mTPI2_simulation_result.csv"
ntrial = 1000
startdose = 1
ncohor = 12
cohortsize = 3
target = 0.30
epsilon1 = 0.05
epsilon2 = 0.05
decison_table_prior_alpha = 1
decison_table_prior_beta = 1
mtd_prior_alpha =10^(-5)
mtd_prior_beta = 10^(-5)
exclusive_certainty = 0.95
max_pts_per_dose = 18
follow_up_time = 30
monthly_pts_enroll = 3

mTPI2_summary(scenarios=scenarios,
              file_dir= file_dir,ntrial= ntrial, startdose = startdose, 
              ncohor = ncohor, 
              cohortsize = cohortsize, target = target, epsilon1 = epsilon1, epsilon2 = epsilon2, 
              decison_table_prior_alpha = decison_table_prior_alpha,
              decison_table_prior_beta = decison_table_prior_beta, 
              mtd_prior_alpha = mtd_prior_alpha, 
              mtd_prior_beta = mtd_prior_beta,
              exclusive_certainty = exclusive_certainty, 
              max_pts_per_dose = max_pts_per_dose, 
              follow_up_time = follow_up_time,
              monthly_pts_enroll=monthly_pts_enroll)
















