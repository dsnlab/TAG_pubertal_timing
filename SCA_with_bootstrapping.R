##########################################################################################################
### Run SCAs and permutations ####
### This scripts takes a dataset
### and runs SCAs using regression (not structural equation modelling)
### It consistes of multiple steps 
### 1. Loading libraries, controlling analyses and loading data
### 2. Running SCAs via three main functions get_data, get_model, get_coef
### 3. Make null models for bootstrapping
### 4. Run bootstrap models
### 5. Analyse bootstrap models 
### The output includes specification curve analyses and one csv file of the bootstrap tests
###
### Adapted from Amy Orben's scripts for the screens, teens, and well-being study https://github.com/OrbenAmy/PS_2019
##########################################################################################################

#########
# Issues:
# -Two participants still have no wave 2 data after imputation because they had no observed values at wave 2
# -the script pulls eta-squared effect sizes but only works for lm, not logistic models. These are not really necessary though


##########################################################################################################
# 1. Loading libraries, controlling analyses and loading data ############################################
##########################################################################################################

#####################################################################################
# a) Load libraries ####
#####################################################################################
library("tidyverse")
library("stringr")
library("foreign")
library("lavaan")
library("heplots")
set.seed(98)
cas_dir="Y:/dsnlab/TAG/"
samplesize <- 174

#####################################################################################
# b) Turn on analyses ####
# 1 = yes, 0 = no
#####################################################################################
run_sca <- 1
run_null <- 1
run_boot <- 1

#####################################################################################
# c) Determine how many bootstraps to run
#####################################################################################
bootstraps <- 500

#####################################################################################
# d) Load data #### changed to direct to our imputed and non-imputed datasets
#####################################################################################
data_all <- read.csv(file=paste0(cas_dir,"projects/W1_W2_pubertal_timing/Final_timing_int.csv"))

##########################################################################################################
# 2. Running SCAs via three main functions, get_data, get_model and get_coef #############################
##########################################################################################################

#####################################################################################
#### Function: get_data
#### This function takes the specification and makes the apporiate dataset
####
#### Input: 
#### results_frame = a data frame showing specifications on each line
#### dataset = a character saying what dataset the specifications are for
#### data = the data of interest
####
#### Output: A list with the folliwng items
#### results_frame = retain results_frame to use in next function
#### dataset = retain dataset to use in next function
#### data = dataset of interest, now with new dv and iv from specifications
####        furthermore for the bootstraps sample the data is bootstrapped
####
#### Method: It goes to the row of interest (i) of the results frame.
#### There it reads what the relevant iv and dv variables are and extracts these
#### from the dataset of interest.
#####################################################################################
get_data <- function(results_frame, dataset, data, boot) {

  data$iv <- data[,as.character(results_frame$predictor[i])]
  data$dv <- data[,as.character(results_frame$outcome[i])]
  data$wave1_control <- data[,as.character(results_frame$wave1_control[i])]
  data$ctq_control <- data[,as.character(results_frame$ctq_control[i])]
  
  # Do bootstrapping
  if(boot == TRUE){
    if (b <= bootstraps){
      k=sample(nrow(data),replace=T)
      data <- data[k, ]
    } else {}
  } else {
    
  }
  
  #make output into list
  get_data_list <- list(results_frame, dataset, data)
  return(get_data_list)
}


#####################################################################################
#### Function: get_model (2 versions)
#### There are two versions: get_model_scale and get_model_noscale
#### _scale standardises the variables before computing a regression
#### _noscale does not standardise the variables before computing the regression (for bootstrapping)
#### Note: I'm not scaling binary variables
####
#### Input: A list from the previous function containing:
#### results_frame = a data frame showing specifications on each line
#### dataset = a character saying what dataset the specifications are for
#### data = the data of interest with relevant iv and dv
####
#### Output: A list with the folliwng items
#### reg = regression model from the specification
#### results_frame = retain results_frame to use in next function
#### data = dataset of interest
####
#### Method: It goes to the row of interest (i) of the results frame.
#### There it reads what the relevant dataset and control variables,
#### this determines the model which is then run on the dataset of interest
#####################################################################################
get_model_scale <- function(get_data_list){
  
  #unpackage list
  results_frame <- get_data_list[[1]]
  dataset <- get_data_list[[2]]
  data <- get_data_list[[3]]
  
  # Model with no controls
  if (results_frame$control[i] == "none") {
    if (results_frame$model[i] == "lm") {
      reg <- lm(
        scale(dv) ~ scale(iv), 
        data = data)
    } else {
      reg <- glm(
        dv ~ scale(iv), 
        data = data, 
        family = binomial)
    }
  # Model with only early-life stress as control
  } else if (results_frame$control[i] == "ctq") {
     if (results_frame$model[i] == "lm") {
        reg <- lm(
          scale(dv) ~ scale(iv) + scale(ctq_control), 
          data = data)
     } else {
        reg <- glm(
          dv ~ scale(iv) + scale(ctq_control), 
          data = data, 
          family = binomial)
     }
  # Model with only the wave 1 mental health variable as control
  } else if (results_frame$control[i] == "mh") {
    if (results_frame$model[i] == "lm") {
      reg <- lm(
        scale(dv) ~ scale(iv) + scale(wave1_control), 
        data = data)
    } else {
      reg <- glm(
        dv ~ scale(iv) + wave1_control, 
        data = data, 
        family = binomial)
    }
  #Model with both controls  
  } else {
    if (results_frame$model[i] == "lm") {
      reg <- lm(
        scale(dv) ~ scale(iv) + scale(wave1_control) + scale(ctq_control), 
        data = data)
    } else {
      reg <- glm(
        dv ~ scale(iv) + wave1_control + scale(ctq_control),
        data = data, 
        family = binomial)
    }
  }  
  # Return model and results_frame
  get_model_list <- list(reg, results_frame, data)
  return(get_model_list)
}


get_model_noscale <- function(get_data_list){
  
  #unpackage list
  results_frame <- get_data_list[[1]]
  dataset <- get_data_list[[2]]
  data <- get_data_list[[3]]
  
  # Model with no controls
  if (results_frame$control[i] == "none") {
    if (results_frame$model[i] == "lm") {
      reg <- lm(
        dv ~ iv, 
        data = data)
    } else {
      reg <- glm(
        dv ~ iv, 
        data = data, 
        family = binomial)
    }
    # Model with only early-life stress as control
  } else if (results_frame$control[i] == "ctq") {
    if (results_frame$model[i] == "lm") {
      reg <- lm(
        dv ~ iv + ctq_control, 
        data = data)
    } else {
      reg <- glm(
        dv ~ iv + ctq_control, 
        data = data, 
        family = binomial)
    }
    # Model with only the wave 1 mental health variable as control
  } else if (results_frame$control[i] == "mh") {
    if (results_frame$model[i] == "lm") {
      reg <- lm(
        dv ~ iv + wave1_control, 
        data = data)
    } else {
      reg <- glm(
        dv ~ iv + wave1_control, 
        data = data, 
        family = binomial)
    }
    #Model with both controls  
  } else {
    if (results_frame$model[i] == "lm") {
      reg <- lm(
        dv ~ iv + wave1_control + ctq_control, 
        data = data)
    } else {
      reg <- glm(
        dv ~ iv + wave1_control + ctq_control,
        data = data, 
        family = binomial)
    }
  }  
  # Return model and results_frame
  get_model_list <- list(reg, results_frame, data)
  return(get_model_list)
}

#####################################################################################
#### Function: get_coef
#### This function extracts variables of interest from specification regression
####
#### Input: A list from the previous function with the following items
#### reg = regression model from the specification
#### results_frame = retain results_frame to use in next function
#### data = dataset of interest
####
#### Output:
#### results_frame[i,] = the relevant line of the results frame with key variables
####
#### Method: The function takes the relevant regression and results frame.
#### It extracts the t_value, regression coefficiant, p_value, SE and N
#### It adds these to the relevant results frame
#####################################################################################
get_coef <- function(get_model_list) {
  
  #unpackage list
  reg <- get_model_list[[1]]
  results_frame <- get_model_list[[2]]
  
  results_frame$t_value[i] <-
    summary(reg)$coef[[2, 3]] %>% {
      ifelse(. == 0, NA, .)
    }
  results_frame$effect[i] <-
    summary(reg)$coef[[2, 1]] %>% {
      ifelse(. == 0, NA, .)
    }
  results_frame$p_value[i] <- 
    summary(reg)$coef[[2, 4]]
  results_frame$standard_error[i] <-
    summary(reg)$coef[[2, 2]] %>% {
      ifelse(. == 0, NA, .)
    }
  results_frame$number[i] <- nobs(reg)
  results_frame$etasqrd[i] <- etasq(reg)["scale(iv)", "Partial eta^2"]
  return(results_frame[i,])
}

#####################################################################################
# a) Setup Results Frame ####
# Furthermore we add an additional results frame (suffix = boot) 
# for our bootstrapping models later in the code
#####################################################################################

###############################
#### Function: get_results_frame
#### Makes a data frame to add
#### results of SCA in
####
#### Input: 
#### outcome = measures of internalizing symptoms/disorders
#### predictor = pubertal timing measures
#### controls = categories to indicate which control variable combo to use 
####
#### Output:
#### results_frame = data frame
####
#### Method: it makes grid with every possible combo of 
#### input variables
#### changed to add a variable indicating the relevant wave 1 mh control 
#### and to add variable indicating whether to run linear or logistic regression
###############################

get_results_frame <- function(outcome, predictor, control_set, ctq){
  results_frame <- expand.grid(outcome, predictor, control_set, ctq)
  names(results_frame) <- c("outcome", "predictor", "control","ctq_control")
  results_frame$wave1_control <- str_replace(outcome,"2","1")
  results_frame$model <- ifelse(str_detect(outcome,"_d"),'logistic',"lm")
  results_frame[, c("t_value", "effect", "p_value", "standard_error", "number", "etasqrd")] <- NA
  return(results_frame)
}

###############################
# i) Define measures  ### changed to reflect our predictors, outcomes and controls
###############################

### define outcome measures for each dataset
outcome_im <- c("int_d_im_wave2","depres_d_im_wave2","anx_d_im_wave2","distress_d_im_wave2","fear_d_im_wave2",
                "CESDC_total_im_wave2","SCARED_anxiety_mean_im_wave2")
outcome_nonim <- c("int_d_wave2","depres_d_wave2","anx_d_wave2","distress_d_wave2","fear_d_wave2",
                   "CESDC_total_wave2","SCARED_anxiety_mean_wave2")

### define predictors 
predictor_im <- c("aam_final","subj_timing_im_wave1","parent_subj_timing_im_wave1","resid_neg_ldstage_im_wave1",
                  "resid_neg_pdsstage_im_wave1","resid_neg_parent_pdsstage_im_wave1","resid_neg_ADRENcomp_im_wave1","resid_neg_GONADcomp_im_wave1",
                  "resid_neg_PUBcomp_im_wave1","resid_neg_DHEA_cor_im_wave1","resid_neg_TEST_cor_im_wave1",
                  "resid_neg_EST_cor_im_wave1","subj_timing_im_wave2","parent_subj_timing_im_wave2",
                  "resid_neg_ldstage_im_wave2","resid_neg_pdsstage_im_wave2","resid_neg_parent_pdsstage_im_wave2",
                  "resid_neg_ADRENcomp_im_wave2","resid_neg_GONADcomp_im_wave2","resid_neg_PUBcomp_im_wave2",
                  "resid_neg_DHEA_cor_im_wave2","resid_neg_TEST_cor_im_wave2","resid_neg_EST_cor_im_wave2")
predictor_nonim <- c("aam_final","subj_timing_wave1","parent_subj_timing_wave1","resid_neg_ldstage_wave1",
                  "resid_neg_pdsstage_wave1","resid_neg_parent_pdsstage_wave1","resid_neg_ADRENcomp_wave1",
                  "resid_neg_GONADcomp_wave1","resid_neg_PUBcomp_wave1","resid_neg_DHEA_cor_wave1",
                  "resid_neg_TEST_cor_wave1","resid_neg_EST_cor_wave1","subj_timing_wave2",
                  "parent_subj_timing_wave2","resid_neg_ldstage_wave2","resid_neg_pdsstage_wave2",
                  "resid_neg_parent_pdsstage_wave2","resid_neg_ADRENcomp_wave2",
                  "resid_neg_GONADcomp_wave2","resid_neg_PUBcomp_wave2","resid_neg_DHEA_cor_wave2",
                  "resid_neg_TEST_cor_wave2","resid_neg_EST_cor_wave2")

### define categories to indicate which control variable combo to use (no controls, only CTQ, only wave1 mh, or both)
control_set <- c("none","ctq","mh","both")
control_im <- "CTQ_threat_im_wave1"
control_nonim <- "CTQ_threat_wave1"

###############################
# ii) Run functions ### changed so that both results_frames are combined into one
###############################
results_frame_im <- get_results_frame(outcome_im, predictor_im, control_set, control_im)
results_frame_nonim <- get_results_frame(outcome_nonim, predictor_nonim, control_set, control_nonim)
results_frame <- rbind(results_frame_im, results_frame_nonim)

###############################
# iii) Make bootstrap version
###############################
results_frame_boot <- results_frame


#####################################################################################
# b) Run Specification Curve Analyses ####
#####################################################################################
sca_full <- list(0)


if (run_sca == 1){
  
  for (b in 1:(bootstraps+1)) {
    print(paste0("bootstrap_",b))  
    for (i in 1:nrow(results_frame)) {
      results_frame[i,] <- get_coef(get_model_scale(get_data(
        results_frame = results_frame, dataset = "full", data = data_all, boot = TRUE
      )))
    }
    sca_full[[b]] <- results_frame
  }

  # Save SCAs - - NOTE: you need to rename or delete these files in the folder if run previously
  write.csv(sca_full[[bootstraps+1]], file=paste0(cas_dir,"projects/W1_W2_pubertal_timing/3_1_SCA.csv"))
  write_rds(sca_full, paste0(cas_dir,"projects/W1_W2_pubertal_timing/3_1_SCA_SCAboot.rds"))
 
} else {
  
  # Read SCAs
  results_frame <- read.csv(file=paste0(cas_dir,"projects/W1_W2_pubertal_timing/3_1_SCA.csv"))
  sca_full <- read_rds(paste0(cas_dir,"projects/W1_W2_pubertal_timing/3_1_SCA_SCAboot.rds"))
}

##########################################################################################################
# 3. Make Null Models for Bootstrapping ##################################################################
##########################################################################################################

#####################################################################################
#### Function: get_ynull
#### This function contrains model under null for each specification
####
#### Input: A list from the previous function with the following items
#### reg = regression model from the specification
#### results_frame = retain results_frame to use in next function
#### data = dataset of interest
####
#### Output:
#### y.null.i = dv data for specification constrained under null
####
#### Method: The function takes the relevant regression and results frame.
#### It extracts the regression coefficiant, and multiplies it by the iv
#### It substracts this from the dv creating a dataset where null is true
#####################################################################################
get_ynull <- function(get_model_list) {
  
  #unpackage list
  reg <- get_model_list[[1]]
  data <- get_model_list[[3]]
  
  if (results_frame$model[i] == "lm") {
    # extract coefficient
    b.i <-
     summary(reg)$coef[[2, 1]] %>% {
        ifelse(. == 0, NA, .)
     }
    #make null model
    y.null.i <- data$dv-(b.i*data$iv)

#We use an alternative approach for binary y variables, predicting each person's probability of having a 1 when the effect of iv is zero, 
#then regenerating observed scores (0, 1) from the probability    
  } else {
    if (results_frame$control[i]=="none") {
      y.null.i = arm::invlogit(coef(reg)[[1]] + 0*data$wave1_control)
    }
    else if (results_frame$control[i]=="mh") {
      control_coef <- replace(coef(reg)[[3]], is.na(coef(reg)[[3]]), 0) 
      y.null.i = arm::invlogit(coef(reg)[[1]] + control_coef*data$wave1_control)
    }
    else if (results_frame$control[i]=="ctq") {
      control_coef <- replace(coef(reg)[[3]], is.na(coef(reg)[[3]]), 0) 
      y.null.i = arm::invlogit(coef(reg)[[1]] + control_coef*data$ctq_control)
    }
    else {
      control_coef <- replace(coef(reg)[[3]], is.na(coef(reg)[[3]]), 0) 
      y.null.i = arm::invlogit(coef(reg)[[1]] + control_coef*data$wave1_control + coef(reg)[[4]]*data$ctq_control)
    }
  }
  
  return(y.null.i)
}

###############################
# i) Run Null Models
###############################
y_null <- list(0)

if (run_null == 1){
  
  for (i in 1:nrow(results_frame)) {
    print(paste0("NULL_", i))
    y_null[[i]] <- get_ynull(get_model_noscale(get_data(
      results_frame = results_frame, dataset = "full", data = data_all, boot = FALSE
    )))
  }
  
  # Save null models - NOTE: you need to rename or delete these files in the folder if run previously
  write_rds(y_null, paste0(cas_dir,"projects/W1_W2_pubertal_timing/3_1_ynull.rds"))

} else {
  
  # Read null models
  y_null <- read_rds(paste0(cas_dir,"projects/W1_W2_pubertal_timing/3_1_ynull.rds"))

}


##########################################################################################################
# 4. Run Bootstrap Models ################################################################################
##########################################################################################################

#####################################################################################
#### Function: get_boot_data
#### This function takes the specification and makes the apporiate bootstrapped dataset
#### Except the last round trough (b >= bootstraps) it runs the original function
####
#### Input: 
#### results_frame = a data frame showing specifications on each line
#### dataset = a character saying what dataset the specifications are for
#### data = the data of interest
#### y.null = list of null datasets for each specification
####
#### Output: A list with the folliwng items
#### results_frame = retain results_frame to use in next function
#### dataset = retain dataset to use in next function
#### data = dataset of interest, now with new dv and iv for specification constrained under null
####
#### Method: It goes to the row of interest (i) of the results frame.
#### There it reads what the relevant iv and dv variables are and extracts these
#### from the dataset of interest.
#### During the bootstrapping tests (b <= bootstraps) it makes a bootstrapped dataset
#### that is constrained under the null. In the last test, it makes non-bootstrapped dataset
####
#### changed to find the iv, dv and control variables in out data
#####################################################################################
get_boot_data <- function(results_frame, dataset, data, y.null) {
  
  # Setup variable names iv and controls
  data$iv <- data[,as.character(results_frame$predictor[i])]
  data$wave1_control <- data[,as.character(results_frame$wave1_control[i])]
  data$ctq_control <- data[,as.character(results_frame$ctq_control[i])]
 
  # Setup variable names dv
  if (b <= bootstraps){
    data$dv <- y.null[[i]]
  } else {
   data$dv <- data[,as.character(results_frame$outcome[i])]
  }

  if (results_frame$model[[i]]=="lm"){
    if (b <= bootstraps){
      k=sample(nrow(data),replace=T)
      data <- data[k, ]
    } else {}
  } else {
    if (b <= bootstraps){
      data$dv = rbinom(size = 1, prob = y.null[[i]], n = samplesize)
      k=sample(nrow(data),replace=T)
      data <- data[k, ]
    } else {}
  }
 
  
  # Do bootstrapping
  
  
  #make output into list
  get_data_list <- list(results_frame, dataset, data)
  return(get_data_list)
}

#####################################################################################
# a) Run Bootstrapping ####
#####################################################################################
bootstraps_full <- list(0)

if (run_boot == 1){
  
  for (b in 1:(bootstraps+1)) {
    print(paste0("bootstrap_",b))
    results_frame_it <- results_frame_boot
    for (i in 1:nrow(results_frame_it)) {
      results_frame_it[i,] <- get_coef(get_model_noscale(get_boot_data(
        results_frame = results_frame_it, dataset = "full", data = data_all, y.null = y_null
      )))
    }
    bootstraps_full[[b]] <- results_frame_it
  }
  # save - NOTE: you need to rename or delete these files in the folder if run previously
  write_rds(bootstraps_full, paste0(cas_dir,"projects/W1_W2_pubertal_timing/3_1_boot.rds"))

} else {
  
  # Read null models
  bootstraps_full <- read_rds(paste0(cas_dir,"projects/W1_W2_pubertal_timing/3_1_boot.rds"))

}
