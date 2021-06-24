##########################################################################################################
### Script 3b: Make Results Table ####
### This scripts takes the three completed SCAs and Bootstraps 
### It creates a table of results for all their results and also aggregate effects.
##########################################################################################################

##########################################################################################################
# Loading libraries, controlling analyses and loading data ############################################
##########################################################################################################

#####################################################################################
# a) Load libraries ####
#####################################################################################
library("plyr")
library("tidyverse")
cas_dir="Y:/dsnlab/TAG/"

## set to how many bootstraps I want to do
bootstraps <- 1000

#####################################################################################
# b) Load SCAs, SCAs with bootstrap and Bootstrapped Models ####
#####################################################################################  
results_frame <- read.csv(file=paste0(cas_dir,"projects/W1_W2_pubertal_timing/3_1_SCA.csv"))

boot_sca <- read_rds(paste0(cas_dir,"projects/W1_W2_pubertal_timing/3_1_SCA_SCAboot.rds"))

bootstraps_full <- read_rds(paste0(cas_dir,"projects/W1_W2_pubertal_timing/3_1_boot.rds"))

#####################################################################################
# b) At the beginning of this file it makes sense to add columns                 ####
#    to the SCA frames that show the sign and significance of each specification ####
#####################################################################################  
make_sign <- function(results_frame){
  results_frame$sign <- ifelse(results_frame$effect > 0, "positive", "negative")
  results_frame$sig_sign <- ifelse(results_frame$p_value > 0.05, NA, ifelse(results_frame$effect > 0, "sig_positive", "sig_negative"))
  results_frame$sign <- as.factor(results_frame$sign)
  results_frame$sig_sign <- as.factor(results_frame$sig_sign)
  return(results_frame)
}

results_frame <- make_sign(results_frame)

##########################################################################################################
# 1. Make Bootstrapped Model Table: ######################################################################
##########################################################################################################

#####################################################################################
#### Function: fill_per
#### This function takes the results of the bootstrapped models (a list of results
#### frames) and summarises them into one permutation frame
####
#### Input: 
#### boot_results = a list of results_frames from the bootstrapping procedure
#### frame = permutation frame which can hold summarised results
####
#### Output:
#### frame = permutation frame with summarised bootstrapped results
####
#### Method: The function takes the bootstrapped results and goes through them 
#### individually (m). For each tech variable and control (defined by filter value, 
#### i) it filters the bootstrapped result frame so that it only examines 
#### specifications with the relevant outcome or predictor values.It first gets the median 
#### size of the effect and then the amount of effects that are positivie or negatve. 
#### Lastly it examines the amount of results that are significant and positive or negative. 
#####################################################################################
fill_per_o <- function(boot_results, frame){
  for (m in 1:(bootstraps+1)){
    
    boot_results_m <- as.data.frame(boot_results[m])
    n <- 2
    
    for (i in 1:8){
      
      if(i == 1){
        results_subset <- boot_results_m
      } else if (i == 2) {
        results_subset <- boot_results_m %>% filter(grepl(pattern="int_d", outcome))
      } else if (i == 3) {
        results_subset <- boot_results_m %>% filter(grepl(pattern="depres_d", outcome))
      } else if (i == 4){
        results_subset <- boot_results_m %>% filter(grepl(pattern="anx_d", outcome))
      } else if (i == 5){
        results_subset <- boot_results_m %>% filter(grepl(pattern="distress_d", outcome))
      } else if (i == 6){
        results_subset <- boot_results_m %>% filter(grepl(pattern="fear_d", outcome))
      } else if (i == 7){
        results_subset <- boot_results_m %>% filter(grepl(pattern="CESDC", outcome))
      } else if (i == 8){
        results_subset <- boot_results_m %>% filter(grepl(pattern="SCARED", outcome))
      }
      
      frame[m, n] <- median(results_subset[["effect"]], na.rm = TRUE)
      n <- n+1
      frame[m, n] <- length(results_subset[results_subset$effect < 0, "effect"])
      n <- n+1
      frame[m, n] <- length(results_subset[results_subset$effect > 0, "effect"])
      n <- n+1
      sig_data <- filter(results_subset, p_value < 0.05)
      
      if(nrow(sig_data) > 0){
        frame[m, n] <- length(sig_data[sig_data$effect < 0, "effect"])
        n <- n+1
        frame[m, n] <- length(sig_data[sig_data$effect > 0, "effect"])
        n <- n+1
      } else {
        frame[m, n] <- 0
        n <- n+1
        frame[m, n] <- 0
        n <- n+1
      }
    }
  }
  
  return(frame)
}

#####################################################################################
# a) Make dataframe to hold the results ####
#####################################################################################
permutation_frame_o <-
  data.frame(matrix(NA, nrow = bootstraps, ncol = 41))
names(permutation_frame_o) <-
  c("permutation_number",
    "effect.tot", "sign.neg.tot", "sign.pos.tot", "sign.sig.neg.tot", "sign.sig.pos.tot",
    
    #Internalizing disorder
    "effect.int","sign.neg.int","sign.pos.int","sign.sig.neg.int","sign.sig.pos.int",
    
    #Depressive disorder
    "effect.depr","sign.neg.depr","sign.pos.depr","sign.sig.neg.depr","sign.sig.pos.depr",
    
    #Anxiety disorder
    "effect.anx","sign.neg.anx","sign.pos.anx","sign.sig.neg.anx","sign.sig.pos.anx",
    
    #Distress disorder HiTOP
    "effect.distr","sign.neg.distr","sign.pos.distr","sign.sig.neg.distr","sign.sig.pos.distr",
    
    #Fear disorder HiTOP
    "effect.fear","sign.neg.fear","sign.pos.fear","sign.sig.neg.fear","sign.sig.pos.fear",
    
    #CESDC depression symptoms
    "effect.cesdc","sign.neg.cesdc","sign.pos.cesdc","sign.sig.neg.cesdc","sign.sig.pos.cesdc",
    
    #SCARED anxziety symptoms
    "effect.scar","sign.neg.scar","sign.pos.scar","sign.sig.neg.scar","sign.sig.pos.scar"
 
    )

permutation_frame_o$permutation_number <- seq.int(bootstraps)
per_frame_o <- permutation_frame_o

#####################################################################################
# b) run fill_per function to fill permutation frames ####
#####################################################################################
per_frame_o <- fill_per_o(bootstraps_full, per_frame_o)

#####################################################################################
# c) anaylse the permutation frame to test significance ####
#####################################################################################

#####################################################################################
#### Function: analyse_boot
#### The function analyses the bootstrapping models to test significance
####
#### Input: 
#### boot_results = the permutation frame summarising the results of the bootstraps
#### 
#### Output:
#### table = a data table summarising the results of the bootstraps for
####
#### Method: For each column of the specified tests we examine how many bootstrap
#### tests have larger effects than the original sca. Then we examine how many have
#### more effects in the direction that most effects go and lastly we examine 
#### whether they have more significant effects in the major direction.
#####################################################################################
analyse_boot_o <- function(boot_results){

  results_table <- data.frame(matrix(NA, nrow = 8, ncol = 3))
  names(results_table) <- c("median_effect","share_neg","share_sig_neg")
  results_table$outc_var <- c("tot","int","depr","anx","distr","fear","cesdc","scar")

  for(i in 1:nrow(results_table)){
    suffix <- results_table$outc_var[i]
    
    ### 1. effect sizes
    results_table[i,1] <-
           mean(abs(boot_results[1:bootstraps, paste0("effect.", suffix)]) >= abs(boot_results[(bootstraps+1), paste0("effect.", suffix)]))
    
    ### 2. sign of effect 
    sign <- pmax(boot_results[1:bootstraps, paste0("sign.neg.", suffix)], boot_results[1:bootstraps, paste0("sign.pos.", suffix)])
    sign_obs <- pmax(boot_results[(bootstraps+1), paste0("sign.neg.", suffix)], boot_results[(bootstraps+1), paste0("sign.pos.", suffix)])
    results_table[i,2] <-
           mean(sign >= sign_obs)
    
    ### 3. sign of significant effects
    sign <- pmax(boot_results[1:bootstraps, paste0("sign.sig.neg.", suffix)], boot_results[1:bootstraps, paste0("sign.sig.pos.", suffix)])
    sign_obs <- pmax(boot_results[(bootstraps+1), paste0("sign.sig.neg.", suffix)], boot_results[(bootstraps+1), paste0("sign.sig.pos.", suffix)])
    results_table[i,3] <-
           mean(sign >= sign_obs)
    
  }

  return(results_table)
}

bootstrap_summary_o <- analyse_boot_o(per_frame_o)
#bootstrap_summary[,1:ncol(bootstrap_summary)] <- sapply(bootstrap_summary[,1:ncol(bootstrap_summary)],as.numeric)


##########################################################################################################
# 2. Make Table  2 #######################################################################################
##########################################################################################################
results_table_o <- data.frame(outcome = NA, sig_measure = NA,
                            observed = NA, lower = NA, upper = NA, p = NA)
empty_rows <- data.frame(matrix(nrow = 20, ncol = ncol(results_table_o)))
names(empty_rows) <- names(results_table_o)
results_table_o <- rbind(results_table_o, empty_rows)

results_table_o$outcome <- c("int_d","depres_d","anx_d",
                           "distress_d", "fear_d", "CESDC_total", "SCARED_anxiety_mean")
results_table_o$sig_measure <- rep(c("median ES", "share of results", "share of sig results"),each= 7)

#make a results frame with variable names collapsed over imp and wave
results_frame_forwide <- results_frame %>% 
  mutate(wave=ifelse(grepl("wave2", results_frame$predictor),"wave2","wave1"),
         imputation=ifelse(grepl("_im_", results_frame$outcome),"imp","nonimp"))

results_frame_forwide$predictor <- sub("_im_wave[1:2]$", "", results_frame_forwide$predictor)
results_frame_forwide$predictor <- sub("_wave[1:2]$", "", results_frame_forwide$predictor)
results_frame_forwide$outcome <- sub("_im_wave2$", "", results_frame_forwide$outcome)
results_frame_forwide$outcome <- sub("_wave2$", "", results_frame_forwide$outcome)

##########################################################################################################
# 3. Populate Table ######################################################################################
##########################################################################################################

### this function populates the results table with the observed variables 
### (observed median point estimate, share of results or share of sig results)
populate_table_observed_o <- function(results_frame, table){
  
  medians <- results_frame %>% group_by(outcome) %>% dplyr::summarise(median = median(effect))
  signs <- results_frame %>% group_by(outcome, sign) %>% dplyr::summarise(counts = n())
  sig_signs <- results_frame %>% group_by(outcome, sig_sign) %>% dplyr::summarise(counts = n()) %>% filter(is.na(sig_sign) == FALSE)
    
    for(i in 1:nrow(table)){
      if((table[i,"sig_measure"]) == "median ES"){
        table[i, "observed"] <- medians %>% filter(outcome == table[i, "outcome"]) %>% pull(median) %>% round(2)
      } else if ((table[i,"sig_measure"]) == "share of results") {
        table[i, "observed"] <- signs %>% filter(outcome == table[i, "outcome"]) %>% dplyr::summarise(max = max(counts)) %>% pull(max) %>% round(1)
      } else if ((table[i,"sig_measure"]) == "share of sig results") {
        sig_results <- results_frame %>% filter(outcome == table[i, "outcome"]) %>% filter(is.na(sig_sign) == FALSE)
        sig_results_sign <- sig_signs %>% filter(outcome == table[i, "outcome"]) %>% dplyr::summarise(max = max(counts)) %>% pull(max) %>% round(1)
        table[i, "observed"] <- ifelse(nrow(sig_results) == 0, 0, sig_results_sign)
      } else {
        table[i, "observed"] <- NA
      }
    }
  
  return(table)
}


### this function populates the results table with the bootstrapped p values for the different significance tests
populate_table_bootstrap_o <- function(bootstrap, table){
  
  for (i in 1:nrow(table)){
    if (table[i, "outcome"] == "anx_d") {
        outc <- "anx"
    } else if (table[i, "outcome"] == "CESDC_total") {
      outc <- "cesdc"
    } else if (table[i, "outcome"] == "depres_d") {
      outc <- "depr"
    } else if (table[i, "outcome"] == "distress_d") {
      outc <- "distr"
    } else if (table[i, "outcome"] == "fear_d") {
      outc <- "fear"
    } else if (table[i, "outcome"] == "int_d") {
      outc <- "int"
    } else if (table[i, "outcome"] == "SCARED_anxiety_mean") {
      outc <- "scar"
    }   
     
    if((table[i,"sig_measure"]) == "median ES"){
        table[i, "p"] <- bootstrap[which(bootstrap$outc_var== outc),"median_effect"] %>% round(3)
    } else if ((table[i,"sig_measure"]) == "share of results") {
        table[i,"p"] <- bootstrap[which(bootstrap$outc_var== outc),"share_neg"] %>% round(3)
    } else if ((table[i,"sig_measure"]) == "share of sig results") {
        table[i, "p"] <- bootstrap[which(bootstrap$outc_var== outc),"share_sig_neg"] %>% round(3)
    } 

  }
  return(table)
}

### this function populates the results table with the bootstrapped confidence intervals for the median point estimates
populate_table_bounds_o <- function(boot, table){

  results_frame_medians <- data.frame("anx_d" = NA, "CESDC_total" = NA,"depres_d" = NA,
                                      "distress_d" = NA, "fear_d" = NA, "int_d" = NA, "SCARED_anxiety_mean" = NA) #needs to be alphabetical order 

  for(b in 1:(length(boot)-1)){
    boot_sh <- boot[[b]] %>% mutate(outcome= sub("_wave[1:2]|_im_wave[1:2]$", "", outcome))  
    medians <- boot_sh %>% group_by(outcome) %>% dplyr::summarise(median = median(effect))
    results_frame_medians <- rbind(results_frame_medians, as.numeric(t(medians)[2,]))
  }
    
  results_frame_medians <- results_frame_medians[2:nrow(results_frame_medians),]
    
  for(i in 1:nrow(table)){
      if((table[i,"sig_measure"]) == "median ES"){
        table[i, c("lower", "upper")] <- round(quantile(results_frame_medians[, paste0(table[i, "outcome"])], probs = c(0.025, 0.975),na.rm=T), 2)
      } else {
        table[i, c("lower", "upper")] <- NA
      }
  }

  return(table)
}


results_table_complete_o <- populate_table_observed_o(results_frame = results_frame_forwide, 
                                                  table = results_table_o)
results_table_complete_o <- populate_table_bootstrap_o(bootstrap_summary_o, results_table_complete_o)
results_table_complete_o <- populate_table_bounds_o(boot = boot_sca, results_table_complete_o)


##########################################################################################################
# 4. Save Table ##########################################################################################
##########################################################################################################
write.csv(results_table_complete_o, file = paste0(cas_dir,"projects/W1_W2_pubertal_timing/table_byoutcome.csv"))


