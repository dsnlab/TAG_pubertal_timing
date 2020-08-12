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
bootstraps <- 500

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
fill_per <- function(boot_results, frame){
  for (m in 1:(bootstraps+1)){
    
    boot_results_m <- as.data.frame(boot_results[m])
    n <- 2
    
    for (i in 1:13){
      
      if(i == 1){
        results_subset <- boot_results_m
      } else if (i == 2) {
        results_subset <- boot_results_m %>% filter(grepl(pattern="PUBcomp", predictor))
      } else if (i == 3) {
        results_subset <- boot_results_m %>% filter(grepl(pattern="ADRENcomp", predictor))
      } else if (i == 4){
        results_subset <- boot_results_m %>% filter(grepl(pattern="GONADcomp", predictor))
      } else if (i == 5){
        results_subset <- boot_results_m %>% filter(grepl(pattern="neg_pdsstage", predictor))
      } else if (i == 6){
        results_subset <- boot_results_m %>% filter(grepl(pattern="parent_pdsstage", predictor))
      } else if (i == 7){
        results_subset <- boot_results_m %>% filter(grepl(pattern="ldstage", predictor))
      } else if (i == 8){
        results_subset <- boot_results_m %>% filter(grepl(pattern="aam", predictor))
      } else if (i == 9){
        results_subset <- boot_results_m %>% filter(grepl(pattern="^subj_timing", predictor))
      } else if (i == 10){
        results_subset <- boot_results_m %>% filter(grepl(pattern="parent_subj_timing", predictor))
      } else if (i == 11){
        results_subset <- boot_results_m %>% filter(grepl(pattern="DHEA", predictor))
      } else if (i == 12){
        results_subset <- boot_results_m %>% filter(grepl(pattern="TEST", predictor))
      } else if (i == 13){
        results_subset <- boot_results_m %>% filter(grepl(pattern="_EST", predictor))
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
permutation_frame <-
  data.frame(matrix(NA, nrow = bootstraps, ncol = 66))
names(permutation_frame) <-
  c("permutation_number",
    "effect.tot", "sign.neg.tot", "sign.pos.tot", "sign.sig.neg.tot", "sign.sig.pos.tot",
    
    #PUBcomp
    "effect.pcomp","sign.neg.pcomp","sign.pos.pcomp","sign.sig.neg.pcomp","sign.sig.pos.pcomp",
    
    #ADRENcomp
    "effect.acomp","sign.neg.acomp","sign.pos.acomp","sign.sig.neg.acomp","sign.sig.pos.acomp",
    
    #GONADcomp
    "effect.gcomp","sign.neg.gcomp","sign.pos.gcomp","sign.sig.neg.gcomp","sign.sig.pos.gcomp",
    
    #pdsstage
    "effect.pds","sign.neg.pds","sign.pos.pds","sign.sig.neg.pds","sign.sig.pos.pds",
    
    #parent pdsstage
    "effect.ppds","sign.neg.ppds","sign.pos.ppds","sign.sig.neg.ppds","sign.sig.pos.ppds",
    
    #ldstage
    "effect.ld","sign.neg.ld","sign.pos.ld","sign.sig.neg.ld","sign.sig.pos.ld",
    
    #age at menarche
    "effect.aam","sign.neg.aam","sign.pos.aam","sign.sig.neg.aam","sign.sig.pos.aam",
    
    #subjective timing
    "effect.subj","sign.neg.subj","sign.pos.subj","sign.sig.neg.subj","sign.sig.pos.subj",
    
    #parent subjective timing
    "effect.psubj","sign.neg.psubj","sign.pos.psubj","sign.sig.neg.psubj","sign.sig.pos.psubj",
    
    #DHEA
    "effect.dhea","sign.neg.dhea","sign.pos.dhea","sign.sig.neg.dhea","sign.sig.pos.dhea",
    
    #TEST
    "effect.test","sign.neg.test","sign.pos.test","sign.sig.neg.test","sign.sig.pos.test",
    
    #EST
    "effect.est","sign.neg.est","sign.pos.est","sign.sig.neg.est","sign.sig.pos.est"
    )

permutation_frame$permutation_number <- seq.int(bootstraps)
per_frame <- permutation_frame

#####################################################################################
# b) run fill_per function to fill permutation frames ####
#####################################################################################
per_frame <- fill_per(bootstraps_full, per_frame)

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
analyse_boot <- function(boot_results){

  results_table <- data.frame(matrix(NA, nrow = 13, ncol = 3))
  names(results_table) <- c("median_effect","share_neg","share_sig_neg")
  results_table$timing_var <- c("tot","pcomp","acomp","gcomp","pds","ppds","ld","aam","subj","psubj","dhea","test","est")

  for(i in 1:nrow(results_table)){
    suffix <- results_table$timing_var[i]
    
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

bootstrap_summary <- analyse_boot(per_frame)
#bootstrap_summary[,1:ncol(bootstrap_summary)] <- sapply(bootstrap_summary[,1:ncol(bootstrap_summary)],as.numeric)


##########################################################################################################
# 2. Make Table  2 #######################################################################################
##########################################################################################################
results_table2 <- data.frame(predictor = NA, sig_measure = NA,
                            observed = NA, lower = NA, upper = NA, p = NA)
empty_rows <- data.frame(matrix(nrow = 35, ncol = ncol(results_table2)))
names(empty_rows) <- names(results_table2)
results_table2 <- rbind(results_table2, empty_rows)

results_table2$predictor <- c("aam_final","parent_subj_timing","resid_neg_ADRENcomp",
                           "resid_neg_DHEA_cor", "resid_neg_EST_cor", "resid_neg_GONADcomp", 
                              "resid_neg_ldstage", "resid_neg_parent_pdsstage" , "resid_neg_pdsstage",
                              "resid_neg_PUBcomp", "resid_neg_TEST_cor", "subj_timing")
results_table2$sig_measure <- rep(c("median ES", "share of results", "share of sig results"),each= 12)

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
populate_table_observed <- function(results_frame, table){
  
  medians <- results_frame %>% group_by(predictor) %>% dplyr::summarise(median = median(effect))
  signs <- results_frame %>% group_by(predictor, sign) %>% dplyr::summarise(counts = n())
  sig_signs <- results_frame %>% group_by(predictor, sig_sign) %>% dplyr::summarise(counts = n()) %>% filter(is.na(sig_sign) == FALSE)
    
    for(i in 1:nrow(table)){
      if((table[i,"sig_measure"]) == "median ES"){
        table[i, "observed"] <- medians %>% filter(predictor == table[i, "predictor"]) %>% pull(median) %>% round(2)
      } else if ((table[i,"sig_measure"]) == "share of results") {
        table[i, "observed"] <- signs %>% filter(predictor == table[i, "predictor"]) %>% dplyr::summarise(max = max(counts)) %>% pull(max) %>% round(1)
      } else if ((table[i,"sig_measure"]) == "share of sig results") {
        sig_results <- results_frame %>% filter(predictor == table[i, "predictor"]) %>% filter(is.na(sig_sign) == FALSE)
        sig_results_sign <- sig_signs %>% filter(predictor == table[i, "predictor"]) %>% dplyr::summarise(max = max(counts)) %>% pull(max) %>% round(1)
        table[i, "observed"] <- ifelse(nrow(sig_results) == 0, 0, sig_results_sign)
      } else {
        table[i, "observed"] <- NA
      }
    }
  
  return(table)
}


### this function populates the results table with the bootstrapped p values for the different significance tests
populate_table_bootstrap <- function(bootstrap, table){
  
  for (i in 1:nrow(table)){
    if (table[i, "predictor"] == "aam_final") {
        pubtim <- "aam"
    } else if (table[i, "predictor"] == "subj_timing") {
      pubtim <- "subj"
    } else if (table[i, "predictor"] == "parent_subj_timing") {
      pubtim <- "psubj"
    } else if (table[i, "predictor"] == "resid_neg_ldstage") {
      pubtim <- "ld"
    } else if (table[i, "predictor"] == "resid_neg_pdsstage") {
      pubtim <- "pds"
    } else if (table[i, "predictor"] == "resid_neg_parent_pdsstage") {
      pubtim <- "ppds"
    } else if (table[i, "predictor"] == "resid_neg_ADRENcomp") {
      pubtim <- "acomp"
    } else if (table[i, "predictor"] == "resid_neg_GONADcomp") {
      pubtim <- "gcomp"
    } else if (table[i, "predictor"] == "resid_neg_PUBcomp") {
      pubtim <- "pcomp"
    } else if (table[i, "predictor"] == "resid_neg_DHEA_cor") {
      pubtim <- "dhea"
    } else if (table[i, "predictor"] == "resid_neg_TEST_cor") {
      pubtim <- "test"
    } else if (table[i, "predictor"] == "resid_neg_EST_cor") {
      pubtim <- "est"
    }   
     
    if((table[i,"sig_measure"]) == "median ES"){
        table[i, "p"] <- bootstrap[which(bootstrap$timing_var== pubtim),"median_effect"] %>% round(3)
    } else if ((table[i,"sig_measure"]) == "share of results") {
        table[i,"p"] <- bootstrap[which(bootstrap$timing_var== pubtim),"share_neg"] %>% round(3)
    } else if ((table[i,"sig_measure"]) == "share of sig results") {
        table[i, "p"] <- bootstrap[which(bootstrap$timing_var== pubtim),"share_sig_neg"] %>% round(3)
    } 

  }
  return(table)
}

### this function populates the results table with the bootstrapped confidence intervals for the median point estimates
populate_table_bounds <- function(boot, table){

  results_frame_medians <- data.frame("aam_final" = NA, "parent_subj_timing" = NA,"resid_neg_ADRENcomp" = NA,
                                      "resid_neg_DHEA_cor" = NA, "resid_neg_EST_cor" = NA, "resid_neg_GONADcomp" = NA, 
                                        "resid_neg_ldstage" = NA, "resid_neg_parent_pdsstage" = NA, "resid_neg_pdsstage" = NA,
                                        "resid_neg_PUBcomp" = NA, "resid_neg_TEST_cor" = NA, "subj_timing" = NA) #needs to be alphabetical order 
    
  for(b in 1:(length(boot)-1)){
    boot_sh <- boot[[b]] %>% mutate(predictor= sub("_wave[1:2]|_im_wave[1:2]$", "", predictor))  
    medians <- boot_sh %>% group_by(predictor) %>% dplyr::summarise(median = median(effect))
    results_frame_medians <- rbind(results_frame_medians, as.numeric(t(medians)[2,]))
  }
    
  results_frame_medians <- results_frame_medians[2:nrow(results_frame_medians),]
    
  for(i in 1:nrow(table)){
      if((table[i,"sig_measure"]) == "median ES"){
        table[i, c("lower", "upper")] <- round(quantile(results_frame_medians[, paste0(table[i, "predictor"])], probs = c(0.025, 0.975),na.rm=T), 2)
      } else {
        table[i, c("lower", "upper")] <- NA
      }
  }

  return(table)
}


results_table_complete <- populate_table_observed(results_frame = results_frame_forwide, 
                                                  table = results_table2)
results_table_complete <- populate_table_bootstrap(bootstrap_summary, results_table_complete)
results_table_complete <- populate_table_bounds(boot = boot_sca, results_table_complete)


##########################################################################################################
# 4. Save Table ##########################################################################################
##########################################################################################################

write.csv(results_table_complete, file = paste0(cas_dir,"projects/W1_W2_pubertal_timing/table_complete.csv"))


##########################################################################################################
# 5. Effect of model decisions
# 
# 
##########################################################################################################


################# imputed or not #################################

results_frame_wide_impu <- results_frame_forwide %>% 
  select (predictor, outcome, control, wave, imputation, effect) %>% 
  spread(imputation, effect)

median(results_frame_wide_impu$imp)
median(results_frame_wide_impu$nonimp)

#paired t test
t.test(results_frame_wide_impu$imp, 
       results_frame_wide_impu$nonimp, 
       paired=TRUE, 
       conf.level=0.95)
#non-parametric ttest
wilcox.test(results_frame_wide_impu$imp, 
       results_frame_wide_impu$nonimp, 
       paired=TRUE, 
       conf.level=0.95)
#standard errors
results_frame_wide_impu2 <- results_frame_forwide %>% 
  select (predictor, outcome, control, wave, imputation, standard_error) %>% 
  spread(imputation, standard_error)
median(results_frame_wide_impu2$imp)
median(results_frame_wide_impu2$nonimp)

#pulling the median difference between models with and without imputation from bootstrapped null models
results_frame_imps <- data.frame("diff_imp_full" = NA)
for(b in 1:(length(bootstraps_full)-1)){
  boot_sh <- bootstraps_full[[b]] 
  boot_sh1 <- boot_sh %>% filter(grepl(pattern="_im", predictor)) %>%
    mutate(predictor= sub("_im", "", predictor),outcome= sub("_im", "", outcome),effect_imp=effect) %>%
    dplyr::select(outcome,predictor,control,effect_imp)
  boot_sh2 <- boot_sh %>% filter(!grepl(pattern="_im", predictor)) %>%
    mutate(effect_ful=effect) %>%
    dplyr::select(outcome,predictor,control,effect_ful)
  boot_2sh <- merge(boot_sh1,boot_sh2,by=c("outcome","predictor","control")) %>%
    mutate(effect_diff=effect_imp - effect_ful)
  median <- median(boot_2sh$effect_diff)
  results_frame_imps <- rbind(results_frame_imps, as.numeric(median))
}
results_frame_imps <- results_frame_imps[2:nrow(results_frame_imps),]

round(quantile(results_frame_imps, probs = c(0.025, 0.975),na.rm=T), 3)
obs_diff_impu <- median(results_frame_wide_impu$imp) - median(results_frame_wide_impu$nonimp)
mean(results_frame_imps <= obs_diff_impu)
#observed difference -0.195, p=.07, CI from bootstrapped null models -0.026 to 0.029


################# testing the effect of including control variables #################################

### find medians
effects_by_control <- results_frame %>% 
  dplyr::group_by(as.factor(control)) %>%
  dplyr::summarize(median_effect = median(effect))

results_frame_wide_con <- results_frame_forwide %>% 
  select (predictor, outcome, control, wave, imputation, effect) %>% 
  spread(control, effect)
results_frame_wide_con$speci <- 1:nrow(results_frame_wide_con)
results_frame_long_con <-results_frame_wide_con %>% 
  pivot_longer(cols=c("both","ctq","mh","none"),names_to = "control",values_to = "effect")

#nonparametric repeated measures anova
friedman.test(effect ~ control | speci, data=results_frame_long_con)

#pairwise comparisons
results_frame_long_con %>% filter(control=="none"|control=="mh") %>% wilcox.test(effect ~ control,data=.,paired=TRUE,conf.level=0.95)
results_frame_long_con %>% filter(control=="none"|control=="ctq") %>% wilcox.test(effect ~ control,data=.,paired=TRUE,conf.level=0.95)
results_frame_long_con %>% filter(control=="none"|control=="both") %>% wilcox.test(effect ~ control,data=.,paired=TRUE,conf.level=0.95)
results_frame_long_con %>% filter(control=="ctq"|control=="both") %>% wilcox.test(effect ~ control,data=.,paired=TRUE,conf.level=0.95)
results_frame_long_con %>% filter(control=="mh"|control=="both") %>% wilcox.test(effect ~ control,data=.,paired=TRUE,conf.level=0.95)
results_frame_long_con %>% filter(control=="ctq"|control=="mh") %>% wilcox.test(effect ~ control,data=.,paired=TRUE,conf.level=0.95)

#pulling the median difference between models with and without the ctq/mh control from bootstrapped null models
results_frame_co <- data.frame("mh" = NA,"ctq" = NA)
for(b in 1:(length(bootstraps_full)-1)){
  boot_sh <- bootstraps_full[[b]] 
  boot_sh0 <- boot_sh %>% filter(control=="none") %>%
    mutate(effect_0=effect) %>%
    dplyr::select(outcome,predictor,effect_0)
  for(c in 1:ncol(results_frame_co)){
      boot_sh1 <- boot_sh %>% filter(control==colnames(results_frame_co[c])) %>%
        mutate(effect_1=effect) %>%
       dplyr::select(outcome,predictor,effect_1)
      boot_2sh <- merge(boot_sh0,boot_sh1,by=c("outcome","predictor")) %>%
        mutate(effect_diff=effect_1 - effect_0)
      median[c] <- median(boot_2sh$effect_diff)
  }    
  results_frame_co <- rbind(results_frame_co, as.numeric(median))
}
results_frame_co <- results_frame_co[2:nrow(results_frame_co),]
#mh versus none #observed difference=0.011, p=.30, and CI from bootstrapped null models = -0.037 to 0.041
round(quantile(results_frame_co["mh"], probs = c(0.025, 0.975),na.rm=T), 3)
obs_diff_mh <- effects_by_control$median_effect[3]-effects_by_control$median_effect[4]
mean(results_frame_co >= obs_diff_mh)
#ctq versus none #observed difference=0.022, p=.14, and CI from bootstrapped null models = -0.038 to 0.04
round(quantile(results_frame_co["ctq"], probs = c(0.025, 0.975),na.rm=T), 3)
obs_diff_ctq <- effects_by_control$median_effect[2]-effects_by_control$median_effect[4]
mean(results_frame_co >= obs_diff_ctq)

############################# wave 1 or wave 2 pubertal timing #################################

results_frame_wide_wave <- results_frame_forwide %>% 
  select (predictor, outcome, control, wave, imputation, effect) %>% 
  filter(!predictor=="aam_final") %>%
  spread(wave, effect)

median(results_frame_wide_wave$wave1)
median(results_frame_wide_wave$wave2)

#paired t test
t.test(results_frame_wide_wave$wave1, 
       results_frame_wide_wave$wave2, 
       paired=TRUE, 
       conf.level=0.95)
#non-parametric ttest
wilcox.test(results_frame_wide_wave$wave1, 
        results_frame_wide_wave$wave2, 
        paired=TRUE, 
        conf.level=0.95)

#pulling the median difference between prospective and cross-sectional models from bootstrapped null models
results_frame_diffs <- data.frame("diff_pros_cross" = NA)
for(b in 1:(length(bootstraps_full)-1)){
  boot_sh <- bootstraps_full[[b]] 
  boot_sh1 <- boot_sh %>% filter(grepl(pattern="wave1", predictor)) %>%
                          mutate(predictor= sub("_wave[1:2]", "", predictor),effect_w1=effect) %>%
                          dplyr::select(outcome,predictor,control,effect_w1)
  boot_sh2 <- boot_sh %>% filter(grepl(pattern="wave2", predictor)) %>%
                          mutate(predictor= sub("_wave[1:2]", "", predictor),effect_w2=effect) %>%
                          dplyr::select(outcome,predictor,control,effect_w2)
  boot_2sh <- merge(boot_sh1,boot_sh2,by=c("outcome","predictor","control")) %>%
              mutate(effect_diff=effect_w2 - effect_w1)
  median <- median(boot_2sh$effect_diff)
  results_frame_diffs <- rbind(results_frame_diffs,  as.numeric(median))
}
results_frame_diffs <- results_frame_diffs[2:nrow(results_frame_diffs),]

obs_diff <- median(results_frame_wide_wave$wave1) - median(results_frame_wide_wave$wave2)
round(quantile(results_frame_diffs, probs = c(0.025, 0.975),na.rm=T), 3)
mean(results_frame_diffs <= obs_diff)
#observed difference = -0.027, bootstrapped p=.02, null models CI = -0.026 to 0.024


########################################################################################
# 6. create inferential statistics separately for each wave 
# Same approach as above but filter only wave 1 or only wave 2
#
########################################################################################

#adjusted fill_per function 
fill_per_w <- function(boot_results, frame, wave){
  for (m in 1:(bootstraps+1)){
    
    boot_results_m <- as.data.frame(boot_results[m])
    n <- 2
    if (wave=="wave1") {

      for (i in 1:13){      
        if(i == 1){
          results_subset <- boot_results_m
        } else if (i == 2) {
          results_subset <- boot_results_m %>% filter(grepl(pattern="PUBcomp.*wave1", predictor))
        } else if (i == 3) {
          results_subset <- boot_results_m %>% filter(grepl(pattern="ADRENcomp.*wave1", predictor))
        } else if (i == 4){
          results_subset <- boot_results_m %>% filter(grepl(pattern="GONADcomp.*wave1", predictor))
        } else if (i == 5){
          results_subset <- boot_results_m %>% filter(grepl(pattern="neg_pdsstage.*wave1", predictor))
        } else if (i == 6){
          results_subset <- boot_results_m %>% filter(grepl(pattern="parent_pdsstage.*wave1", predictor))
        } else if (i == 7){
          results_subset <- boot_results_m %>% filter(grepl(pattern="ldstage.*wave1", predictor))
        } else if (i == 8){
          results_subset <- boot_results_m %>% filter(grepl(pattern="aam", predictor))
        } else if (i == 9){
          results_subset <- boot_results_m %>% filter(grepl(pattern="^subj_timing.*wave1", predictor))
        } else if (i == 10){
          results_subset <- boot_results_m %>% filter(grepl(pattern="parent_subj_timing.*wave1", predictor))
        } else if (i == 11){
          results_subset <- boot_results_m %>% filter(grepl(pattern="DHEA.*wave1", predictor))
        } else if (i == 12){
          results_subset <- boot_results_m %>% filter(grepl(pattern="TEST.*wave1", predictor))
        } else if (i == 13){
          results_subset <- boot_results_m %>% filter(grepl(pattern="_EST.*wave1", predictor))
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
    } else if (wave=="wave2") {
        
        for (i in 1:13){      
          if(i == 1){
            results_subset <- boot_results_m
          } else if (i == 2) {
            results_subset <- boot_results_m %>% filter(grepl(pattern="PUBcomp.*wave2", predictor))
          } else if (i == 3) {
            results_subset <- boot_results_m %>% filter(grepl(pattern="ADRENcomp.*wave2", predictor))
          } else if (i == 4){
            results_subset <- boot_results_m %>% filter(grepl(pattern="GONADcomp.*wave2", predictor))
          } else if (i == 5){
            results_subset <- boot_results_m %>% filter(grepl(pattern="neg_pdsstage.*wave2", predictor))
          } else if (i == 6){
            results_subset <- boot_results_m %>% filter(grepl(pattern="parent_pdsstage.*wave2", predictor))
          } else if (i == 7){
            results_subset <- boot_results_m %>% filter(grepl(pattern="ldstage.*wave2", predictor))
          } else if (i == 8){
            results_subset <- boot_results_m %>% filter(grepl(pattern="aam", predictor))
          } else if (i == 9){
            results_subset <- boot_results_m %>% filter(grepl(pattern="^subj_timing.*wave2", predictor))
          } else if (i == 10){
            results_subset <- boot_results_m %>% filter(grepl(pattern="parent_subj_timing.*wave2", predictor))
          } else if (i == 11){
            results_subset <- boot_results_m %>% filter(grepl(pattern="DHEA.*wave2", predictor))
          } else if (i == 12){
            results_subset <- boot_results_m %>% filter(grepl(pattern="TEST.*wave2", predictor))
          } else if (i == 13){
            results_subset <- boot_results_m %>% filter(grepl(pattern="_EST.*wave2", predictor))
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
  }
  return(frame)
}

#adjusted function to get confidence intervals
populate_table_bounds_w <- function(boot, table, wave){
  
  results_frame_medians <- data.frame("aam_final" = NA, "parent_subj_timing" = NA,"resid_neg_ADRENcomp" = NA,
                                      "resid_neg_DHEA_cor" = NA, "resid_neg_EST_cor" = NA, "resid_neg_GONADcomp" = NA, 
                                      "resid_neg_ldstage" = NA, "resid_neg_parent_pdsstage" = NA, "resid_neg_pdsstage" = NA,
                                      "resid_neg_PUBcomp" = NA, "resid_neg_TEST_cor" = NA, "subj_timing" = NA) #needs to be alphabetical order 
  if (wave=="wave1") {
    for(b in 1:(length(boot)-1)){
      boot_sh <- boot[[b]] %>% 
        filter(grepl(pattern="wave1|aam", predictor)) %>% 
        mutate(predictor= sub("_wave[1:2]|_im_wave[1:2]$", "", predictor))  
      medians <- boot_sh %>% group_by(predictor) %>% dplyr::summarise(median = median(effect))
      results_frame_medians <- rbind(results_frame_medians, as.numeric(t(medians)[2,]))
    }
  } else {
      for(b in 1:(length(boot)-1)){
        boot_sh <- boot[[b]] %>% 
          filter(grepl(pattern="wave2|aam", predictor)) %>% 
          mutate(predictor= sub("_wave[1:2]|_im_wave[1:2]$", "", predictor))  
        medians <- boot_sh %>% group_by(predictor) %>% dplyr::summarise(median = median(effect))
        results_frame_medians <- rbind(results_frame_medians, as.numeric(t(medians)[2,]))
    }
  }    
  
  results_frame_medians <- results_frame_medians[2:nrow(results_frame_medians),]
  
  for(i in 1:nrow(table)){
    if((table[i,"sig_measure"]) == "median ES"){
      table[i, c("lower", "upper")] <- round(quantile(results_frame_medians[, paste0(table[i, "predictor"])], probs = c(0.025, 0.975),na.rm=T), 2)
    } else {
      table[i, c("lower", "upper")] <- NA
    }
  }
  
  return(table)
}

per_frame_w1 <- fill_per_w(bootstraps_full, permutation_frame, wave="wave1")
per_frame_w2 <- fill_per_w(bootstraps_full, permutation_frame, wave="wave2")
bootstrap_summary_w1 <- analyse_boot(per_frame_w1)
bootstrap_summary_w2 <- analyse_boot(per_frame_w2)

results_frame_wave1pred <- results_frame_forwide %>% filter(wave=="wave1"|predictor=="aam_final")
results_table_wave1 <- populate_table_observed(results_frame = results_frame_wave1pred, table = results_table2)
results_table_wave1 <- populate_table_bootstrap(bootstrap_summary_w1, results_table_wave1)
results_table_wave1 <- populate_table_bounds_w(boot = boot_sca, results_table_wave1, wave="wave1")
results_table_wave1$predictor_wave  <- "wave1"

results_frame_wave2pred <- results_frame_forwide %>% filter(wave=="wave2"|predictor=="aam_final")
results_table_wave2 <- populate_table_observed(results_frame = results_frame_wave2pred, table = results_table2)
results_table_wave2 <- populate_table_bootstrap(bootstrap_summary_w2, results_table_wave2)
results_table_wave2 <- populate_table_bounds_w(boot = boot_sca, results_table_wave2, wave="wave2")
results_table_wave2$predictor_wave  <- "wave2"

results_table_bywave <- rbind(results_table_wave1,results_table_wave2)
write.csv(results_table_bywave, file = paste0(cas_dir,"projects/W1_W2_pubertal_timing/table_bywave.csv"))
