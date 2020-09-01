#load packages and data 
require(specr)
require(dplyr)
require(plyr)
cas_dir="Y:/dsnlab/TAG/"
bootstraps <- 500
set.seed(88)

results_frame <- read.csv(file=paste0(cas_dir,"projects/W1_W2_pubertal_timing/3_1_SCA.csv"))
boot_sca <- read_rds(paste0(cas_dir,"projects/W1_W2_pubertal_timing/3_1_SCA_SCAboot.rds"))
bootstraps_full <- read_rds(paste0(cas_dir,"projects/W1_W2_pubertal_timing/3_1_boot.rds"))

#set up conf intervals
get_boot_results <- function(boot_sca){
  all_data_frames <- do.call("rbind", boot_sca[1:bootstraps])
  results_frameb <- all_data_frames %>% select(outcome, predictor, control, effect) %>%
    group_by(outcome, predictor, control) %>% 
    dplyr::summarise_all(funs(mean))
  
  results_frameb[, c("conf.low", "conf.high")] <- all_data_frames %>% 
    dplyr::group_by(outcome, predictor, control) %>%
    dplyr::summarise(`conf.low`=quantile(effect, probs=0.025,na.rm=T),
                     `conf.high`=quantile(effect, probs=0.975,na.rm=T)) %>% ungroup() %>% select(conf.low, conf.high)
  return(results_frameb)
}

sc_boot <- get_boot_results(boot_sca)
results_frame_sc <- left_join(results_frame, sc_boot[, c("outcome", "predictor", "control","conf.low", "conf.high")], by = c("outcome", "predictor", "control"))

#make a results frame with variable names collapsed over imp and wave
results_frame_curve <- results_frame_sc %>% 
  mutate(timepoint=ifelse(grepl("wave2", results_frame_sc$predictor),"cross-sectional","prospective"),
         imputation=ifelse(grepl("_im_", results_frame_sc$outcome),"yes","no"))

results_frame_curve$predictor <- sub("_im_wave[1:2]$", "", results_frame_curve$predictor)
results_frame_curve$predictor <- sub("_wave[1:2]$", "", results_frame_curve$predictor)
results_frame_curve$predictor <- revalue(results_frame_curve$predictor, 
                                         c("subj_timing"="l self-reported subjective timing",
                                         "parent_subj_timing"="k parent-reported subjective timing",
                                         "resid_neg_PUBcomp"="j residualized puberty composite",
                                        "resid_neg_ADRENcomp"="i residualized adrenal composite",
                                        "resid_neg_GONADcomp"="h residualized gonadal composite",
                                        "resid_neg_ldstage"="g residualized line drawings stage",
                                        "resid_neg_pdsstage"="f residualized self-reported PDS stage",
                                        "resid_neg_parent_pdsstage"="e residualized parent-reported PDS stage",
                                        "resid_neg_DHEA_cor"="d residualized DHEA level",
                                        "resid_neg_TEST_cor"="c residualized testosterone level",
                                        "resid_neg_EST_cor"="b residualized estradiol level",
                                        "aam_final"="a age at menarche"))
results_frame_curve$outcome <- sub("_im_wave2$", "", results_frame_curve$outcome)
results_frame_curve$outcome <- sub("_wave2$", "", results_frame_curve$outcome)
results_frame_curve$outcome <- as.factor(revalue(results_frame_curve$outcome, 
                                         c("CESDC_total"="g depressive symptoms",
                                           "SCARED_mean"="f anxiety symptoms",
                                           "depres_d"="d depressive disorder",
                                           "anx_d"="c anxiety disorder",
                                           "int_d"="e internalizing disorder",
                                           "distress_d"="b distress disorder",
                                           "fear_d"="a fear disorder")))
results_frame_curve$estimate <- results_frame_curve$effect
results_frame_curve$controls <- revalue(results_frame_curve$control, 
                                        c("mh"="time 1 psychopathology",
                                          "ctq"="early life stress",
                                          "both"="both controls",
                                          "none"="without controls"))

#Base the colors on p-value instead of conf.interval 
results_frame_curve <- results_frame_curve %>% mutate(color2 = ifelse(results_frame_curve$p_value<0.05, "#e41a1c", "darkgrey"))

# Plot specification curve
p1 <- plot_curve(results_frame_curve, ci=F, ribbon=T) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  ylim(-2.5, 2.5) + labs(y = "regression coefficient",size=.5)

p1[["data"]][["color"]] <- p1[["data"]][["color2"]] 

p2 <- plot_choices(results_frame_curve, choices = c("predictor", "outcome", "timepoint","controls")) +
  labs(x = "specifications (ranked)")   # element_text(angle = 0) 

p2[["data"]][["color"]] <- p2[["data"]][["color2"]] 

plot_specs(plot_a = p1, plot_b = p2, rel_height = c(1, 2))

