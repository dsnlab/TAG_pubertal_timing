#load packages and data 
require(specr)
require(dplyr)
require(plyr)
cas_dir="Y:/dsnlab/TAG/"
bootstraps <- 500

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
                                         c("subj_timing"="self-reported subjective timing",
                                         "parent_subj_timing"="parent-reported subjective timing",
                                         "resid_neg_PUBcomp"="residualized puberty composite",
                                        "resid_neg_ADRENcomp"="residualized adrenal composite",
                                        "resid_neg_GONADcomp"="residualized gonadal composite",
                                        "resid_neg_ldstage"="residualized line drawings stage",
                                        "resid_neg_pdsstage"="residualized self-reported PDS stage",
                                        "resid_neg_parent_pdsstage"="residualized parent-reported PDS stage",
                                        "resid_neg_DHEA_cor"="residualized DHEA level",
                                        "resid_neg_TEST_cor"="residualized testosterone level",
                                        "resid_neg_EST_cor"="residualized estradiol level",
                                        "aam_final"="age at menarche"))
results_frame_curve$outcome <- sub("_im_wave2$", "", results_frame_curve$outcome)
results_frame_curve$outcome <- sub("_wave2$", "", results_frame_curve$outcome)
results_frame_curve$outcome <- revalue(results_frame_curve$outcome, 
                                         c("CESDC_total"="depressive symptoms",
                                           "SCARED_mean"="anxiety symptoms",
                                           "depres_d"="depressive disorder",
                                           "anx_d"="anxiety disorder",
                                           "int_d"="internalizing disorder",
                                           "distress_d"="distress disorder",
                                           "fear_d"="fear disorder"))
results_frame_curve$estimate <- results_frame_curve$effect
results_frame_curve$controls <- revalue(results_frame_curve$control, 
                                        c("mh"="time 1 psychopathology",
                                          "ctq"="early life stress",
                                          "both"="both controls",
                                          "none"="without controls"))


# Plot specification curve
p1 <- plot_curve(results_frame_curve, ci=F, ribbon=T) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  ylim(-2.5, 2.5) 

p2 <- plot_choices(results_frame_curve, choices = c("predictor", "outcome", "timepoint","controls")) +
  labs(x = "specifications (ranked)") + theme(strip.text.y.right = element_blank()) # element_text(angle = 0)


plot_specs(plot_a = p1, plot_b = p2, rel_height = c(1, 2))

