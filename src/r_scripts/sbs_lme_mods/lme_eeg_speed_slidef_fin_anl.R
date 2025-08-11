
knitr::opts_chunk$set(echo = TRUE)
# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console
cat("\014")
# Clean workspace
rm(list=ls())

# Packages & Setup
# install.packages(c("tidyverse","purrr","R.matlab","readxl","dplyr"))
#%% PACKAGES FOR STATS
library(readxl);
library(purrr);
library(tidyverse);
library(tibble);
library(knitr);
library(gtsummary);
library(kableExtra);
library(lme4);
library(MuMIn);
library(car);
library(effectsize)
library(sjPlot);
library(emmeans);
library(knitr)
library(lmerTest)
#%% PACKAGES FOR PLOTS & HTML HANDLING
# library(effects);
# library(sjPlot);
# library(plotly);
# library(webshot)
# library(reshape2);
# library(htmltools)
# library(Polychrome);
# library(htmlwidgets);
# library(shiny)
# library(webshot)
library(circlize)
library(purrr)
library(scatterplot3d)
library(RColorBrewer)
library(openxlsx)
library(kableExtra) # devtools::install_github("haozhu233/kableExtra")
library(R.matlab)

## GTSUMMARY THEME
#%% ANOVA OPTIONS
options(contrasts = c("contr.sum","contr.poly")) # suggested for type III
gtsummary::set_gtsummary_theme(theme_gtsummary_journal("jama"))

#%% FIGURE & TABLE OPTIONS
calc_cohensf2 <- function(mod_main,mod_alt){
  r2_out = r.squaredGLMM(mod_main);
  r2_outalt = r.squaredGLMM(mod_alt);
  r2m = r2_out[1] # input your R2
  f2m = r2m/(1 - r2m)
  r2c = r2_out[2] # input your R2
  f2c = r2c/(1 - r2c)
  f2m = (r2_out[1]-r2_outalt[1])/(1-r2_out[1]);
  f2c = (r2_out[2]-r2_outalt[2])/(1-r2_out[2]);
  print(str_glue("r2m: {round(r2m,4)},\tr2c: {round(r2c,4)}\n\n"))
  print(str_glue("f2m: {round(f2m,4)},\tf2c: {round(f2c,4)}\n\n"))
  vals = data.frame(r2m, r2c, f2m, f2c);
  return (vals)
}
#--
ispc <- function() {
  sys_name <- Sys.info()["sysname"]
  if (sys_name == "Windows") {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#%% CLUSTERS TO PLOT
#%% CLUSTERS TO PLOT
clusters = c(3,4,5,6,7,8,9,10,11,12,13) # RSup/RSM, PreC, LSM, Mid Cing, LSup, LPPA, RPPA
cluster_chars = c('RPP','RSM','LPreC','LSM','RPre','LPP','LSupM','ROcc','MCing','LTemp','LOcc')
#%% KIN PARAMS
#--
# kin_measures = c('mu_ml_exc_mm_gc','mu_step_dur',
#                  'cov_ml_exc_mm_gc','cov_step_dur_mm_gc',
#                  'mu_double_sup_dur','mu_single_sup_dur',
#                  'cov_double_sup_dur','cov_single_sup_dur',
#                  'mu_stance_dur','mu_swing_dur',
#                  'cov_stance_dur','cov_swing_dur');
# kin_title_chars = c('mu_ml_exc_mm_gc','mu_step_dur',
#                     'cov_ml_exc_mm_gc','cov_step_dur_mm_gc',
#                     'mu_double_sup_dur','mu_single_sup_dur',
#                     'cov_double_sup_dur','cov_single_sup_dur',
#                     'mu_stance_dur','mu_swing_dur',
#                     'cov_stance_dur','cov_swing_dur')
kin_measures = c('mu_ml_exc_mm_gc_fn1','mu_step_dur_fn1',
                 'cov_ml_exc_mm_gc_fn1','cov_step_dur_fn1',
                 'mu_double_sup_dur_fn1','mu_single_sup_dur_fn1',
                 'cov_double_sup_dur_fn1','cov_single_sup_dur_fn1',
                 'mu_stance_dur_fn1','mu_swing_dur_fn1',
                 'cov_stance_dur_fn1','cov_swing_dur_fn1',
                 'std_ml_exc_mm_gc_fn1','std_step_dur_fn1');
kin_title_chars = c('mu_ml_exc_mm_gc','mu_step_dur',
                    'cov_ml_exc_mm_gc','cov_step_dur',
                    'mu_double_sup_dur','mu_single_sup_dur',
                    'cov_double_sup_dur','cov_single_sup_dur',
                    'mu_stance_dur','mu_swing_dur',
                    'cov_stance_dur','cov_swing_dur')
#%% EEG PARAMS
#--
# eeg_measures = c('mu_avg_theta','mu_avg_alpha','mu_avg_beta',
#                  'std_avg_theta','std_avg_alpha','std_avg_beta',
#                  'cov_i_avg_theta','cov_i_avg_alpha','cov_i_avg_beta');
# eeg_title_chars = c("**THETA** Mean","**ALPHA** Mean","**BETA** Mean",
#                     "**THETA** Std. Dev.","**ALPHA** Std. Dev.","**BETA** Std. Dev.",
                    # "**THETA** COVi","**ALPHA** COVi","**BETA** COVi");
#--
eeg_measures = c('mu_avg_theta_fn1','mu_avg_alpha_fn1','mu_avg_beta_fn1',
                 'std_avg_theta_fn1','std_avg_alpha_fn1','std_avg_beta_fn1',
                 'cov_avg_theta_calc_fn1','cov_avg_alpha_calc_fn1','cov_avg_beta_calc_fn1');
eeg_title_chars = c("**THETA** Mean","**ALPHA** Mean","**BETA** Mean",
                    "**THETA** Std. Dev.","**ALPHA** Std. Dev.","**BETA** Std. Dev.",
                    "**THETA** COVi","**ALPHA** COVi","**BETA** COVi");

#%% LOADING PARAMETERS
#--
FEXT_KINT = 'allcond_perstridefb_apfix_std_mi_nfslidingb36';

#%% ======================================================================== %%#
#%% LOAD
EXCEL_FNAME = sprintf("%s_lmers.xlsx",FEXT_KINT)
DATA_FNAME = sprintf("sbs_eeg_psd_kin_%s.xlsx",FEXT_KINT)
#--
dpath = paste0("/jsalminen/GitHub/MIND_IN_MOTION_PRJ/_data/MIM_dataset/_studies/02202025_mim_yaoa_powpow0p3_crit_speed/__iclabel_cluster_allcond_rb3/icrej_5/11/kin_eeg_step_to_step");
excel_dir <- file.path(dpath,DATA_FNAME);

#--
if(ispc()){
  excel_dir <- paste0("M:",excel_dir)
}else{
  excel_dir <- paste0("/blue/dferris",excel_dir);
}
orig_eegt <- read_excel(excel_dir,sheet="Sheet1")

#%% SAVE DIR
if(ispc()){
  save_dir <- paste0("M:",dpath)
}else{
  save_dir <- paste0("/blue/dferris",dpath);
}

#%% SUBSET TABLE
orig_eegt <- orig_eegt %>%
  select(subj_char,cond_char,speed_n,group_char,cluster_n,model_n,group_n,
         std_avg_theta,std_avg_beta,std_avg_alpha,
         std_ml_exc_mm_gc,std_ap_exc_mm_gc,std_ud_exc_mm_gc,
         std_stance_dur,std_step_dur,std_single_sup_dur,std_double_sup_dur,std_swing_dur,
         mu_avg_theta,mu_avg_beta,mu_avg_alpha,
         mu_ml_exc_mm_gc,mu_ap_exc_mm_gc,mu_ud_exc_mm_gc,
         mu_stance_dur,mu_step_dur,mu_single_sup_dur,mu_double_sup_dur,mu_swing_dur,
         cov_i_avg_theta,cov_i_avg_alpha,cov_i_avg_beta,
         cov_avg_theta,cov_avg_beta,cov_avg_alpha,
         cov_ml_exc_mm_gc,cov_ap_exc_mm_gc,cov_ud_exc_mm_gc,
         cov_stance_dur,cov_step_dur,cov_single_sup_dur,cov_double_sup_dur,cov_swing_dur);

orig_eegt <- transform(orig_eegt,cov_avg_theta_calc = std_avg_theta/mu_avg_theta)
orig_eegt <- transform(orig_eegt,cov_avg_alpha_calc = std_avg_alpha/mu_avg_alpha)
orig_eegt <- transform(orig_eegt,cov_avg_beta_calc = std_avg_beta/mu_avg_beta)
# orig_eegt <- transform(orig_eegt,cov_i_avg_theta = log(cov_i_avg_theta))
# orig_eegt <- transform(orig_eegt,cov_i_avg_alpha = log(cov_i_avg_alpha))
# orig_eegt <- transform(orig_eegt,cov_i_avg_beta = log(cov_i_avg_beta))
# orig_eegt <- transform(orig_eegt,qcv_i_avg_theta = log10(abs(qcv_i_avg_theta)))
# orig_eegt <- transform(orig_eegt,qcv_i_avg_alpha = log10(abs(qcv_i_avg_alpha)))
# orig_eegt <- transform(orig_eegt,qcv_i_avg_beta = log10(abs(qcv_i_avg_beta)))
#--
orig_eegt <- transform(orig_eegt,cov_ml_exc_mm_gc = 100*cov_ml_exc_mm_gc)
orig_eegt <- transform(orig_eegt,cov_step_dur = 100*cov_step_dur)
orig_eegt <- transform(orig_eegt,cov_stance_dur = 100*cov_stance_dur)
orig_eegt <- transform(orig_eegt,cov_single_sup_dur = 100*cov_single_sup_dur)
orig_eegt <- transform(orig_eegt,cov_double_sup_dur = 100*cov_double_sup_dur)
orig_eegt <- transform(orig_eegt,cov_swing_dur = 100*cov_swing_dur)
#--
dtbl <- orig_eegt;

#%% CREATE MEAN & STD MEASURES PER CONDITION
dtbl <- filter_at(dtbl,vars('cond_char'), any_vars(. %in% c('0p25','0p5','0p75','1p0')))
dtbl <- dtbl %>%
  group_by(subj_char,cond_char,speed_n,model_n,group_char,cluster_n,group_n) %>%
  summarise_all(list("mean","sd"))

#%% CHANGE GROUP LABS
dtbl$speed_cond_num <- as.numeric(dtbl$speed_n);
dtbl$subj_char <- as.factor(dtbl$subj_char);
dtbl$group_char <- as.factor(dtbl$group_char);
dtbl$group_n <- as.factor(dtbl$group_n); 
dtbl$model_n <- as.factor(dtbl$model_n);

#%% MUTATE VARIABLES
dtbl$speed_cond_num <- as.numeric(dtbl$speed_n);
write.xlsx(dtbl,file.path(save_dir,paste0("lek_musd_",FEXT_KINT,"_tbl.xlsx")))

#%% EXCEL DATAFRAME
excel_df <- data.frame(cluster_n=double(),
                       group_char=character(),
                       model_char=character(),
                       kinematic_char=character(),
                       freq_band_char=character(),
                       mod_num_obs=character(),
                       coeff_chars=character(),
                       coeffs=character(),
                       coeff_desm=character(),
                       lme_pval=character(),
                       lme_est=character(),
                       lme_stat=character(),
                       lme_cc=character(),
                       coeff_m=character(),
                       coeff_b=character(),
                       pwc_cc=character(),
                       pwc_est=character(),
                       pwc_pval=character(),
                       confint_chars=character(),
                       emmeans=character(),
                       emmeans_se=character(),
                       confint_lwr=character(),
                       confint_upr=character(),
                       ttn0_cc=character(),
                       ttn0_est=character(),
                       ttn0_pval=character(),
                       ttna0_cc=character(),
                       ttna0_est=character(),
                       ttna0_pval=character(),
                       emm_slp_cc=character(),
                       emm_slp_est=character(),
                       emm_slp_se=character(),
                       emm_slp_df=character(),
                       emm_slpcon_cc=character(),
                       emm_slpcon_est=character(),
                       emm_slpcon_pval=character(),
                       emm_slpcon_se=character(),
                       emm_slpcon_trat=character(),
                       emm_slptt_cc=character(),
                       emm_slptt_est=character(),
                       emm_slptt_pval=character(),
                       emm_slptt_se=character(),
                       emm_slptt_trat=character(),
                       anv_chars=character(),
                       anv_pvals=character(),
                       anv_stats=character(),
                       anv_dfs=character(),
                       r2_m_int=double(),
                       r2_c_int=double(),
                       f2_m_int=double(),
                       f2_c_int=double(),
                       fsq_chars=character(),
                       fsq_vals=character(),
                       etasq_chars=character(),
                       etasq_vals=character(),
                       ran_effs_char=character(),
                       ran_effs_n=character())


#%% ======================================================================== %%#
#%% HTML TAG
cat(paste0("\n\n# Interaction Models\n"))
# SPEED INTERACITON
theme_set(theme_classic()) # This sets the default ggplot theme
#%% LOOP VARS
tmp <- length(clusters)
cis <- rep(NA,tmp*length(eeg_measures))
eis <- rep(NA,tmp*length(eeg_measures))
cnt = 1;
for (i in 1:length(clusters)) {
  for(k in 1:length(eeg_measures)){
    cis[cnt] = i;
    eis[cnt] = k;
    cnt = cnt + 1;
  }
}

for (i in 1:length(cis)){
  #%% SET LOOP PARAMS
  ci = clusters[cis[i]];
  ei = eeg_measures[eis[i]];
  
  #%% HTML TAG
  cat(paste0("\n\n## ",ei,',',cluster_chars[cis[i]],"\n"))
  
  #%% SUBSET BY CLUSTER
  tmpt <- filter_at(dtbl,vars('cluster_n'), any_vars(. %in% ci));
  # -- remove NAN cases
  tmpt = tmpt[complete.cases(tmpt[[ei]]),];
  
  #%% COMPUTE LME MODEL
  mod_char = paste(ei," ~ 1 + speed_n + group_char + speed_n:group_char + (1|subj_char)");
  # mod <- lme4::lmer(str_glue(mod_char), data=tmpt);
  mod <- lmerTest::lmer(str_glue(mod_char), data=tmpt);
  sm = summary(mod)
  
  #%% PRINT TABLE
  # tbl <- mod %>%
  #   tbl_regression(tidy_fun = broom.mixed::tidy,
  #                  pvalue_fun = purrr::partial(style_pvalue, digits = 2),
  #                  intercept=TRUE) %>%
  #   add_global_p() %>%
  #   add_q()
  # 
  # t_out <- as_gt(tbl) %>%
  #   gt::tab_header(title=c(ei,") for Connection ",ci,"-",cj)) %>%
  #   gt::tab_options(table.layout = "fixed") %>%
  #   gt::as_raw_html()
  # print(t_out)
  
  # #-- extract fixed & random effects coeffs
  anv_vals <- car::Anova(mod,type=3)
  fix_effs <- fixef(mod)
  ran_effs <- ranef(mod)$subj_char
  
  #%% EFFECT SIZES
  cat(paste0("\n\n### model effect sizes\n"))
  #-- Intercept Model for Cohen's f^2
  print('intercept model\n');
  mod_char_int = paste(ei," ~ 1 + (1|subj_char)");
  mod <- lmerTest::lmer(str_glue(mod_char), data=tmpt);
  sfit <- summary(mod) 
  lme_est <- as.numeric(sfit$coefficients[, "Estimate"])
  lme_pval <- as.numeric(sfit$coefficients[, "Pr(>|t|)"])
  lme_stat <- as.numeric(sfit$coefficients[,"t value"])
  lme_cc <- rownames(sfit$coefficients[,])
  #-- extract fixed & random effects coeffs
  # anv_vals <- car::Anova(mod,type=3)
  # anv_chars=paste0(unlist(row.names(anv_vals)),collapse=',')
  # anv_pvals=paste0(unlist(anv_vals$`Pr(>Chisq)`),collapse=',')
  # anv_stats=paste0(unlist(anv_vals$Chisq),collapse=',')
  # anv_dfs=paste0(unlist(anv_vals$Df),collapse=',')
  #--
  anv_vals <- anova(mod)
  anv_chars=paste0(unlist(row.names(anv_vals)),collapse=',')
  anv_pvals=paste0(unlist(anv_vals$`Pr(>F)`),collapse=',')
  anv_stats=paste0(unlist(anv_vals$`F value`),collapse=',')
  anv_dfs=paste0(unlist(anv_vals$DenDF),collapse=',')
  
  fix_effs <- fixef(mod)
  ran_effs <- ranef(mod)$subj_char
  
  #%% EFFECT SIZES
  cat(paste0("\n\n### model effect sizes\n"))
  #-- Intercept Model for Cohen's f^2
  print('intercept model\n');
  mod_char_int = paste0(ei," ~ 1 + (1|subj_char)");
  tmp_mod = lme4::lmer(str_glue(mod_char_int), data=tmpt)
  cf_alt = calc_cohensf2(mod, tmp_mod)
  #-- calculate effect sizes using r-package
  anv_vals_f <- car::Anova(mod,type=3,test.statistic = "F")
  etasq_res <- effectsize::eta_squared(anv_vals_f, partial = TRUE)
  cohens_fsq_res <- effectsize::cohens_f_squared(anv_vals_f, partial = TRUE)
  #-- confidence intervals
  summ <- summary(mod);
  coeff_chars = row.names(summ$coefficients)
  emm = emmeans::emmeans(mod,spec=c('speed_n','group_char'),level=0.95);
  bfc = emmeans::contrast(emm,adjutst="bonferroni")
  pwc = emmeans::contrast(emm,method="pairwise")
  ttz = data.frame(emmeans::test(emm,null=0))
  #--
  emma = emmeans::emmeans(mod,spec=c('speed_n'),level=0.95);
  ttza = data.frame(emmeans::test(emma,null=0))
  #-- design mat
  desm = unique(mod@pp[[".->X"]]);
  desm = format_csv(data.frame(desm));
  #--
  cfis <- data.frame(emm)
  pwcs = data.frame(pwc)
  #--
  emms = emtrends(mod, pairwise ~ group_char, var = "speed_n")
  emm_slp = data.frame(emms$emtrends)
  emm_slpcon = data.frame(emms$contrasts)
  emm_slptt = test(emtrends(mod, specs = ~ group_char, var = "speed_n"), side = "two_sided")
  emm_slpcontt = data.frame(emm_slptt)
  # print(kable(cfis));
  # print(kable(anv_vals));
  # print(kable(data.frame(bfc)));
  # print(kable(data.frame(pwc)));
  
  #%% EXCEL PRINTS
  new_row = data.frame(cluster_n=ci,
                       group_char=c('all'),
                       model_char=c('speed_group_intact_all'),
                       kinematic_char='none',
                       freq_band_char=c(ei),
                       mod_num_obs=nrow(tmpt),
                       coeff_chars=paste0(unlist(coeff_chars),collapse=','),
                       coeffs=paste0(unlist(mod@beta),collapse=','),
                       coeff_desm=desm,
                       lme_pval=paste0(unlist(lme_pval),collapse=','),
                       lme_est=paste0(unlist(lme_est),collapse=','),
                       lme_stat=paste0(unlist(lme_stat),collapse=','),
                       lme_cc=paste0(unlist(lme_cc),collapse=','),
                       pwc_cc=paste0(unlist(pwcs$contrast),collapse=','),
                       pwc_est=paste0(unlist(pwcs$estimate),collapse=','),
                       pwc_pval=paste0(unlist(pwcs$p.value),collapse=','),
                       confint_chars=paste0(unlist(cfis$group_char),collapse=','),
                       emmeans=paste0(unlist(as.numeric(cfis$emmean)),collapse=','),
                       emmeans_se=paste0(unlist(as.numeric(cfis$SE)),collapse=','),
                       confint_lwr=paste0(unlist(as.numeric(cfis$lower.CL)),collapse=','),
                       confint_upr=paste0(unlist(as.numeric(cfis$upper.CL)),collapse=','),
                       ttn0_cc=paste0(unlist(as.numeric(ttz$group_char)),collapse=','),
                       ttn0_est=paste0(unlist(as.numeric(ttz$emmean)),collapse=','),
                       ttn0_pval=paste0(unlist(as.numeric(ttz$p.value)),collapse=','),
                       ttna0_cc="all",
                       ttna0_est=paste0(unlist(as.numeric(ttza$emmean)),collapse=','),
                       ttna0_pval=paste0(unlist(as.numeric(ttza$p.value)),collapse=','),
                       emm_slp_cc=paste0(unlist(emm_slp$group_char),collapse=','),
                       emm_slp_est=paste0(unlist(as.numeric(emm_slp[,'speed_n.trend'])),collapse=','),
                       emm_slp_se=paste0(unlist(as.numeric(emm_slp$SE)),collapse=','),
                       emm_slp_df=paste0(unlist(as.numeric(emm_slp$df)),collapse=','),
                       emm_slpcon_cc=paste0(unlist(emm_slpcon$contrast),collapse=','),
                       emm_slpcon_est=paste0(unlist(as.numeric(emm_slpcon$estimate)),collapse=','),
                       emm_slpcon_pval=paste0(unlist(as.numeric(emm_slpcon$p.value)),collapse=','),
                       emm_slpcon_se=paste0(unlist(as.numeric(emm_slpcon$SE)),collapse=','),
                       emm_slpcon_trat=paste0(unlist(as.numeric(emm_slpcon$t.ratio)),collapse=','),
                       emm_slptt_cc=paste0(unlist(emm_slpcontt$group_char),collapse=','),
                       emm_slptt_est=paste0(unlist(as.numeric(emm_slpcontt[,'speed_n.trend'])),collapse=','),
                       emm_slptt_pval=paste0(unlist(as.numeric(emm_slpcontt[,'p.value'])),collapse=','),
                       emm_slptt_se=paste0(unlist(as.numeric(emm_slpcontt[,'SE'])),collapse=','),
                       emm_slptt_trat=paste0(unlist(as.numeric(emm_slpcontt[,'t.ratio'])),collapse=','),
                       anv_chars=anv_chars,
                       anv_pvals=anv_pvals,
                       anv_stats=anv_stats,
                       anv_dfs=anv_dfs,
                       r2_m_int=cf_alt$r2m,
                       r2_c_int=cf_alt$r2c,
                       f2_m_int=cf_alt$f2m,
                       f2_c_int=cf_alt$f2c,
                       fsq_chars=paste0(unlist(cohens_fsq_res$Parameter),collapse=','),
                       fsq_vals=paste0(unlist(cohens_fsq_res$Cohens_f2_partial),collapse=','),
                       etasq_chars=paste0(unlist(etasq_res$Parameter),collapse=','),
                       etasq_vals=paste0(unlist(etasq_res$Eta2_partial),collapse=','),
                       ran_effs_char=paste0(unlist(row.names(ran_effs)),collapse=','),
                       ran_effs_n=paste0(unlist(ran_effs$`(Intercept)`),collapse=','))
  excel_df <- rbind(excel_df,new_row)
  
  # #%% MODEL VALIDATION PLOTS
  # cat(paste0("\n\n### model validations\n"))
  # print(plot_model(mod, type = 'diag'))
  # cat("\n")
  #
  # tmpt$fit.m <- predict(mod, re.form = NA)
  # tmpt$fit.c <- predict(mod, re.form = NULL)
  # 
  # #%% VISUALIZATION OF MODELS
  # cat(paste0("\n\n### data plot\n"))
  # theme_set(theme_classic()) # This sets the default ggplot theme
  # print(
  #   tmpt %>%
  #     ggplot() +
  #     geom_jitter(aes(x = speed_n, y = .data[[ei]], color = subj_char), pch=20, size=10, alpha=0.2) +
  #     geom_line(aes(x= speed_n, y=fit.c, group = subj_char),linetype = "dashed", linewidth = .5) +
  #     geom_line(aes(x= speed_n, y=fit.m), linewidth = 2)+
  #     geom_errorbar(data=cfis,aes(x = speed_n, ymin=lower.CL, ymax=upper.CL), width=0.33, linewidth=1.5, colour='red') +
  #     facet_wrap(~group_char) +
  #     ggtitle(paste0('Connection ',ci,'-',cj))+
  #     xlab("Treadmill Speed (m/s)") +
  #     ylab(eeg_title_chars[eis[i]]) +
  #     guides(group="none")
  # )
}


#%% HTML TAG
cat(paste0("\n\n# Group Models\n"))
# SPEED GROUP
theme_set(theme_classic()) # This sets the default ggplot theme
#%% LOOP VARS
tmp <- length(clusters)
cis <- rep(NA,tmp*length(eeg_measures))
eis <- rep(NA,tmp*length(eeg_measures))
cnt = 1;
for (i in 1:length(clusters)) {
  for(k in 1:length(eeg_measures)){
    cis[cnt] = i;
    eis[cnt] = k;
    cnt = cnt + 1;
  }
}

for (i in 1:length(cis)){
  #%% SET LOOP PARAMS
  ci = clusters[cis[i]];
  ei = eeg_measures[eis[i]];
  
  #%% HTML TAG
  cat(paste0("\n\n## ",ei,',',cluster_chars[cis[i]],"\n"))
  
  #%% SUBSET BY CLUSTER
  tmpt <- filter_at(dtbl,vars('cluster_n'), any_vars(. %in% ci));
  # -- remove NAN cases
  tmpt = tmpt[complete.cases(tmpt[[ei]]),];
  
  #%% COMPUTE LME MODEL
  mod_char = paste(ei," ~ 1 + speed_n + group_char + (1|subj_char)");
  # mod <- lme4::lmer(str_glue(mod_char), data=tmpt);
  mod <- lmerTest::lmer(str_glue(mod_char), data=tmpt);
  sfit <- summary(mod)
  
  #%% PRINT TABLE
  # tbl <- mod %>%
  #   tbl_regression(tidy_fun = broom.mixed::tidy,
  #                  pvalue_fun = purrr::partial(style_pvalue, digits = 2),
  #                  intercept=TRUE) %>%
  #   add_global_p() %>%
  #   add_q()
  # 
  # t_out <- as_gt(tbl) %>%
  #   gt::tab_header(title=c(ei,") for Connection ",ci,"-",cj)) %>%
  #   gt::tab_options(table.layout = "fixed") %>%
  #   gt::as_raw_html()
  # print(t_out)
  
  lme_est <- as.numeric(sfit$coefficients[, "Estimate"])
  lme_pval <- as.numeric(sfit$coefficients[, "Pr(>|t|)"])
  lme_stat <- as.numeric(sfit$coefficients[,"t value"])
  lme_cc <- rownames(sfit$coefficients[,])
  #-- extract fixed & random effects coeffs
  # anv_vals <- car::Anova(mod,type=3)
  # anv_chars=paste0(unlist(row.names(anv_vals)),collapse=',')
  # anv_pvals=paste0(unlist(anv_vals$`Pr(>Chisq)`),collapse=',')
  # anv_stats=paste0(unlist(anv_vals$Chisq),collapse=',')
  # anv_dfs=paste0(unlist(anv_vals$Df),collapse=',')
  #--
  anv_vals <- anova(mod)
  anv_chars=paste0(unlist(row.names(anv_vals)),collapse=',')
  anv_pvals=paste0(unlist(anv_vals$`Pr(>F)`),collapse=',')
  anv_stats=paste0(unlist(anv_vals$`F value`),collapse=',')
  anv_dfs=paste0(unlist(anv_vals$DenDF),collapse=',')
  
  fix_effs <- fixef(mod)
  ran_effs <- ranef(mod)$subj_char
  
  #%% EFFECT SIZES
  cat(paste0("\n\n### model effect sizes\n"))
  #-- Intercept Model for Cohen's f^2
  print('intercept model\n');
  mod_char_int = paste0(ei," ~ 1 + (1|subj_char)");
  tmp_mod = lme4::lmer(str_glue(mod_char_int), data=tmpt)
  cf_alt = calc_cohensf2(mod, tmp_mod)
  #-- calculate effect sizes using r-package
  anv_vals_f <- car::Anova(mod,type=3,test.statistic = "F")
  etasq_res <- effectsize::eta_squared(anv_vals_f, partial = TRUE)
  cohens_fsq_res <- effectsize::cohens_f_squared(anv_vals_f, partial = TRUE)
  #-- confidence intervals
  summ <- summary(mod);
  coeff_chars = row.names(summ$coefficients)
  emm = emmeans::emmeans(mod,spec=c('speed_n','group_char'),level=0.95);
  bfc = emmeans::contrast(emm,adjutst="bonferroni")
  pwc = emmeans::contrast(emm,method="pairwise")
  ttz = data.frame(emmeans::test(emm,null=0))
  #--
  emma = emmeans::emmeans(mod,spec=c('speed_n'),level=0.95);
  ttza = data.frame(emmeans::test(emma,null=0))
  #-- design mat
  desm = unique(mod@pp[[".->X"]]);
  desm = format_csv(data.frame(desm));
  #--
  cfis <- data.frame(emm)
  pwcs = data.frame(pwc)
  #--
  # print(kable(cfis));
  # print(kable(anv_vals));
  # print(kable(data.frame(bfc)));
  # print(kable(data.frame(pwc)));
  
  #%% EXCEL PRINTS
  new_row = data.frame(cluster_n=ci,
                       group_char=c('all'),
                       model_char=c('speed_group_all'),
                       kinematic_char='none',
                       freq_band_char=c(ei),
                       mod_num_obs=nrow(tmpt),
                       coeff_chars=paste0(unlist(coeff_chars),collapse=','),
                       coeffs=paste0(unlist(mod@beta),collapse=','),
                       coeff_desm=desm,
                       lme_pval=paste0(unlist(lme_pval),collapse=','),
                       lme_est=paste0(unlist(lme_est),collapse=','),
                       lme_stat=paste0(unlist(lme_stat),collapse=','),
                       lme_cc=paste0(unlist(lme_cc),collapse=','),
                       pwc_cc=paste0(unlist(pwcs$contrast),collapse=','),
                       pwc_est=paste0(unlist(pwcs$estimate),collapse=','),
                       pwc_pval=paste0(unlist(pwcs$p.value),collapse=','),
                       confint_chars=paste0(unlist(cfis$group_char),collapse=','),
                       emmeans=paste0(unlist(as.numeric(cfis$emmean)),collapse=','),
                       emmeans_se=paste0(unlist(as.numeric(cfis$SE)),collapse=','),
                       confint_lwr=paste0(unlist(as.numeric(cfis$lower.CL)),collapse=','),
                       confint_upr=paste0(unlist(as.numeric(cfis$upper.CL)),collapse=','),
                       ttn0_cc=paste0(unlist(as.numeric(ttz$group_char)),collapse=','),
                       ttn0_est=paste0(unlist(as.numeric(ttz$emmean)),collapse=','),
                       ttn0_pval=paste0(unlist(as.numeric(ttz$p.value)),collapse=','),
                       ttna0_cc="all",
                       ttna0_est=paste0(unlist(as.numeric(ttza$emmean)),collapse=','),
                       ttna0_pval=paste0(unlist(as.numeric(ttza$p.value)),collapse=','),
                       emm_slp_cc="none",
                       emm_slp_est="none",
                       emm_slp_se="none",
                       emm_slp_df="none",
                       emm_slpcon_cc="none",
                       emm_slpcon_est="none",
                       emm_slpcon_pval="none",
                       emm_slpcon_se="none",
                       emm_slpcon_trat="none",
                       emm_slptt_cc="none",
                       emm_slptt_est="none",
                       emm_slptt_pval="none",
                       emm_slptt_se="none",
                       emm_slptt_trat="none",
                       anv_chars=anv_chars,
                       anv_pvals=anv_pvals,
                       anv_stats=anv_stats,
                       anv_dfs=anv_dfs,
                       r2_m_int=cf_alt$r2m,
                       r2_c_int=cf_alt$r2c,
                       f2_m_int=cf_alt$f2m,
                       f2_c_int=cf_alt$f2c,
                       fsq_chars=paste0(unlist(cohens_fsq_res$Parameter),collapse=','),
                       fsq_vals=paste0(unlist(cohens_fsq_res$Cohens_f2_partial),collapse=','),
                       etasq_chars=paste0(unlist(etasq_res$Parameter),collapse=','),
                       etasq_vals=paste0(unlist(etasq_res$Eta2_partial),collapse=','),
                       ran_effs_char=paste0(unlist(row.names(ran_effs)),collapse=','),
                       ran_effs_n=paste0(unlist(ran_effs$`(Intercept)`),collapse=','))
  excel_df <- rbind(excel_df,new_row)
  
  # #%% MODEL VALIDATION PLOTS
  # cat(paste0("\n\n### model validations\n"))
  # print(plot_model(mod, type = 'diag'))
  # cat("\n")
  # 
  # tmpt$fit.m <- predict(mod, re.form = NA)
  # tmpt$fit.c <- predict(mod, re.form = NULL)
  # 
  # #%% VISUALIZATION OF MODELS
  # cat(paste0("\n\n### data plot\n"))
  # theme_set(theme_classic()) # This sets the default ggplot theme
  # print(
  #   tmpt %>%
  #     ggplot() +
  #     geom_jitter(aes(x = speed_n, y = .data[[ei]], color = subj_char), pch=20, size=10, alpha=0.2) +
  #     geom_line(aes(x= speed_n, y=fit.c, group = subj_char),linetype = "dashed", linewidth = .5) +
  #     geom_line(aes(x= speed_n, y=fit.m), linewidth = 2)+
  #     geom_errorbar(data=cfis,aes(x = speed_n, ymin=lower.CL, ymax=upper.CL), width=0.33, linewidth=1.5, colour='red') +
  #     facet_wrap(~group_char) +
  #     ggtitle(paste0('Connection from',ci,'->',cj))+
  #     xlab("Treadmill Speed (m/s)") +
  #     ylab(eeg_title_chars[eis[i]]) +
  #     guides(group="none")
  # )
}

#%% ======================================================================== %%#
# KINEMATIC CALCULATIONS
#%% HTML TAG
cat(paste0("\n\n# Kin Interaction Models\n"))
# SPEED INTERACITON
theme_set(theme_classic()) # This sets the default ggplot theme
#%% LOOP VARS
kis <- rep(NA,length(kin_measures))
cnt = 1;
for (i in 1:length(kin_measures)) {
  kis[cnt] = i;
  cnt = cnt + 1;
}

for (i in 1:length(kis)){
  #%% SET LOOP PARAMS
  ki = kin_measures[kis[i]];
  
  #%% HTML TAG
  cat(paste0("\n\n## ", ki,"\n"))
  
  #%% CREATE MEAN & STD MEASURES PER CONDITION
  tmpt = dtbl;
  tmpt <- tmpt %>%
    group_by(subj_char,cond_char,speed_n,model_n,group_char,group_n) %>%
    summarise_all(list("mean"))
  
  #%% COMPUTE LME MODEL
  mod_char = paste(ki," ~ 1 + speed_n + group_char + speed_n:group_char + (1|subj_char)");
  # mod <- lme4::lmer(str_glue(mod_char), data=tmpt);
  mod <- lmerTest::lmer(str_glue(mod_char), data=tmpt);
  sm = summary(mod)
  
  #%% PRINT TABLE
  # tbl <- mod %>%
  #   tbl_regression(tidy_fun = broom.mixed::tidy,
  #                  pvalue_fun = purrr::partial(style_pvalue, digits = 2),
  #                  intercept=TRUE) %>%
  #   add_global_p() %>%
  #   add_q()
  # 
  # t_out <- as_gt(tbl) %>%
  #   gt::tab_header(title=c(ei,") for Connection ",ci,"-",cj)) %>%
  #   gt::tab_options(table.layout = "fixed") %>%
  #   gt::as_raw_html()
  # print(t_out)
  
  # #-- extract fixed & random effects coeffs
  anv_vals <- car::Anova(mod,type=3)
  fix_effs <- fixef(mod)
  ran_effs <- ranef(mod)$subj_char
  
  #%% EFFECT SIZES
  cat(paste0("\n\n### model effect sizes\n"))
  #-- Intercept Model for Cohen's f^2
  print('intercept model\n');
  mod_char_int = paste(ki," ~ 1 + (1|subj_char)");
  mod <- lmerTest::lmer(str_glue(mod_char), data=tmpt);
  sfit <- summary(mod) 
  lme_est <- as.numeric(sfit$coefficients[, "Estimate"])
  lme_pval <- as.numeric(sfit$coefficients[, "Pr(>|t|)"])
  lme_stat <- as.numeric(sfit$coefficients[,"t value"])
  lme_cc <- rownames(sfit$coefficients[,])
  #-- extract fixed & random effects coeffs
  # anv_vals <- car::Anova(mod,type=3)
  # anv_chars=paste0(unlist(row.names(anv_vals)),collapse=',')
  # anv_pvals=paste0(unlist(anv_vals$`Pr(>Chisq)`),collapse=',')
  # anv_stats=paste0(unlist(anv_vals$Chisq),collapse=',')
  # anv_dfs=paste0(unlist(anv_vals$Df),collapse=',')
  #--
  anv_vals <- anova(mod)
  anv_chars=paste0(unlist(row.names(anv_vals)),collapse=',')
  anv_pvals=paste0(unlist(anv_vals$`Pr(>F)`),collapse=',')
  anv_stats=paste0(unlist(anv_vals$`F value`),collapse=',')
  anv_dfs=paste0(unlist(anv_vals$DenDF),collapse=',')
  
  fix_effs <- fixef(mod)
  ran_effs <- ranef(mod)$subj_char
  
  #%% EFFECT SIZES
  cat(paste0("\n\n### model effect sizes\n"))
  #-- Intercept Model for Cohen's f^2
  print('intercept model\n');
  mod_char_int = paste0(ei," ~ 1 + (1|subj_char)");
  tmp_mod = lme4::lmer(str_glue(mod_char_int), data=tmpt)
  cf_alt = calc_cohensf2(mod, tmp_mod)
  #-- calculate effect sizes using r-package
  anv_vals_f <- car::Anova(mod,type=3,test.statistic = "F")
  etasq_res <- effectsize::eta_squared(anv_vals_f, partial = TRUE)
  cohens_fsq_res <- effectsize::cohens_f_squared(anv_vals_f, partial = TRUE)
  #-- confidence intervals
  summ <- summary(mod);
  coeff_chars = row.names(summ$coefficients)
  emm = emmeans::emmeans(mod,spec=c('speed_n','group_char'),level=0.95);
  bfc = emmeans::contrast(emm,adjutst="bonferroni")
  pwc = emmeans::contrast(emm,method="pairwise")
  ttz = data.frame(emmeans::test(emm,null=0))
  #--
  emma = emmeans::emmeans(mod,spec=c('speed_n'),level=0.95);
  ttza = data.frame(emmeans::test(emma,null=0))
  #-- design mat
  desm = unique(mod@pp[[".->X"]]);
  desm = format_csv(data.frame(desm));
  #--
  cfis <- data.frame(emm)
  pwcs = data.frame(pwc)
  #--
  emms = emtrends(mod, pairwise ~ group_char, var = "speed_n")
  emm_slp = data.frame(emms$emtrends)
  emm_slpcon = data.frame(emms$contrasts)
  emm_slptt = test(emtrends(mod, specs = ~ group_char, var = "speed_n"), side = "two_sided")
  emm_slpcontt = data.frame(emm_slptt)
  # print(kable(cfis));
  # print(kable(anv_vals));
  # print(kable(data.frame(bfc)));
  # print(kable(data.frame(pwc)));
  
  #%% EXCEL PRINTS
  new_row = data.frame(cluster_n=0,
                       group_char=c('all'),
                       model_char=c('speed_group_intact_all'),
                       kinematic_char=ki,
                       freq_band_char="none",
                       mod_num_obs=nrow(tmpt),
                       coeff_chars=paste0(unlist(coeff_chars),collapse=','),
                       coeffs=paste0(unlist(mod@beta),collapse=','),
                       coeff_desm=desm,
                       lme_pval=paste0(unlist(lme_pval),collapse=','),
                       lme_est=paste0(unlist(lme_est),collapse=','),
                       lme_stat=paste0(unlist(lme_stat),collapse=','),
                       lme_cc=paste0(unlist(lme_cc),collapse=','),
                       pwc_cc=paste0(unlist(pwcs$contrast),collapse=','),
                       pwc_est=paste0(unlist(pwcs$estimate),collapse=','),
                       pwc_pval=paste0(unlist(pwcs$p.value),collapse=','),
                       confint_chars=paste0(unlist(cfis$group_char),collapse=','),
                       emmeans=paste0(unlist(as.numeric(cfis$emmean)),collapse=','),
                       emmeans_se=paste0(unlist(as.numeric(cfis$SE)),collapse=','),
                       confint_lwr=paste0(unlist(as.numeric(cfis$lower.CL)),collapse=','),
                       confint_upr=paste0(unlist(as.numeric(cfis$upper.CL)),collapse=','),
                       ttn0_cc=paste0(unlist(as.numeric(ttz$group_char)),collapse=','),
                       ttn0_est=paste0(unlist(as.numeric(ttz$emmean)),collapse=','),
                       ttn0_pval=paste0(unlist(as.numeric(ttz$p.value)),collapse=','),
                       ttna0_cc="all",
                       ttna0_est=paste0(unlist(as.numeric(ttza$emmean)),collapse=','),
                       ttna0_pval=paste0(unlist(as.numeric(ttza$p.value)),collapse=','),
                       emm_slp_cc=paste0(unlist(emm_slp$group_char),collapse=','),
                       emm_slp_est=paste0(unlist(as.numeric(emm_slp[,'speed_n.trend'])),collapse=','),
                       emm_slp_se=paste0(unlist(as.numeric(emm_slp$SE)),collapse=','),
                       emm_slp_df=paste0(unlist(as.numeric(emm_slp$df)),collapse=','),
                       emm_slpcon_cc=paste0(unlist(emm_slpcon$contrast),collapse=','),
                       emm_slpcon_est=paste0(unlist(as.numeric(emm_slpcon$estimate)),collapse=','),
                       emm_slpcon_pval=paste0(unlist(as.numeric(emm_slpcon$p.value)),collapse=','),
                       emm_slpcon_se=paste0(unlist(as.numeric(emm_slpcon$SE)),collapse=','),
                       emm_slpcon_trat=paste0(unlist(as.numeric(emm_slpcon$t.ratio)),collapse=','),
                       emm_slptt_cc=paste0(unlist(emm_slpcontt$group_char),collapse=','),
                       emm_slptt_est=paste0(unlist(as.numeric(emm_slpcontt[,'speed_n.trend'])),collapse=','),
                       emm_slptt_pval=paste0(unlist(as.numeric(emm_slpcontt[,'p.value'])),collapse=','),
                       emm_slptt_se=paste0(unlist(as.numeric(emm_slpcontt[,'SE'])),collapse=','),
                       emm_slptt_trat=paste0(unlist(as.numeric(emm_slpcontt[,'t.ratio'])),collapse=','),
                       anv_chars=anv_chars,
                       anv_pvals=anv_pvals,
                       anv_stats=anv_stats,
                       anv_dfs=anv_dfs,
                       r2_m_int=cf_alt$r2m,
                       r2_c_int=cf_alt$r2c,
                       f2_m_int=cf_alt$f2m,
                       f2_c_int=cf_alt$f2c,
                       fsq_chars=paste0(unlist(cohens_fsq_res$Parameter),collapse=','),
                       fsq_vals=paste0(unlist(cohens_fsq_res$Cohens_f2_partial),collapse=','),
                       etasq_chars=paste0(unlist(etasq_res$Parameter),collapse=','),
                       etasq_vals=paste0(unlist(etasq_res$Eta2_partial),collapse=','),
                       ran_effs_char=paste0(unlist(row.names(ran_effs)),collapse=','),
                       ran_effs_n=paste0(unlist(ran_effs$`(Intercept)`),collapse=','))
  excel_df <- rbind(excel_df,new_row)
  
  # #%% MODEL VALIDATION PLOTS
  # cat(paste0("\n\n### model validations\n"))
  # print(plot_model(mod, type = 'diag'))
  # cat("\n")
  #
  # tmpt$fit.m <- predict(mod, re.form = NA)
  # tmpt$fit.c <- predict(mod, re.form = NULL)
  # 
  # #%% VISUALIZATION OF MODELS
  # cat(paste0("\n\n### data plot\n"))
  # theme_set(theme_classic()) # This sets the default ggplot theme
  # print(
  #   tmpt %>%
  #     ggplot() +
  #     geom_jitter(aes(x = speed_n, y = .data[[ei]], color = subj_char), pch=20, size=10, alpha=0.2) +
  #     geom_line(aes(x= speed_n, y=fit.c, group = subj_char),linetype = "dashed", linewidth = .5) +
  #     geom_line(aes(x= speed_n, y=fit.m), linewidth = 2)+
  #     geom_errorbar(data=cfis,aes(x = speed_n, ymin=lower.CL, ymax=upper.CL), width=0.33, linewidth=1.5, colour='red') +
  #     facet_wrap(~group_char) +
  #     ggtitle(paste0('Connection ',ci,'-',cj))+
  #     xlab("Treadmill Speed (m/s)") +
  #     ylab(eeg_title_chars[eis[i]]) +
  #     guides(group="none")
  # )
}


#%% HTML TAG
cat(paste0("\n\n# Kin Group Models\n"))
# SPEED GROUP
theme_set(theme_classic()) # This sets the default ggplot theme
#%% LOOP VARS
kis <- rep(NA,length(kin_measures))
cnt = 1;
for (i in 1:length(kin_measures)) {
  kis[cnt] = i;
  cnt = cnt + 1;
}
for (i in 1:length(kis)){
  #%% SET LOOP PARAMS
  ki = kin_measures[kis[i]];

  #%% HTML TAG
  cat(paste0("\n\n## ", ki,"\n"))
  
  #%% CREATE MEAN & STD MEASURES PER CONDITION
  tmpt = dtbl;
  tmpt <- tmpt %>%
    group_by(subj_char,cond_char,speed_n,model_n,group_char,group_n) %>%
    summarise_all(list("mean"))
  
  #%% SUBSET BY CLUSTER
  #-- remove NAN cases
  tmpt = tmpt[complete.cases(tmpt[[ki]]),];
  
  #%% COMPUTE LME MODEL
  mod_char = paste(ki," ~ 1 + speed_n + group_char + (1|subj_char)");
  # mod <- lme4::lmer(str_glue(mod_char), data=tmpt);
  mod <- lmerTest::lmer(str_glue(mod_char), data=tmpt);
  sfit <- summary(mod)
  
  #%% PRINT TABLE
  # tbl <- mod %>%
  #   tbl_regression(tidy_fun = broom.mixed::tidy,
  #                  pvalue_fun = purrr::partial(style_pvalue, digits = 2),
  #                  intercept=TRUE) %>%
  #   add_global_p() %>%
  #   add_q()
  # 
  # t_out <- as_gt(tbl) %>%
  #   gt::tab_header(title=c(ei,") for Connection ",ci,"-",cj)) %>%
  #   gt::tab_options(table.layout = "fixed") %>%
  #   gt::as_raw_html()
  # print(t_out)
  
  lme_est <- as.numeric(sfit$coefficients[, "Estimate"])
  lme_pval <- as.numeric(sfit$coefficients[, "Pr(>|t|)"])
  lme_stat <- as.numeric(sfit$coefficients[,"t value"])
  lme_cc <- rownames(sfit$coefficients[,])
  #-- extract fixed & random effects coeffs
  # anv_vals <- car::Anova(mod,type=3)
  # anv_chars=paste0(unlist(row.names(anv_vals)),collapse=',')
  # anv_pvals=paste0(unlist(anv_vals$`Pr(>Chisq)`),collapse=',')
  # anv_stats=paste0(unlist(anv_vals$Chisq),collapse=',')
  # anv_dfs=paste0(unlist(anv_vals$Df),collapse=',')
  #--
  anv_vals <- anova(mod)
  anv_chars=paste0(unlist(row.names(anv_vals)),collapse=',')
  anv_pvals=paste0(unlist(anv_vals$`Pr(>F)`),collapse=',')
  anv_stats=paste0(unlist(anv_vals$`F value`),collapse=',')
  anv_dfs=paste0(unlist(anv_vals$DenDF),collapse=',')
  
  fix_effs <- fixef(mod)
  ran_effs <- ranef(mod)$subj_char
  
  #%% EFFECT SIZES
  cat(paste0("\n\n### model effect sizes\n"))
  #-- Intercept Model for Cohen's f^2
  print('intercept model\n');
  mod_char_int = paste0(ki," ~ 1 + (1|subj_char)");
  tmp_mod = lme4::lmer(str_glue(mod_char_int), data=tmpt)
  cf_alt = calc_cohensf2(mod, tmp_mod)
  #-- calculate effect sizes using r-package
  anv_vals_f <- car::Anova(mod,type=3,test.statistic = "F")
  etasq_res <- effectsize::eta_squared(anv_vals_f, partial = TRUE)
  cohens_fsq_res <- effectsize::cohens_f_squared(anv_vals_f, partial = TRUE)
  #-- confidence intervals
  summ <- summary(mod);
  coeff_chars = row.names(summ$coefficients)
  emm = emmeans::emmeans(mod,spec=c('speed_n','group_char'),level=0.95);
  bfc = emmeans::contrast(emm,adjutst="bonferroni")
  pwc = emmeans::contrast(emm,method="pairwise")
  ttz = data.frame(emmeans::test(emm,null=0))
  #--
  emma = emmeans::emmeans(mod,spec=c('speed_n'),level=0.95);
  ttza = data.frame(emmeans::test(emma,null=0))
  #-- design mat
  desm = unique(mod@pp[[".->X"]]);
  desm = format_csv(data.frame(desm));
  #--
  cfis <- data.frame(emm)
  pwcs = data.frame(pwc)
  #--
  # print(kable(cfis));
  # print(kable(anv_vals));
  # print(kable(data.frame(bfc)));
  # print(kable(data.frame(pwc)));
  
  #%% EXCEL PRINTS
  new_row = data.frame(cluster_n=0,
                       group_char=c('all'),
                       model_char=c('speed_group_all'),
                       kinematic_char=ki,
                       freq_band_char="none",
                       mod_num_obs=nrow(tmpt),
                       coeff_chars=paste0(unlist(coeff_chars),collapse=','),
                       coeffs=paste0(unlist(mod@beta),collapse=','),
                       coeff_desm=desm,
                       lme_pval=paste0(unlist(lme_pval),collapse=','),
                       lme_est=paste0(unlist(lme_est),collapse=','),
                       lme_stat=paste0(unlist(lme_stat),collapse=','),
                       lme_cc=paste0(unlist(lme_cc),collapse=','),
                       pwc_cc=paste0(unlist(pwcs$contrast),collapse=','),
                       pwc_est=paste0(unlist(pwcs$estimate),collapse=','),
                       pwc_pval=paste0(unlist(pwcs$p.value),collapse=','),
                       confint_chars=paste0(unlist(cfis$group_char),collapse=','),
                       emmeans=paste0(unlist(as.numeric(cfis$emmean)),collapse=','),
                       emmeans_se=paste0(unlist(as.numeric(cfis$SE)),collapse=','),
                       confint_lwr=paste0(unlist(as.numeric(cfis$lower.CL)),collapse=','),
                       confint_upr=paste0(unlist(as.numeric(cfis$upper.CL)),collapse=','),
                       ttn0_cc=paste0(unlist(as.numeric(ttz$group_char)),collapse=','),
                       ttn0_est=paste0(unlist(as.numeric(ttz$emmean)),collapse=','),
                       ttn0_pval=paste0(unlist(as.numeric(ttz$p.value)),collapse=','),
                       ttna0_cc="all",
                       ttna0_est=paste0(unlist(as.numeric(ttza$emmean)),collapse=','),
                       ttna0_pval=paste0(unlist(as.numeric(ttza$p.value)),collapse=','),
                       emm_slp_cc="none",
                       emm_slp_est="none",
                       emm_slp_se="none",
                       emm_slp_df="none",
                       emm_slpcon_cc="none",
                       emm_slpcon_est="none",
                       emm_slpcon_pval="none",
                       emm_slpcon_se="none",
                       emm_slpcon_trat="none",
                       emm_slptt_cc="none",
                       emm_slptt_est="none",
                       emm_slptt_pval="none",
                       emm_slptt_se="none",
                       emm_slptt_trat="none",
                       anv_chars=anv_chars,
                       anv_pvals=anv_pvals,
                       anv_stats=anv_stats,
                       anv_dfs=anv_dfs,
                       r2_m_int=cf_alt$r2m,
                       r2_c_int=cf_alt$r2c,
                       f2_m_int=cf_alt$f2m,
                       f2_c_int=cf_alt$f2c,
                       fsq_chars=paste0(unlist(cohens_fsq_res$Parameter),collapse=','),
                       fsq_vals=paste0(unlist(cohens_fsq_res$Cohens_f2_partial),collapse=','),
                       etasq_chars=paste0(unlist(etasq_res$Parameter),collapse=','),
                       etasq_vals=paste0(unlist(etasq_res$Eta2_partial),collapse=','),
                       ran_effs_char=paste0(unlist(row.names(ran_effs)),collapse=','),
                       ran_effs_n=paste0(unlist(ran_effs$`(Intercept)`),collapse=','))
  excel_df <- rbind(excel_df,new_row)
  
  # #%% MODEL VALIDATION PLOTS
  # cat(paste0("\n\n### model validations\n"))
  # print(plot_model(mod, type = 'diag'))
  # cat("\n")
  # 
  # tmpt$fit.m <- predict(mod, re.form = NA)
  # tmpt$fit.c <- predict(mod, re.form = NULL)
  # 
  # #%% VISUALIZATION OF MODELS
  # cat(paste0("\n\n### data plot\n"))
  # theme_set(theme_classic()) # This sets the default ggplot theme
  # print(
  #   tmpt %>%
  #     ggplot() +
  #     geom_jitter(aes(x = speed_n, y = .data[[ei]], color = subj_char), pch=20, size=10, alpha=0.2) +
  #     geom_line(aes(x= speed_n, y=fit.c, group = subj_char),linetype = "dashed", linewidth = .5) +
  #     geom_line(aes(x= speed_n, y=fit.m), linewidth = 2)+
  #     geom_errorbar(data=cfis,aes(x = speed_n, ymin=lower.CL, ymax=upper.CL), width=0.33, linewidth=1.5, colour='red') +
  #     facet_wrap(~group_char) +
  #     ggtitle(paste0('Connection from',ci,'->',cj))+
  #     xlab("Treadmill Speed (m/s)") +
  #     ylab(eeg_title_chars[eis[i]]) +
  #     guides(group="none")
  # )
}

#%% ========================================================================= %%#
write.xlsx(excel_df,file.path(save_dir,EXCEL_FNAME))
