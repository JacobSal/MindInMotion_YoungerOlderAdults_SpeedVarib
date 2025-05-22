# install.packages(c("parallel","purrr","dplyr","tibble","Matrix",
#                    "lattice","gridExtra","R.devices","R.matlab","ggplot2",
#                    "gridExtra","genlasso"));
print('Clearing Workspace');
# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console
cat("\014") 
# Clean workspace
rm(list=ls())

print('Loading packages...')
#%% UTILITY
library(parallel);
library(purrr);
library(dplyr);
library(tibble);
library(Matrix)
library(lattice);
library(R.devices);
library(ggplot2)
library(gridExtra);
library(grid)
library(lmerTest)
#%% PACKAGES FOR STATS
library(genlasso);
library(R.matlab);

#%% CUSTOM FUNCTIONS
# curr_dir = "/blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/r_scripts/sbs_lme_mods";
curr_dir = getwd();
setwd(curr_dir)
source(file.path(curr_dir,"tf_flasso_funcs.R"));
print(curr_dir);

#%% LOAD DATA ============================================================== %%#
print("Loading data...");
# clusters = c(3,4,5,6,7,8,9,10,11,12,13) # RSup/RSM, PreC, LSM, Mid Cing, LSup, LPPA, RPPA
clusters = c(3,4,6,8)
conds = c(1,2,3,4)
fext = 'itc_rdata_table_phasec_notw_mw_based'

#%% CREATE SAVE DIR
curr_dir <- getwd();
# fname = "itc_rdata_table_phasec_notw_mw_based_flasso_based_results";
# fname = "itc_rdata_table_phasec_notw_mw_based_fl_res_bsz5_nob";
fname = "itc_rdata_flasso_125f_out_bsz5_nob_sliding";
save_dir <- file.path(curr_dir,fname)
dir.create(save_dir);

#%% LOAD
# fname = "itc_rdata_struct_phasec_notw_mw_based.mat";
fname = "itc_rdata_struct_125f_phasec_notw_mw_based.mat";
# fname = "itc_rdata_cell_phasec_notw_mw_based.mat";
dpath = "/jsalminen/GitHub/MIND_IN_MOTION_PRJ/_data/MIM_dataset/_studies/02202025_mim_yaoa_powpow0p3_crit_speed/__iclabel_cluster_allcond_rb3/icrej_5/11/kin_eeg_step_to_step"
mat_fpath <- file.path(dpath,fname)
#--
if(ispc()){
  mat_fpath <- paste0("M:",mat_fpath)
}else{
  mat_fpath <- paste0("/blue/dferris",mat_fpath);
}
#-- load mat data and manipulate
tmp_mat = R.matlab::readMat(mat_fpath)
#-- extract data
dato = get_mat_dat(tmp_mat);
indexl <- dato$indexl;
itc_datl <- dato$itc_datl;
freqs <- dato$freqs;
times <- dato$times;
nfreqs = length(freqs)
ntimes = length(times);
#-- filter
indexl <- filter_at(indexl,vars('cond_n'), any_vars(. %in% conds));
indexl <- filter_at(indexl,vars('cluster_n'), any_vars(. %in% clusters));

#%% LOOP VARS ============================================================== %%#
clusters=unique(indexl$cluster_n);
conds=unique(indexl$cond_n);
subjs=unique(indexl$subj_n);
datl <- ntimes*nfreqs;

#%% ======================================================================== %%#
X1 <- array(rnorm(107*60*4), dim=c(107,60,4))
X2 <- array(rnorm(107*60*4), dim=c(107,60,3))

speed = c(0.25,0.5,0.75,1,0.25,0.5,1)

subid = c(1,1,1,1,2,2,2)
library(lmerTest)

estimate <- matrix(0, 107, 60)
pvalue   <- matrix(0, 107, 60)

for(i in 1:107){
  for(j in 1:60){
    y <- c(X1[i,j,], X2[i,j,])
    
    fit <- lmer(y~speed + (1|subid) )
    
    summary_fit <- summary(fit)
    estimate[i,j] <- summary_fit$coefficients["speed", "Estimate"]
    pvalue[i,j] <- summary_fit$coefficients["speed", "Pr(>|t|)"]
  }
}