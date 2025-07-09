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
library(Matrix);
library(lattice);
library(R.devices);
library(ggplot2);
library(gridExtra);
library(grid);
library(lmerTest);
library(emmeans);
# library(lme4)
#%% PACKAGES FOR STATS
library(genlasso);
# library(car);
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
clusters = c(3,4,6,8,5,9,11)
conds = c(1,2,3,4)

#%% CREATE SAVE DIR
curr_dir <- getwd();

#%% LOAD
# fname = "itc_rdata_struct_phasec_notw_mw_based.mat";
# fname = "itc_rdata_struct_125f_phasec_notw_mw_based";
# fname = "itc_rdata_cell_phasec_notw_mw_based.mat";
# fname = "rdata_struct_ersp_spca";
# fname = "itc_rdata_struct_ext_crop_phasec_notw_based"
fname = "rdata_struct_ext_crop_phasec_notw_based"
fname_save = "rdata_extc_phasec_notw_condb"
dpath = "/jsalminen/GitHub/MIND_IN_MOTION_PRJ/_data/MIM_dataset/_studies/02202025_mim_yaoa_powpow0p3_crit_speed/__iclabel_cluster_allcond_rb3/icrej_5/11/kin_eeg_step_to_step"
#--
mat_fpath <- file.path(dpath,sprintf('%s.mat',fname))
save_dir <- file.path(curr_dir,fname_save)
dir.create(save_dir);
#--
if(ispc()){
  mat_fpath <- paste0("M:",mat_fpath)
}else{
  mat_fpath <- paste0("/blue/dferris",mat_fpath);
}
#-- load mat data and manipulate
tmp_mat = R.matlab::readMat(mat_fpath)
#-- extract data
# FREQ_BOUND = c(3,60);
# TIME_BOUND = c(0,1450);
FREQ_BOUND = c(3,80);
TIME_BOUND = c(0,1500);
# FREQ_BOUND = c(3,13); # test
# TIME_BOUND = c(0,200); # test
#--
dato = get_mat_dat(tmp_mat);
indexl <- dato$indexl;
itc_datl <- dato$itc_datl;
freqs <- dato$freqs;
times <- dato$times;
# finds = !vector(mode="logical",length=length(freqs));
# tinds = !vector(mode="logical",length=length(times));
nfreqs = length(freqs)
ntimes = length(times);
#--
finds = freqs > FREQ_BOUND[1] & freqs < FREQ_BOUND[2];
tinds = times > TIME_BOUND[1] & times < TIME_BOUND[2];
# freqs = freqso[finds];
# times = timeso[tinds];
# nfreqs = length(freqs);
# ntimes = length(times);
#-- filter
indexl <- filter_at(indexl,vars('cond_n'), any_vars(. %in% conds));
indexl <- filter_at(indexl,vars('cluster_n'), any_vars(. %in% clusters));
#--
clusters=unique(indexl$cluster_n);
conds=unique(indexl$cond_n);
subjs=unique(indexl$subj_n);
datl <- ntimes*nfreqs;
cond_vals = c(0.25,0.50,0.75,1.0);

#%% LOOP =================================================================== %%#
#-- RUN MPI
print("Running MPI");
if(ispc()){
  numCores <- detectCores()
}else{
  numCores = as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
}
print(numCores)
clust <- makeCluster(numCores)
set.seed(42069)

for(i in 1:length(clusters)){
  #-- get dirs and clust
  cli = clusters[i];
  tmp_save_dir <- file.path(save_dir,sprintf("cl%i_plots",cli))
  dir.create(tmp_save_dir);
  #--
  vals <- get_stat_vecs(indexl,cli,freqs,times,finds,tinds,
                                    save_dir)
  #%%
  loop_dat <- get_stat_dat_mat(vals$dat_mat,vals$nfreqs,vals$ntimes,
                               vals$speed_vals,vals$group_vals,vals$subj_vals);
  #--
  system.time(saves_grp <- mclapply(loop_dat,function(x) lmer_fl_intstat_ij(x,DO_INT_MOD=FALSE),
                                mc.preschedule = FALSE,
                                mc.set.seed = FALSE,
                                mc.cores=numCores));
  #--
  system.time(saves_int <- mclapply(loop_dat,function(x) lmer_fl_intstat_ij(x,DO_INT_MOD=TRUE),
                                mc.preschedule = FALSE,
                                mc.set.seed = FALSE,
                                mc.cores=numCores));
  
  #%% RUN LOCAL
  # loop_dat <- get_stat_dat_mat(vals$dat_mat,vals$nfreqs,vals$ntimes,
  #                              vals$speed_vals,vals$group_vals,vals$subj_vals);
  #-- unit test
  # saves_int <- lmer_fl_intstat_ij(loop_dat[[1]]);
  #--
  # saves_grp <- lapply(loop_dat[1:10],function(x) lmer_fl_intstat_ij(x,DO_INT_MOD=FALSE))
  # dato <- flstat_int_agg(saves_grp,cli,freqs,times,pinds=1:4,apinds=1:2,pwcinds=1:3)
  # saves_int <- lapply(loop_dat[1:10],function(x) lmer_fl_intstat_ij(x,DO_INT_MOD=TRUE))
  # dato <- flstat_int_agg(saves_int,cli,freqs,times,pinds=1:6,apinds=1:3,pwcinds=1:3)

  #%% EXTRACT DATA
   
  esto <- array(dato$estimate,dim=c(vals$nfreqs,vals$ntimes,length(dato$pinds)));
  fdrp <- array(dato$fdrp,dim=c(vals$nfreqs,vals$ntimes,length(dato$pinds)));
  astato <- array(dato$astat,dim=c(vals$nfreqs,vals$ntimes,length(dato$apinds)));
  afdrp <- array(dato$afdrp,dim=c(vals$nfreqs,vals$ntimes,length(dato$apinds)));
  pwcfdrp <- array(dato$pwcfdrp,dim=c(vals$nfreqs,vals$ntimes,length(dato$pwcinds)));
  
  #%% SAVE
  fname = sprintf("statmatint_cl%i.mat",cli);
  writeMat(con=file.path(save_dir,fname),stat_mat=dato)
  fname = sprintf("statrdsint_cl%i.RData",cli);
  saveRDS(dato, file=file.path(save_dir,fname))
  
  #%% PLOT
  zlim_in = range(esto,na.rm=TRUE)
  p1 <- tf_plot(esto[,,1],vals$times,vals$freqs,"intercept",zlim_in,
                do_cbar=FALSE)
  p2 <- tf_plot(esto[,,2],vals$times,vals$freqs,"speed",zlim_in,
                do_cbar=FALSE)
  p3 <- tf_plot(esto[,,3],vals$times,vals$freqs,"grp2-grp1",zlim_in,
                do_cbar=FALSE)
  p4 <- tf_plot(esto[,,4],vals$times,vals$freqs,"grp3-grp1",zlim_in,
                do_cbar=TRUE)
  p5 <- tf_plot(esto[,,5],vals$times,vals$freqs,"sp:grp2",zlim_in,
                do_cbar=FALSE)
  p6 <- tf_plot(esto[,,6],vals$times,vals$freqs,"sp:grp3",zlim_in,
                do_cbar=TRUE)
  zlim_in = range(afdrp,na.rm=TRUE)
  p11 <- tf_plot(astato[,,1],vals$times,vals$freqs,"aSpeed",zlim_in,
                 do_cbar=FALSE)
  p22 <- tf_plot(astato[,,2],vals$times,vals$freqs,"aGroup",zlim_in,
                 do_cbar=FALSE)
  p33 <- tf_plot(astato[,,3],vals$times,vals$freqs,"aInt",zlim_in,
                 do_cbar=TRUE)
  #-- join them
  pl <- list(p1,p2,p3,p4,p5,p6,p11,p22,p33)
  grid.arrange(grid::rectGrob(),grid::rectGrob())
  ml <- marrangeGrob(pl, nrow=1, ncol=2);
  #-- resave tf-plots?
  fname = sprintf("allestint_cl%i_tfplot.pdf",cli);
  ggsave(file.path(tmp_save_dir,fname), ml)
  
  #%%
  zlim_in = range(fdrp,na.rm=TRUE)
  p1 <- tf_plot(fdrp[,,1],vals$times,vals$freqs,"intercept",zlim_in,
                do_cbar=FALSE)
  p2 <- tf_plot(fdrp[,,2],vals$times,vals$freqs,"speed",zlim_in,
                do_cbar=FALSE)
  p3 <- tf_plot(fdrp[,,3],vals$times,vals$freqs,"grp2-grp1",zlim_in,
                do_cbar=FALSE)
  p4 <- tf_plot(fdrp[,,4],vals$times,vals$freqs,"grp3-grp1",zlim_in,
                do_cbar=TRUE)
  p5 <- tf_plot(fdrp[,,5],vals$times,vals$freqs,"sp:grp2",zlim_in,
                do_cbar=FALSE)
  p6 <- tf_plot(fdrp[,,6],vals$times,vals$freqs,"sp:grp3",zlim_in,
                do_cbar=TRUE)
  p11 <- tf_plot(afdrp[,,1],vals$times,vals$freqs,"aSpeed",zlim_in,
                 do_cbar=FALSE)
  p22 <- tf_plot(afdrp[,,2],vals$times,vals$freqs,"aGroup",zlim_in,
                 do_cbar=FALSE)
  p33 <- tf_plot(afdrp[,,3],vals$times,vals$freqs,"aInt",zlim_in,
                 do_cbar=TRUE)
  p111 <- tf_plot(pwcfdrp[,,1],vals$times,vals$freqs,dato$pwc_c[1],zlim_in,
                  do_cbar=FALSE)
  p222 <- tf_plot(pwcfdrp[,,2],vals$times,vals$freqs,dato$pwc_c[2],zlim_in,
                  do_cbar=FALSE)
  p333 <- tf_plot(pwcfdrp[,,3],vals$times,vals$freqs,dato$pwc_c[3],zlim_in,
                  do_cbar=TRUE)
  #-- join them
  pl <- list(p1,p2,p3,p4,p5,p6,p11,p22,p33,p111,p222,p333)
  grid.arrange(grid::rectGrob(),grid::rectGrob())
  ml <- marrangeGrob(pl, nrow=1, ncol=2);
  #-- resave tf-plots?
  fname = sprintf("allpvalint_cl%i_tfplot.pdf",cli);
  ggsave(file.path(tmp_save_dir,fname), ml)
  
  
  #%% EXTRACT DATA
  dato <- flstat_int_agg(saves_grp,cli,vals$freqs,vals$times,pinds=1:4,apinds=1:2,pwcinds=1:3)
  esto <- array(dato$estimate,dim=c(vals$nfreqs,vals$ntimes,length(dato$pinds)));
  fdrp <- array(dato$fdrp,dim=c(vals$nfreqs,vals$ntimes,length(dato$pinds)));
  afdrp <- array(dato$afdrp,dim=c(vals$nfreqs,vals$ntimes,length(dato$apinds)));
  astato <- array(dato$astat,dim=c(vals$nfreqs,vals$ntimes,length(dato$apinds)));
  pwcfdrp <- array(dato$pwcfdrp,dim=c(vals$nfreqs,vals$ntimes,length(dato$pwcinds)));
  
  
  #%% SAVE
  fname = sprintf("statmatgrp_cl%i.mat",cli);
  writeMat(con=file.path(save_dir,fname),stat_mat=dato)
  fname = sprintf("statrdsgrp_cl%i.RData",cli);
  saveRDS(dato, file=file.path(save_dir,fname))
  
  #%% PLOT
  zlim_in = range(esto,na.rm=TRUE)
  p1 <- tf_plot(esto[,,1],vals$times,vals$freqs,"intercept",zlim_in,
                do_cbar=FALSE)
  p2 <- tf_plot(esto[,,2],vals$times,vals$freqs,"speed",zlim_in,
                do_cbar=FALSE)
  p3 <- tf_plot(esto[,,3],vals$times,vals$freqs,"grp2-grp1",zlim_in,
                do_cbar=FALSE)
  p4 <- tf_plot(esto[,,4],vals$times,vals$freqs,"grp3-grp1",zlim_in,
                do_cbar=TRUE)
  zlim_in = range(afdrp,na.rm=TRUE)
  p11 <- tf_plot(afdrp[,,1],vals$times,vals$freqs,"aSpeed",zlim_in,
                 do_cbar=FALSE)
  p22 <- tf_plot(afdrp[,,2],vals$times,vals$freqs,"aGroup",zlim_in,
                 do_cbar=TRUE)
  #-- join them
  pl <- list(p1,p2,p3,p4,p11,p22)
  grid.arrange(grid::rectGrob(),grid::rectGrob())
  ml <- marrangeGrob(pl, nrow=1, ncol=2);
  #-- resave tf-plots?
  fname = sprintf("allestgrp_cl%i_tfplot.pdf",cli);
  ggsave(file.path(tmp_save_dir,fname), ml)
  
  #%%
  zlim_in = range(fdrp,na.rm=TRUE)
  p1 <- tf_plot(fdrp[,,1],vals$times,vals$freqs,"intercept",zlim_in,
                do_cbar=FALSE)
  p2 <- tf_plot(fdrp[,,2],vals$times,vals$freqs,"speed",zlim_in,
                do_cbar=FALSE)
  p3 <- tf_plot(fdrp[,,3],vals$times,vals$freqs,"grp2-grp1",zlim_in,
                do_cbar=FALSE)
  p4 <- tf_plot(fdrp[,,4],vals$times,vals$freqs,"grp3-grp1",zlim_in,
                do_cbar=TRUE)
  p11 <- tf_plot(afdrp[,,1],vals$times,vals$freqs,"aSpeed",zlim_in,
                 do_cbar=FALSE)
  p22 <- tf_plot(afdrp[,,2],vals$times,vals$freqs,"aGroup",zlim_in,
                 do_cbar=TRUE)
  p111 <- tf_plot(pwcfdrp[,,1],vals$times,vals$freqs,dato$pwc_c[1],zlim_in,
                  do_cbar=FALSE)
  p222 <- tf_plot(pwcfdrp[,,2],vals$times,vals$freqs,dato$pwc_c[2],zlim_in,
                  do_cbar=FALSE)
  p333 <- tf_plot(pwcfdrp[,,3],vals$times,vals$freqs,dato$pwc_c[3],zlim_in,
                  do_cbar=TRUE)
  #-- join them
  pl <- list(p1,p2,p3,p4,p11,p22,p111,p222,p333)
  grid.arrange(grid::rectGrob(),grid::rectGrob())
  ml <- marrangeGrob(pl, nrow=1, ncol=2);
  #-- resave tf-plots?
  fname = sprintf("allpvalgrp_cl%i_tfplot.pdf",cli);
  ggsave(file.path(tmp_save_dir,fname), ml)


}