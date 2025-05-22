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
# library(lmerTest)
library(lme4)
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
fname = "itc_rdata_flasso_125f_out_bsz5_nob_sliding";
save_dir <- file.path(curr_dir,fname)
dir.create(save_dir);

#%% LOAD
fname = "itc_rdata_struct_125f_phasec_notw_mw_based.mat";
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
FREQ_BOUND = c(3,60);
TIME_BOUND = c(0,1450);
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

#%% LOOP VARS ============================================================== %%#
clusters=unique(indexl$cluster_n);
conds=unique(indexl$cond_n);
subjs=unique(indexl$subj_n);
datl <- ntimes*nfreqs;
cond_vals = c(0.25,0.50,0.75,1.0);

for(i in 1:length(clusters)){
  #-- get dirs and clust
  cli = clusters[i];
  tmp_save_dir <- file.path(save_dir,sprintf("cl%i_plots",cli))
  dir.create(tmp_save_dir);
  #--
  vals <- get_stat_vecs(indexl,cli,freqs,times,finds,tinds,
                                    save_dir)
  
  #%% RUN MPI
  # print("Running MPI");
  # if(ispc()){
  #   numCores <- detectCores()
  # }else{
  #   numCores = as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
  # }
  # print(numCores)
  # clust <- makeCluster(numCores)
  # set.seed(42069)
  #%%
  loop_dat <- get_stat_dat_mat(vals$dat_mat,vals$nfreqs,vals$ntimes,
                               vals$speed_vals,vals$group_vals,vals$subj_vals);
  #--
  system.time(saves <- mclapply(loop_dat,function(x) lmer_fl_grpstat_ij(x),
                                mc.preschedule = FALSE,
                                mc.set.seed = FALSE,
                                mc.cores=numCores));
  #--
  system.time(saves <- mclapply(loop_dat,function(x) lmer_fl_intstat_ij(x),
                                mc.preschedule = FALSE,
                                mc.set.seed = FALSE,
                                mc.cores=numCores));
  
  #%% RUN LOCAL
  # loop_dat <- get_stat_dat_mat(vals$dat_mat,vals$nfreqs,vals$ntimes,
  #                              vals$speed_vals,vals$group_vals,vals$subj_vals);
  # saves <- lapply(loop_dat,function(x) lmer_fl_grpstat_ij(x))
  # saves <- lapply(loop_dat,function(x) lmer_fl_intstat_ij(x))

  #--
  tmp <- lapply(saves,function(x) x$i)
  is <- as.numeric(tmp)
  tmp <- lapply(saves,function(x) x$j)
  js <- as.numeric(tmp)
  # tmp <- lapply(saves,function(x) x$cnt)
  # cnts <- as.numeric(tmp)
  tmp <- lapply(saves,function(x) x$statnames)
  statnames <- unique(as.character(tmp))
  # unm <- sort(cnts,index.return=TRUE)
  #--
  pvo <- array(0,dim=c(vals$nfreqs,vals$ntimes,3))
  esto <- array(0,dim=c(vals$nfreqs,vals$ntimes,3))
  cnt = 1;
  for(cnt in 1:length(saves)){
    i = as.numeric(saves[[cnt]]$i);
    j = as.numeric(saves[[cnt]]$j);
    pv = as.numeric(saves[[cnt]]$pv);
    est = as.numeric(saves[[cnt]]$est);
    pvo[i,j,1] = pv[1];
    pvo[i,j,2] = pv[2];
    pvo[i,j,3] = pv[3];
    pvo[i,j,4]
    pvo[i,j,5]
    esto[i,j,1] = est[1];
    esto[i,j,2] = est[2];
    esto[i,j,3] = est[3];
    
  }

  #%% FDR CORRECTION
  fdrp <- p.adjust(matrix(pvo,nrow=vals$nfreqs*vals$ntimes*3),
                   method="fdr",
                   n=vals$nfreqs*vals$ntimes*3);
  fdrp <- array(fdrp,c(vals$nfreqs,vals$ntimes,3));

  #%% PLOT
  zlim_in = range(esto)
  p1 <- tf_plot(esto[,,1],times[tinds],freqs[finds],"speed",zlim_in,
                do_cbar=FALSE)
  p2 <- tf_plot(esto[,,2],times[tinds],freqs[finds],"grp2-grp1",zlim_in,
                do_cbar=FALSE)
  p3 <- tf_plot(esto[,,3],times[tinds],freqs[finds],"grp3-grp1",zlim_in,
                do_cbar=TRUE)
  #-- join them
  pl <- list(p1,p2,p3)
  grid.arrange(grid::rectGrob(),grid::rectGrob())
  ml <- marrangeGrob(pl, nrow=1, ncol=2);
  #-- resave tf-plots?
  fname = sprintf("allest_cl%i_tfplot.png",cli);
  ggsave(file.path(tmp_save_dir,fname), ml)

  #%%
  zlim_in = range(fdrp)
  p1 <- tf_plot(fdrp[,,1],times[tinds],freqs[finds],"speed",zlim_in,
                do_cbar=FALSE)
  p2 <- tf_plot(fdrp[,,2],times[tinds],freqs[finds],"grp2-grp1",zlim_in,
                do_cbar=FALSE)
  p3 <- tf_plot(fdrp[,,3],times[tinds],freqs[finds],"grp3-grp1",zlim_in,
                do_cbar=TRUE)
  #-- join them
  pl <- list(p1,p2,p3)
  grid.arrange(grid::rectGrob(),grid::rectGrob())
  ml <- marrangeGrob(pl, nrow=1, ncol=2);
  #-- resave tf-plots?
  fname = sprintf("allpval_cl%i_tfplot.png",cli);
  ggsave(file.path(tmp_save_dir,fname), ml)

  #%% SAVE
  dato = list(fdrp=matrix(fdrp,nrow=vals$nfreqs*vals$ntimes*3),
              estimate=matrix(esto,nrow=vals$nfreqs*vals$ntimes*3),
              freqs=freqs[finds],
              times=times[tinds],
              dim3=statnames);
  
  fname = sprintf("statmat_cl%i.mat",cli);
  writeMat(con=file.path(save_dir,fname),stat_mat=dato)
  fname = sprintf("statrds_cl%i.RData",cli);
  saveRDS(dato, file=file.path(save_dir,fname))
  
  #%% SINGLE SPEED TEST
  # speed_i = 1.0;
  # cond_i = 4;
  # loop_dat <- get_statdat_onesubj(vals$dat_mat,vals$nfreqs,vals$ntimes,speed_i,cond_i,
  #                                 vals$speed_vals,vals$group_vals,vals$subj_vals)
  # saves <- lapply(loop_dat,function(x) flstat_onecond_ij(x))
  # system.time(saves <- mclapply(loop_dat,function(x) flstat_onecond_ij(x),
  #                               mc.preschedule = FALSE,
  #                               mc.set.seed = FALSE,
  #                               mc.cores=numCores));
  # 
  # #--
  # tmp <- lapply(saves,function(x) x$i)
  # is <- as.numeric(tmp)
  # tmp <- lapply(saves,function(x) x$j)
  # js <- as.numeric(tmp)
  # tmp <- lapply(saves,function(x) x$cnt)
  # cnts <- as.numeric(tmp)
  # unm <- sort(cnts,index.return=TRUE)
  # #--
  # pvo <- array(0,dim=c(vals$nfreqs,vals$ntimes,2))
  # esto <- array(0,dim=c(vals$nfreqs,vals$ntimes,2))
  # cnt = 1;
  # for(cnt in 1:length(saves)){
  #   i = as.numeric(saves[[cnt]]$i);
  #   j = as.numeric(saves[[cnt]]$j);
  #   pv = as.numeric(saves[[cnt]]$pv);
  #   est = as.numeric(saves[[cnt]]$est);
  #   pvo[i,j,1] = pv[1];
  #   pvo[i,j,2] = pv[2];
  #   esto[i,j,1] = est[1];
  #   esto[i,j,2] = est[2]
  # }
  # 
  # #%% FDR CORRECTION
  # fdrp <- p.adjust(matrix(pvo,nrow=vals$nfreqs*vals$ntimes*2),
  #                  method="fdr",
  #                  n=vals$nfreqs*vals$ntimes*2);
  # fdrp <- array(fdrp,c(vals$nfreqs,vals$ntimes,2));
  # # fdrp <- array(fdrp,c(vals$ntimes,vals$nfreqs,2));
  # 
  # #%% PLOT
  # zlim_in = range(esto)
  # p2 <- tf_plot(esto[,,1],times,freqs,"grp1-grp3",zlim_in,
  #               do_cbar=FALSE)
  # p3 <- tf_plot(esto[,,2],times,freqs,"grp2-grp3",zlim_in,
  #               do_cbar=TRUE)
  # #-- join them
  # pl <- list(p2,p3)
  # grid.arrange(grid::rectGrob(),grid::rectGrob())
  # ml <- marrangeGrob(pl, nrow=1, ncol=2);
  # 
  # #%%
  # zlim_in = range(fdrp)
  # p2 <- tf_plot(fdrp[,,1],times,freqs,"grp1-grp3",zlim_in,
  #               do_cbar=FALSE)
  # p3 <- tf_plot(fdrp[,,2],times,freqs,"grp2-grp3",zlim_in,
  #               do_cbar=TRUE)
  # #-- join them
  # pl <- list(p2,p3)
  # grid.arrange(grid::rectGrob(),grid::rectGrob())
  # ml <- marrangeGrob(pl, nrow=1, ncol=2);
}