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
curr_dir = "/blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/r_scripts/sbs_lme_mods";
# curr_dir = getwd();
setwd(curr_dir)
source(file.path(curr_dir,"tf_flasso_funcs.R"));
print(curr_dir);

#%% LOAD DATA ============================================================== %%#
print("Loading data...");
# clusters = c(3,4,5,6,7,8,9,10,11,12,13) # RSup/RSM, PreC, LSM, Mid Cing, LSup, LPPA, RPPA
clusters = c(3,4,6,8,5,9,11)
conds = c(1,2,3,4)
# conds = c(5,6,7,8)

#%% CREATE SAVE DIR
curr_dir <- getwd();

#%% LOAD
# fname = "itc_rdata_struct_phasec_notw_mw_based.mat";
# fname = "itc_rdata_struct_125f_phasec_notw_mw_based";
# fname = "itc_rdata_cell_phasec_notw_mw_based.mat";
fname = "rdata_struct_ext_crop_phasec_notw_based"
fname_save = "rdata_extc_phasec_notw_condb"
# fname = "rdata_struct_ersp_spca";
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

#%%
for(i in 1:length(clusters)){
  cli = clusters[i];
  tt <- filter_at(indexl,vars('cluster_n'), any_vars(. %in% cli));
  #--
  subjs = unique(tt$subj_n);
  #--
  tmp_save_dir <- file.path(save_dir,sprintf("cl%i_plots",cli))
  dir.create(tmp_save_dir);
  #--
  for(l in 1:length(conds)){
    cnt = 1;
    orig_matrix <- matrix(0L,nrow=datl,ncol=length(subjs));
    mask_matrix <- matrix(0L,nrow=datl,ncol=length(subjs));
    lamb_vec <- vector(mode="double",length(subjs))
    subj_vec <- vector(mode="double",length(subjs))
    cond_vec <- vector(mode="double",length(subjs))
    clust_vec <- vector(mode="double",length(subjs))
    
    for(k in 1:length(subjs)){
      ci = conds[l];
      si = subjs[k];
      
      #%% LOAD RESULTS
      fname = sprintf("cl%i-s%i-c%i_flmoddat.RData",cli,si,ci);
      tfd <- readRDS(file=file.path(save_dir,fname)) #list(beta_hat,best_lambda)
      beta_hat <- tfd$bbeta;
      best_lambda <- tfd$blamb;
      lambda_values = tfd$lamb_vals;
      cv_errors = tfd$cv_errors;
      #--
      tmp <- get_tf_dat_mat(cli,ci,si,indexl,itc_datl)
      # tmp <- get_tf_dat(cli,ci,si,dtbl)
      tf_dat = tmp$tf_dat;
      #%% PLOT RESULTS
      # ggo <- ggplot(data.frame(lambda_values,cv_errors),
      #               aes(x=lambda_values,y=cv_errors),
      #               xlab="Lambda",
      #               ylab="Cross-Validation Error") +
      #   geom_point(pch=19) +
      #   geom_line(col="blue") +
      #   geom_vline(xintercept=best_lambda,col="red",lty=2)
      # #-- save plot
      # fname = sprintf("cl%i-s%i-c%i_cv.png",cli,si,ci);
      # ggsave(file.path(tmp_save_dir,fname),ggo)
      
      #%% TF PLOT
      # zlim_in = range(tf_dat)
      # op <- tf_plot(tf_dat,times,freqs,"Original",zlim_in)
      # #-- plot fusedlasso mask
      # tf_dat = beta_hat$beta[,1];
      # zlim_in = range(tf_dat)
      # tit_in = bquote(lambda==.(sprintf("%.3f",beta_hat$lambda[1])));
      # np <- tf_plot(tf_dat,times,freqs,tit_in,zlim_in)
      # #-- join them
      # pl <- list(op,np)
      # grid.arrange(grid::rectGrob(),grid::rectGrob())
      # ml <- marrangeGrob(pl, nrow=2, ncol=1);
      # #-- resave tf-plots?
      # fname = sprintf("cl%i-s%i-c%i_tfplot.png",cli,si,ci);
      # ggsave(file.path(tmp_save_dir,fname), ml)
      
      #%% ADD MASK TO MATRIX
      subj_vec[cnt] <- si;
      cond_vec[cnt] <- ci;
      clust_vec[cnt] <- cli;
      orig_matrix[,cnt] <- tf_dat;
      mask_matrix[,cnt] <- beta_hat$beta;
      lamb_vec[cnt] <- best_lambda;
      
      #%% ITER
      cnt = cnt + 1;
    }
    #-- save mask-mat for each condition
    fname = sprintf("cl%i-c%i_flmaskagg.RData",cli,ci); 
    saveRDS(list(orig_mat=orig_matrix,mask_mat=mask_matrix,lamb_vec=lamb_vec),file=file.path(save_dir,fname))
    
    #%% TF PLOT
    tf_dat = rowMeans(orig_matrix)
    zlim_in = range(tf_dat)
    op <- tf_plot(tf_dat,times,freqs,"Original",zlim_in)
    #-- plot fusedlasso mask
    tit_in = bquote(lambda==.(sprintf("%.3f",mean(lamb_vec))));
    tf_dat = rowMeans(mask_matrix)
    zlim_in = range(tf_dat)
    np <- tf_plot(tf_dat,times,freqs,tit_in,zlim_in)
    #-- join them
    pl <- list(op,np)
    grid.arrange(grid::rectGrob(),grid::rectGrob())
    ml <- marrangeGrob(pl, nrow=2, ncol=1);
    #-- resave tf-plots?
    fname = sprintf("allmean_cl%i-c%i_tfplot.png",cli,ci);
    ggsave(file.path(tmp_save_dir,fname), ml)
    #--
    dato = list(odat=orig_matrix,
                mmat=mask_matrix,
                lvec=lamb_vec,
                svec=subj_vec,
                cvec=cond_vec,
                clvec=clust_vec,
                times=times,
                freqs=freqs)
    fname = sprintf("allmat_cl%i-c%i.mat",cli,ci);
    R.matlab::writeMat(file.path(save_dir,fname), flasso_tbl = dato)
    # file.info(file.path(save_dir,fname))
    fname = sprintf("allmat_cl%i-c%i.RData",cli,ci);
    saveRDS(dato, file=file.path(save_dir,fname))
  }
}
