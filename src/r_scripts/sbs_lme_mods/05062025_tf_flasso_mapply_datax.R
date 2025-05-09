# install.packages(c("parallel","purrr","dplyr","tibble","Matrix",
#                    "lattice","gridExtra","R.devices","R.matlab","ggplot2",
#                    "gridExtra","genlasso"));
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
#%% PACKAGES FOR STATS
library(genlasso);

#%% CUSTOM FUNCTIONS
curr_dir <- getwd();
source(file.path(curr_dir,"tf_flasso_funcs.R"));

#%% LOAD DATA ============================================================== %%#
print("Loading data...");
# clusters = c(3,4,5,6,7,8,9,10,11,12,13) # RSup/RSM, PreC, LSM, Mid Cing, LSup, LPPA, RPPA
clusters = c(3,4,6,8)
fext = 'itc_rdata_table_phasec_notw_mw_based'

#%% CREATE SAVE DIR
curr_dir <- getwd();
save_dir <- paste0(curr_dir,paste0("/",fext,"_tables_figs"))
dir.create(save_dir);

#%% LOAD
mat_fpath <- paste0("/jsalminen/GitHub/MIND_IN_MOTION_PRJ/_data/MIM_dataset/_studies/02202025_mim_yaoa_powpow0p3_crit_speed/__iclabel_cluster_allcond_rb3/icrej_5/11/kin_eeg_step_to_step/",fext,".csv")

if(ispc()){
  mat_fpath <- paste0("M:",mat_fpath)
}else{
  mat_fpath <- paste0("/blue/dferris",mat_fpath);
}
#--
dtbl <- read.csv(mat_fpath)
dtbl <- filter_at(dtbl,vars('cond_n'), any_vars(. %in% c(5,6,7,8)));
dtbl <- filter_at(dtbl,vars('cluster_n'), any_vars(. %in% clusters));

#%% LOOP VARS ============================================================== %%#
# clusters=unique(dtbl$cluster_n)
tmp = get_loop_vals(dtbl,clusters)
clis = tmp$clis;
cis = tmp$cis;
sis = tmp$sis;

#%% GET SAVED DATA
cli = clusters[1]
ci = conds[1];
si = subjs[1];
tmp <- get_tf_dat(cli,ci,si,dtbl) # tmp load
datl <- length(tmp$tf_dat);
cnt = 1;
for(i in 1:length(clusters)){
  cli = clusters[i];
  tt <- filter_at(dtbl,vars('cluster_n'), any_vars(. %in% cli));
  #--
  subjs = unique(tt$subj_n);
  #--
  for(l in 1:length(conds)){
    mask_matrix <- matrix(0L,nrow=datl,ncol=length(subjs));
    lamb_vec <- vector(mode="double",length(subjs))
    for(k in 1:length(subjs)){
      ci = conds[l];
      si = subjs[k];
      
      #%% LOAD RESULTS
      fname = sprintf("cl%i-s%i-c%i_flmoddat.RData",cli,si,ci);
      tfd <- readRDS(file=file.path(save_dir,fname)) #list(beta_hat,best_lambda)
      beta_hat <- tfd$bbeta;
      best_lambda <- tfd$blamb;
      #--
      tmp <- get_tf_dat(cli,ci,si,dtbl)
      tf_dat = tmp$tf_dat;
      times = tmp$times;
      freqs = tmp$freqs;
      
      #%% TF PLOT
      # zlim_in = range(tf_dat)
      # op <- tf_plot(tf_dat,times,freqs,"Original",zlim_in)
      # #-- plot fusedlasso mask
      # zlim_in = range(beta_hat$beta)
      # tit_in = bquote(lambda==.(sprintf("%.3f",beta_hat$lambda[1])));
      # np <- tf_plot(beta_hat$beta[,1],times,freqs,tit_in,zlim_in)
      # #-- join them
      # pl <- list(op,np)
      # grid.arrange(grid::rectGrob(),grid::rectGrob())
      # ml <- marrangeGrob(pl, nrow=2, ncol=1);
      # #-- resave tf-plots?
      # fname = sprintf("cl%i-s%i-c%i_tfplot.png",clis[i],sis[i],cis[i]); 
      # ggsave(file.path(save_dir,fname), ml)
      
      #%% ADD MASK TO MATRIX
      mask_matrix[,cnt] <- beta_hat$beta;
      lamb_vec[cnt] <- best_lambda;
      
      #%% ITER
      cnt = cnt + 1;
    }
    #-- save mask-mat for each condition
    fname = sprintf("cl%i-c%i_flmask.RData",clis[i],sis[i],cis[i]); 
    saveRDS(list(fl,cv_errors,lambda_values),file=file.path(save_dir,fname))
  }
}