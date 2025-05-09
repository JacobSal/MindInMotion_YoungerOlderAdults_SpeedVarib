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
loop_items = get_mcl_dat(dtbl,clusters)
# loop_items = get_mcl_dat(dtbl,3)

#%% TEST
# saves <- lapply(loop_items,function(x) mfusedl2d(x,save_dir))

#%% RUN MPI
print("Running MPI");
if(ispc()){
  numCores <- detectCores()
}else{
  numCores = as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
}
print(numCores)
clust <- makeCluster(numCores)
system.time(saves <- mclapply(loop_items,function(x) mfusedl2d(x,save_dir),
                              mc.preschedule = FALSE,
                              mc.set.seed = FALSE))

