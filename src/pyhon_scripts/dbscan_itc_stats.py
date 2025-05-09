import os
import sklearn as sk
import scipy
# pip install numpy; pip install scikit-learn; pip install scipy
#%%
itc_fname = "itc_table_phasec_notw_mw_based_new.mat"
itc_fpath = "M:\\jsalminen\\GitHub\\MIND_IN_MOTION_PRJ\\_data\\MIM_dataset\\_studies\\02202025_mim_yaoa_powpow0p3_crit_speed\\__iclabel_cluster_allcond_rb3\\icrej_5\\11\\kin_eeg_step_to_step";
#--
data = scipy.io.loadmat(os.path.join(itc_fpath,itc_fname))


sk.cluster.DBSCAN()
