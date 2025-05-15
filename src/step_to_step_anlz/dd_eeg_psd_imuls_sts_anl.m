%   Project Title: MIM OA & YA SPEED & KINETICS ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/step_to_step_anlz/run_sts_dd_eeg_psd_imuls_sts_anl.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% SET WORKSPACE ======================================================= %%
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
tic
% global ADD_CLEANING_SUBMODS STUDY_DIR SCRIPT_DIR %#ok<GVMIS>
ADD_CLEANING_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        SRC_DIR = fileparts(SCRIPT_DIR); % change this if in sub folder
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        SCRIPT_DIR = getenv('SCRIPT_DIR');
        SRC_DIR = getenv('SRC_DIR');
    end
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
    end
    SRC_DIR = fileparts(SCRIPT_DIR); % change this if in scond_iub folder
end
%## Add Study, Src, & Script Paths
addpath(SCRIPT_DIR);
addpath(SRC_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
% [SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca_speed');
%% (PARAMETERS) ======================================================== %%

%% (PATHS)
%- datset name
DATA_SET = 'MIM_dataset';
%- study name
STUDY_DNAME = '02202025_mim_yaoa_powpow0p3_crit_speed';
STUDY_FNAME = 'kin_eeg_epoch_study';
ANALYSIS_DNAME = 'kin_eeg_step_to_step';
%-
studies_fpath = [PATHS.data_dir filesep DATA_SET filesep '_studies'];
%- load cluster
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
% cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb3'];
cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_allcond_rb3'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
cluster_k_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
%-
save_dir = [cluster_k_dir filesep ANALYSIS_DNAME];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
% if ~ispc
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME '_UNIX.study'],'filepath',save_dir);
% else
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME '.study'],'filepath',save_dir);
% end
%## LOAD CLUSTER STUDY
if ~ispc
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_NAME)]);
    CL_STUDY = tmp.STUDY;
else
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s.study',CLUSTER_STUDY_NAME)]);
    CL_STUDY = tmp.STUDY;
end
%## LOAD STUDY
if ~ispc
    tmp = load('-mat',[studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep sprintf('%s_UNIX.study',STUDY_FNAME)]);
    SBS_STUDY = tmp.STUDY;
else
    tmp = load('-mat',[studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep sprintf('%s.study',STUDY_FNAME)]);
    SBS_STUDY = tmp.STUDY;
end
%% (STUDY DESIGN) ====================================================== %%
ERSP_PARAMS = struct('subbaseline','off',...
    'timerange',[],...
    'ersplim',[-2,2],...
    'freqfac',4,...
    'cycles',[3,0.8],...
    'freqrange',[1,200]);
ERSP_STAT_PARAMS_COND = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','cluster',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
STUDY_DESI_PARAMS = {{'subjselect',{},...
    'variable1','cond','values1',{'flat','low','med','high'},'pairing','on',...
    'variable2','group','values2',{'H1000''s','H2000''s','H3000''s'}},...
    {'subjselect',{},...
    'variable1','cond','values1',{'0p25','0p5','0p75','1p0'},'pairing','on',...
    'variable2','group','values2',{'H1000''s','H2000''s','H3000''s'}}};
% condition_gait = unique({STUDY.datasetinfo(1).trialinfo.cond}); 
% % (09/18/2024) JS, reorders it weirdly. Manually override
%## ersp plot per cluster per condition
args = eeglab_struct2args(ERSP_STAT_PARAMS_COND);
SBS_STUDY = pop_statparams(SBS_STUDY,args{:});
args = eeglab_struct2args(ERSP_PARAMS);
SBS_STUDY = pop_erspparams(SBS_STUDY,args{:});
SBS_STUDY.cache = [];
for des_i = 1:length(STUDY_DESI_PARAMS)
    [SBS_STUDY] = std_makedesign(SBS_STUDY,[],des_i,STUDY_DESI_PARAMS{des_i}{:});
end

%## REASSIGN CLUSTER
cl_struct = par_load(cluster_k_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
SBS_STUDY.cluster = cl_struct;
[comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(SBS_STUDY);
CLUSTER_PICS = main_cl_inds;
%% (SUBJECT LOOP)
%## PARAMETERS
des_i = 2;
speed_n_chars = {'0.25','0.50','0.75','1.0'};
conds = SBS_STUDY.design(des_i).variable(1).value; %event.cond});
% conds = conds(1:4);
theta_band_lims = [4, 8];
alpha_band_lims = [8 12];
beta_band_lims  = [12 30];
alpha1_band_lims = [8,10.5];
alpha2_band_lims = [10.5,13];
beta1_band_lims = [13,20];
beta2_band_lims = [20,30];
TIME_FACTOR = 1/1000;
def_step_structs = struct('subj_char',{''}, ...
        'cond_char',{''}, ...
        'speed_char',{''}, ...
        'speed_n',[], ...
        'model_char',{''}, ...
        'model_n',[], ...
        'comp_n',[], ...
        'cluster_n',[], ...
        'comp_n_old',[], ...
        'stride_n',[], ...
        'rhs_1',[], ...
        'lto_1',[], ...
        'lhs_1',[], ...
        'rto_1',[], ...
        'rhs_2',[], ...
        'eeg_psd',[], ...
        'eeg_freqs',[], ...
        'step_dur',[], ...
        'gait_cycle_dur',[], ...
        'stance_dur',[], ...
        'swing_dur',[], ...
        'double_sup_dur',[], ...
        'single_sup_dur',[], ...
        'avg_theta',[], ...
        'avg_alpha',[], ...
        'avg_beta',[]);
cmaps_speed = linspecer(4*3);
cmaps_speed = [cmaps_speed(1,:);cmaps_speed(2,:);cmaps_speed(3,:);cmaps_speed(4,:)];
%##
% NUM_STRIDES_AVG = 36;
% STRIDES_AVG = [3,6,12,24,48];
STRIDES_AVG = [100];
% tmp_alleeg = cell(length(STUDY.datasetinfo));
conds_keep = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
% xtick_label_g = {'0.25','0.50','0.75','1.0'}; %{'0.25','0.50','0.75','1.0'};
fname_ext_store = cell(length(SBS_STUDY.datasetinfo),1);
%% SETUP PYTHON
if ~ispc
    tmp = strsplit(PATHS.src_dir,filesep);
    chk = find(strcmp(tmp,'GitHub'));
    github_dir = strjoin(tmp(1:chk),filesep);
    pyenv(Version=[github_dir filesep 'fooof_venv' filesep 'bin' filesep 'python']);
end
%% ===================================================================== %%
%## WEIRD 0-OUT SUBJECTS FOR STD_AVG_THETA
% {'H1010' },{'H1012' },{'H1013'},{'H1020' },{'H1027' },{'H1029' },{'H1035' }, ...
% {'H1048' },{'H2017' },{'H2020' },{'H3029' },{'H3072' },{'H3103' ,{'H3120' }, ...
% {'NH3008'},{'NH3021'},{'NH3043'},{'NH3058'},{'NH3066'},{'NH3070'},{'NH3112'}
%(02/04/2025) JS, these have been fixed by removing a bug from
%b_epoch_eeg_kin.m script where the cleaning mask was not appropriately
%being applied to the data.

%## MAIN LOOP
% for subj_i = 1:length(STUDY.datasetinfo)
parfor subj_i = 1:length(SBS_STUDY.datasetinfo)      
    tmp_sbs_study = SBS_STUDY;
    tmp_cl_study = CL_STUDY;
    tmp_fname_exts = cell(length(STRIDES_AVG),1);
    fooof_freqs = [];
    ff = 1;
    log_s = false;
    fr = struct.empty;
    % fname_ext = cell(1,5);
    %## LOAD EEG DATA
    try
        EEG = pop_loadset('filepath',tmp_sbs_study.datasetinfo(subj_i).filepath,'filename',tmp_sbs_study.datasetinfo(subj_i).filename);
        subj_ind = strcmp(EEG.subject,{tmp_cl_study.datasetinfo.subject});
        fprintf('Running subject %s...\n',EEG.subject);     
        %-- load icaact
        % tmpf = strsplit(tmp_sbs_study.datasetinfo(subj_i).filename,'.');
        % tmpf{2} = 'fdt';
        % tmp_dat = strjoin(tmpf,'.');
        % fprintf('Running Subject %s\n',EEG.subject);
        % EEG = eeg_checkset(EEG,'loaddata');
        % if isempty(EEG.icaact)
        %     fprintf('%s) Recalculating ICA activations\n',EEG.subject);
        %     EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        %     EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        % end        
    catch e
        fprintf('%s',getReport(e));
        continue
    end
    if ~any(subj_ind)
        continue
    end 
    %- sanity check
    %{
    pop_eegplot( EEG, 1, 1, 1);
    %}
    %## LOAD PSD DATA
    try
        epoched_fPath = strsplit(tmp_sbs_study.datasetinfo(subj_i).filepath,filesep);
        % epoched_fPath = strsplit(tmp_cl_study.datasetinfo(subj_i).filepath,filesep);
        icatimef_f = [strjoin(epoched_fPath,filesep) filesep sprintf('%s.icaspec',EEG.subject)];    
        for nn = 1:length(STRIDES_AVG)
            NUM_STRIDES_AVG = STRIDES_AVG(nn);
            fprintf('Running number of strides averaging: %i...\n',NUM_STRIDES_AVG);
            %- load .icatimef load-in parameters
            fprintf('Loading Time-Frequency Data...\n');
            tmpf = load(icatimef_f,'-mat');
            %-
            fn = fieldnames(tmpf);
            inds = find(contains(fieldnames(tmpf),'comp'));
            test = tmpf.(fn{inds(1)});
            %-
            eeg_psd = zeros(size(test,1),size(test,2),length(inds),'double');
            for i = 1:length(inds)
                eeg_psd(:,:,i) = tmpf.(fn{inds(i)}); % freq x epoch x chan
            end
            %-    
            freqs_orig = tmpf.freqs;
            trialinfo = tmpf.trialinfo;
            
            %## SUBSET DATA TO CONDS_KEEP
            inds_cond = cellfun(@(x) any(strcmp(x,conds_keep)),{trialinfo.cond});
            trialinfo = trialinfo(inds_cond);
            fprintf('Using conditions: %s\n',strjoin(unique({trialinfo.cond}),','));
            eeg_psd = eeg_psd(:,inds_cond,:);
            nolog_eeg_psd = 10.^(eeg_psd/10);
        
            %## RUN FOOOF
            fprintf('==== running FOOOF ====\n');
            settings = struct('peak_width_limits',[1,8],...
                'min_peak_height',0.05,...
                'max_n_peaks',5);
            %-- general notes on fooof
            % *** fr.gaussian_params is a 2d array where each row corresponds to a peak
            % and the columns are [mean, hight, standard deviation] of that peak.
            % *** fr.peak_params is a 2d array where each row corresponds to a
            % peak's fitting parameters [CF, PW, BW] (see.
            % https://www.nature.com/articles/s41593-020-00744-x, methods for more
            % details)
            f_range = [3, 40];
            f_ind = find(freqs_orig > f_range(1) & freqs_orig < f_range(2));
            f_ind = sort([f_ind; min(f_ind)-1; max(f_ind)+1]);
            return_model = true;    
            %-- stores
            tmp_psd = zeros(length(f_ind),size(eeg_psd,2),size(eeg_psd,3));    
            tmp_psd_std = zeros(length(f_ind),size(eeg_psd,2),size(eeg_psd,3));
            tmp_ap_mu = zeros(2,size(eeg_psd,2),size(eeg_psd,3));
            tmp_ap_std = zeros(2,size(eeg_psd,2),size(eeg_psd,3));
            tmp_ap_cov = zeros(2,size(eeg_psd,2),size(eeg_psd,3));
            tmp_covt = zeros(1,size(eeg_psd,2),size(eeg_psd,3));
            tmp_cova = zeros(1,size(eeg_psd,2),size(eeg_psd,3));
            tmp_covb = zeros(1,size(eeg_psd,2),size(eeg_psd,3));
            tmp_qcvt = zeros(1,size(eeg_psd,2),size(eeg_psd,3));
            tmp_qcva = zeros(1,size(eeg_psd,2),size(eeg_psd,3));
            tmp_qcvb = zeros(1,size(eeg_psd,2),size(eeg_psd,3));
            %-
            % tmp_conds = unique({trialinfo.trial_num_code});
            % tmp_conds = tmp_conds(~cellfun(@isempty,tmp_conds));
            %-
            tmp_conds = unique({trialinfo.cond});
            %-
            fname_ext = cell(1,5);
            def_cond_struct = struct('ind',[], ...
                'cond',{''}, ...
                'indices',[]);
            cond_struct = def_cond_struct;
            ttim = tic();
            for ct = 1:size(nolog_eeg_psd,3)
                %## GET DATA
                tmpp_raw = zeros(length(f_ind),size(eeg_psd,2),1); 
                tmpap_raw = zeros(2,size(eeg_psd,2),1); 
                tmpp_mu = zeros(length(f_ind),size(eeg_psd,2),1); 
                tmpp_std = zeros(length(f_ind),size(eeg_psd,2),1); 
                tmpap_mu = zeros(2,size(eeg_psd,2),1); 
                tmpap_std = zeros(2,size(eeg_psd,2),1);
                tmpap_cov = zeros(2,size(eeg_psd,2),1);
                tmpp_cov_theta = zeros(1,size(eeg_psd,2),1);
                tmpp_cov_alpha = zeros(1,size(eeg_psd,2),1);
                tmpp_cov_beta = zeros(1,size(eeg_psd,2),1);
                tmpp_qcv_theta = zeros(1,size(eeg_psd,2),1);
                tmpp_qcv_alpha = zeros(1,size(eeg_psd,2),1);
                tmpp_qcv_beta = zeros(1,size(eeg_psd,2),1);
                fprintf('Processing channel index %i...\n',ct)
        
                %## REMOVE FOOOF OF EACH STRIDE
                for i = 1:size(nolog_eeg_psd,2)
                    fr = fooof(freqs_orig(f_ind),squeeze(nolog_eeg_psd(f_ind,i,ct)),f_range,settings,return_model);
                    %-- assign data
                    tmpp_raw(:,i) = 10*(log10(squeeze(nolog_eeg_psd(f_ind,i,ct)))) - 10*(fr.ap_fit');   
                    tmpap_raw(1:2,i) = fr.aperiodic_params'; 
                end
                if ct == 1
                    fname_ext{ff} = 'perstridefb';
                    ff = ff + 1;
                end
                fooof_freqs = fr.freqs;
                % %-- cond tracking struct
                % cond_struct = repmat(def_cond_struct,[1,size(ext_tmp,2)]);        
                % for i = 1:length(trialinfo)
                %     cond_struct(i).ind = i;
                %     cond_struct(i).cond = trialinfo(i).cond;
                %     cond_struct(i).indices = i;
                % end    
        
                %## SLIDIN AVERAGE (NO FOOOF)
                %-- band defs
                tband = (fooof_freqs > theta_band_lims(1) & fooof_freqs < theta_band_lims(2));
                aband = (fooof_freqs > alpha_band_lims(1) & fooof_freqs < alpha_band_lims(2));
                bband = (fooof_freqs > beta_band_lims(1) & fooof_freqs < beta_band_lims(2));
                %--
                cnt = 1;
                cond_struct = repmat(def_cond_struct,[1,size(tmpp_raw,2)]);
                for c_i = 1:length(tmp_conds)
                    %-
                    inds_cond = find(cellfun(@(x) strcmp(x,tmp_conds{c_i}),{trialinfo.cond}));
                    fooof_tmp = tmpp_raw(:,inds_cond);
                    ap_tmp = tmpap_raw(:,inds_cond);
                    %-
                    slides = 1:NUM_STRIDES_AVG:size(fooof_tmp,2);
                    % slides = unique([slides,size(fooof_tmp,2)]);
                    %(01/31/2025) JS, removing the edge to prevent zeroed out
                    %standard deviation values due to singular extractions.
                    %- fooof
                    for i = 1:length(slides)-1  
                        spec_in = squeeze(fooof_tmp(:,slides(i):(slides(i+1)-1)));
                        %-- mean std calcs
                        tmpp_mu(:,cnt) = mean(spec_in,2);
                        % tmpp_mu(:,cnt) = median(spec_in,2);
                        % tmpp_std(:,cnt) = std(spec_in,[],2)';
                        tmpp_std(:,cnt) = std(spec_in,[],2)'/sqrt(NUM_STRIDES_AVG);
                        tmpap_mu(1,cnt) = mean(ap_tmp(1,slides(i):(slides(i+1)-1)),2);
                        % tmpap_std(1,cnt) = std(ap_tmp(1,slides(i):(slides(i+1)-1)),[],2);
                        tmpap_std(1,cnt) = std(ap_tmp(1,slides(i):(slides(i+1)-1)),[],2)/sqrt(NUM_STRIDES_AVG);
                        tmpap_mu(2,cnt) = mean(ap_tmp(2,slides(i):(slides(i+1)-1)),2);
                        % tmpap_std(2,cnt) = std(ap_tmp(2,slides(i):(slides(i+1)-1)),[],2);
                        tmpap_std(2,cnt) = std(ap_tmp(2,slides(i):(slides(i+1)-1)),[],2)/sqrt(NUM_STRIDES_AVG);
                        
                        %## COV
                        %-- theta
                        tmp = squeeze(mean(spec_in(tband,:),1));
                        % tmp = squeeze(median(spec_in(tband,:),1));
                        tmpp_cov_theta(1,cnt) = 100*(std(tmp,[],2)/abs(mean(tmp,2))); % theta
                        %-- beta
                        tmp = squeeze(mean(spec_in(bband,:),1));
                        % tmp = squeeze(median(spec_in(bband,:),1));
                        tmpp_cov_beta(1,cnt) = 100*(std(tmp,[],2)/abs(mean(tmp,2))); % beta
                        %-- alpha
                        tmp = squeeze(mean(spec_in(aband,:),1));
                        % tmp = squeeze(median(spec_in(aband,:),1));
                        tmpp_cov_alpha(1,cnt) = 100*(std(tmp,[],2)/abs(mean(tmp,2))); % alpha
                        %-- ap meas.
                        tmpap_cov(2,cnt) = tmpap_std(2,cnt)/tmpap_mu(2,cnt);
                        tmpap_cov(1,cnt) = tmpap_std(1,cnt)/tmpap_mu(1,cnt); 

                        %## QCV
                        tmp = squeeze(median(spec_in(tband,:),1));
                        tmp50 = squeeze(prctile(tmp,50));                        
                        tmp25 = squeeze(prctile(tmp,25));
                        tmp75 = squeeze(prctile(tmp,75));  
                        % tmpp_qcv_theta(1,cnt) = ((tmp75-tmp25)/tmp50)*100; % theta
                        tmpp_qcv_theta(1,cnt) = ((tmp75-tmp25)/(tmp75+tmp25))*100; % theta
                        tmp = squeeze(median(spec_in(bband,:),1));
                        tmp50 = squeeze(prctile(tmp,50));                       
                        tmp25 = squeeze(prctile(tmp,25));
                        tmp75 = squeeze(prctile(tmp,75)); 
                        % tmpp_qcv_beta(1,cnt) = ((tmp75-tmp25)/tmp50)*100; % beta
                        tmpp_qcv_beta(1,cnt) = ((tmp75-tmp25)/(tmp75+tmp25))*100; % beta
                        tmp = squeeze(median(spec_in(aband,:),1));
                        tmp50 = squeeze(prctile(tmp,50));                        
                        tmp25 = squeeze(prctile(tmp,25));
                        tmp75 = squeeze(prctile(tmp,75)); 
                        % tmpp_qcv_alpha(1,cnt) = ((tmp75-tmp25)/tmp50)*100; % alpha
                        tmpp_qcv_alpha(1,cnt) = ((tmp75-tmp25)/(tmp75+tmp25))*100; % alpha
                        %--
                        if ct == 1 && cnt == 1
                            fname_ext{ff} = sprintf('apfix_std_mi_nfslidingb%i',NUM_STRIDES_AVG);
                            ff = ff + 1;
                        end
                        cond_struct(cnt).cond = tmp_conds{c_i};
                        cond_struct(cnt).indices = inds_cond(slides(i):slides(i+1)-1);
                        cond_struct(cnt).ind = cnt;
                        cnt = cnt + 1;
                    end
                end
                %(01/17/2025) JS, This baselining seems most similar to the group
                %average speed results
                %(01/29/2025) JS, average sliding fooof as recommended by
                %arkaprava. Increases homogeneity of the PSD measure. I think
                %this also greatly reduces the variability of the PSDs without
                %sacrificing averages
                
                %## ASSIGN DATA
                chk = all(tmpp_mu==0);
                tmpp_mu = tmpp_mu(:,~chk);        
                tmp_psd(:,1:size(tmpp_mu,2),ct) = tmpp_mu;
                %--
                chk = all(tmpp_std==0);
                tmpp_std = tmpp_std(:,~chk);  
                tmp_psd_std(:,1:size(tmpp_std,2),ct) = tmpp_std;
                %--
                chk = all(tmpap_mu==0);
                tmpap_mu = tmpap_mu(:,~chk); 
                tmp_ap_mu(:,1:size(tmpap_mu,2),ct) = tmpap_mu;
                %--
                chk = all(tmpap_std==0);
                tmpap_std = tmpap_std(:,~chk); 
                tmp_ap_std(:,1:size(tmpap_std,2),ct) = tmpap_std;
                %--
                chk = all(tmpap_cov==0);
                tmpap_cov = tmpap_cov(:,~chk); 
                tmp_ap_cov(:,1:size(tmpap_cov,2),ct) = tmpap_cov;
                %--
                chk = (tmpp_cov_alpha==0);
                tmpp_cov_alpha = tmpp_cov_alpha(:,~chk); 
                tmp_cova(:,1:size(tmpp_cov_alpha,2),ct) = tmpp_cov_alpha;
                %--
                chk = (tmpp_cov_beta==0);
                tmpp_cov_beta = tmpp_cov_beta(:,~chk); 
                tmp_covb(:,1:size(tmpp_cov_beta,2),ct) = tmpp_cov_beta;
                %--
                chk = (tmpp_cov_theta==0);
                tmpp_cov_theta = tmpp_cov_theta(:,~chk); 
                tmp_covt(:,1:size(tmpp_cov_theta,2),ct) = tmpp_cov_theta;
                %--
                chk = (tmpp_qcv_alpha==0);
                tmpp_qcv_alpha = tmpp_qcv_alpha(:,~chk); 
                tmp_qcva(:,1:size(tmpp_qcv_alpha,2),ct) = tmpp_qcv_alpha;
                %--
                chk = (tmpp_qcv_beta==0);
                tmpp_qcv_beta = tmpp_qcv_beta(:,~chk); 
                tmp_qcvb(:,1:size(tmpp_qcv_beta,2),ct) = tmpp_qcv_beta;
                %--
                chk = (tmpp_qcv_theta==0);
                tmpp_qcv_theta = tmpp_qcv_theta(:,~chk); 
                tmp_qcvt(:,1:size(tmpp_qcv_theta,2),ct) = tmpp_qcv_theta;
            end
        
            %## CLEAN UP DATA ARRAYS
            fprintf('FOOOF Process done: %0.2g.\n',toc(ttim));
            inds_store = cell(size(nolog_eeg_psd,3),1);
            for tci = 1:size(nolog_eeg_psd,3)
                inds_store{tci} = find(~all(squeeze(tmp_psd(:,:,tci))==0,1));
                fprintf('%i) number of epochs retained: %i\n',tci,length(inds_store{tci}));
            end
            
            %-- delete 0'd events
            tmp_psd = tmp_psd(:,inds_store{1},:);
            tmp_psd_std = tmp_psd_std(:,inds_store{1},:);
            tmp_ap_mu = tmp_ap_mu(:,inds_store{1},:);
            tmp_ap_std = tmp_ap_std(:,inds_store{1},:);
            tmp_ap_cov = tmp_ap_cov(:,inds_store{1},:);
            tmp_cova = tmp_cova(:,inds_store{1},:);
            tmp_covb = tmp_covb(:,inds_store{1},:);
            tmp_covt = tmp_covt(:,inds_store{1},:);
            tmp_qcva = tmp_qcva(:,inds_store{1},:);
            tmp_qcvb = tmp_qcvb(:,inds_store{1},:);
            tmp_qcvt = tmp_qcvt(:,inds_store{1},:);
            cond_struct = cond_struct(inds_store{1});        
            %-
            % inds_keep = find(~all(squeeze(tmp_psd(:,:,3))==0,1));
            % tmp_psd = tmp_psd(:,inds_keep,:);
            % trialinfo = trialinfo(inds_keep);        
            fprintf('considering %i trials, spanning conditions: %s\n',length(trialinfo),strjoin(unique({trialinfo.cond}),','));
            eeg_psd = tmp_psd;
            freqs_orig = fooof_freqs;
    
            %## GET CLUSTER DATA
            fprintf('Getting cluster information...\n');
            tmp_subjs = {tmp_cl_study.datasetinfo.subject};
            comp_arr = zeros(3,length(CLUSTER_PICS)); %(1, comps; 2, old_comps; 3, cluster
            for i = 1:length(CLUSTER_PICS)
                cl_i = CLUSTER_PICS(i);
                ind = find(strcmp(EEG.subject,tmp_subjs));
                ind = tmp_sbs_study.cluster(cl_i).sets == ind;
                if any(ind)
                    comp_arr(1,i) = tmp_sbs_study.cluster(cl_i).comps(ind);
                    comp_arr(3,i) = cl_i;
                    comp_arr(2,i) = EEG.etc.urreject.ic_keep(comp_arr(1,i));
                end
            end
            comp_arr = comp_arr(:,all(comp_arr,1));
        
            %## PSD DATA SAVE
            psd_struct = struct('psd_mean',tmp_psd, ...
                'psd_std',tmp_psd_std, ...
                'ap_mean',tmp_ap_mu, ...
                'ap_std',tmp_ap_std, ...
                'ap_cov',tmp_ap_cov, ...
                'cond_struct',cond_struct, ...
                'freqs',fooof_freqs, ...
                'indx_to_comp',comp_arr, ...
                'indx_to_comp_dims',{{'new_comp_num','old_comp_num','clust_num'}});
            tmp = fname_ext(~cellfun(@isempty,fname_ext));
            par_save(psd_struct,EEG.filepath,sprintf('psd_output_%s.mat',strjoin(tmp,'_')));
            psd_struct = struct.empty;
    
            %## PSD VALIDATION PLOT
            %{
            %- colors
            cmaps_terrain = linspecer(4);
            custom_yellow = [254,223,0]/255;
            cmaps_terrain = [cmaps_terrain(3,:);custom_yellow;cmaps_terrain(4,:);cmaps_terrain(2,:)];
            cmaps_speed = linspecer(4*3);
            cmaps_speed = [cmaps_speed(1,:);cmaps_speed(2,:);cmaps_speed(3,:);cmaps_speed(4,:)];
            params = struct('cmaps',[cmaps_speed;cmaps_terrain], ...
                'freqs',fooof_freqs, ...
                'xtick_label',{{'0.25','0.50','0.75','1.0','flat','low','med','high'}});
            % conds = {'0p25','0p5','0p75','1p0','flat','low','med','high'}; %unique({cond_struct.cond})
            clusts = (3:13);
            tmp_conds = {'0p25','0p5','0p75','1p0'}; %unique({cond_struct.cond});
            params.xtick_label = tmp_conds;
            %-- loop through comps
            tmp_fext = fname_ext(~cellfun(@isempty,fname_ext));
            mkdir([EEG.filepath filesep strjoin(tmp_fext,'_')])
            for i = 1:3 %length(clusts)
                cl_i = clusts(i);
                comp_n = comp_arr(1,comp_arr(3,:) == cl_i);
                if ~isempty(comp_n)
                    % psd_dat_in = squeeze(tmp_psd_std(:,:,comp_n));
                    % fig = figure;
                    % ax = axes();
                    % local_psd_plot_valid(ax,psd_dat_in,cond_struct,tmp_conds,params);
                    % exportgraphics(fig,[EEG.filepath filesep strjoin(tmp_fext,'_') filesep sprintf('%i_stdpsd_plot.png',cl_i)])
                    % close(fig);
                    %--
                    psd_dat_in = squeeze(tmp_psd(:,:,comp_n));
                    fig = figure;
                    ax = axes();
                    local_psd_plot_valid(ax,psd_dat_in,cond_struct,tmp_conds,params);
                    exportgraphics(fig,[EEG.filepath filesep strjoin(tmp_fext,'_') filesep sprintf('%i_mupsd_plot.png',cl_i)])
                    % close(fig);
                end
            end
            %}
            
            %## EXAMPLE PIPELINE FIGURE
            %{
            PLOT_STRUCT = struct( ...
                'ylim',[],...
                'y_label',{'10*log_{10}(PSD)'},...
                'do_set_ax_props',true, ...
                'y_label_props',{{ ...
                    'Units','Normalized', ...
                    'FontSize',8, ...
                    'FontWeight','bold', ...
                    }}, ...
                'xlim',[fooof_freqs(1),fooof_freqs(end)],...
                'x_label',{'Frequency (Hz)'},...
                'x_label_props',{{ ...
                    'Units','Normalized', ...
                    'FontSize',8, ...
                    'FontWeight','bold', ...
                    }}, ...
                'x_label_yoffset',-0.1,...
                'xtick_labs',{{}}, ...
                'xticks',[], ...
                'xtick_angle',45, ...
                'title',{{''}},...
                'title_props',{{ ...
                    'Units','Normalized', ...
                    'FontSize',8, ...
                    'FontWeight','bold', ...
                    }}, ...
                'ax_props',{{...
                    'box','off', ...
                    'LineWidth',2, ...
                    'FontName','Arial', ...
                    'FontSize',8, ...
                    'OuterPosition',[0 0 1 1], ...
                    'Position',[0.125,0.125,0.7,0.7]}});
            LINE_STRUCT = struct('do_line_avg',false, ...
                'line_props',{{ ...
                    'LineWidth',2, ...
                    'LineStyle','-', ...
                    'DisplayName','line', ...
                    'Color',[0.5,0.5,0.5,0.5], ...
                    }}, ...
                'line_avg_fcn',@(x) mean(x,1), ...
                'do_err_shading',true, ...
                'err_props',{{ ...
                    'LineStyle',':', ...
                    'LineWidth',2, ...
                    'FaceAlpha',0.6, ...
                    'EdgeColor','none', ...
                    'FaceColor',[0.5,0.5,0.5]}}, ...
                'err_bnd_vec',[]);
            %--
            i = 1;
            c_i = 1;
            cl_i = clusts(i);
            comp_n = size(nolog_eeg_psd,3); %comp_arr(1,comp_arr(3,:) == cl_i);
            % %--
            % fig = figure;
            % ax = axes();
            for c_i = 1 %:length(tmp_conds)
                %%
                fig = figure;
                ax = axes();
                tmp_plot_struct = PLOT_STRUCT;
                tmp_line_struct = LINE_STRUCT;
                inds_cond = cellfun(@(x) any(strcmp(x,tmp_conds{c_i})),{cond_struct.cond});
                tmp_c = cond_struct(inds_cond);
                inds_cond = find(inds_cond);
                %--
                psd_in_c = tmp_psd(:,inds_cond(1),comp_n);                
                % psd_in_c = tmp_psd_std(:,inds_cond(1),comp_n);
                %--
                % psd_in_c = squeeze(tmpp_raw(:,[tmp_c.indices]));
                %--
                ap_cos = squeeze(tmpap_raw(:,[tmp_c.indices]));
                chansi = randi(size(ap_cos,2),[1,1]);
                ap_cos = ap_cos(:,chansi);
                %--
                psd_in_c = 10.*log10(squeeze(nolog_eeg_psd(f_ind,[tmp_c.indices],comp_n)));
                psd_in_c = psd_in_c(:,chansi,:);
                %--
                psd_ap = zeros(size(psd_in_c,1),3);
                % expm = 1;
                % km = 1;
                % offm = 1;
                % for cc = 1:size(ap_cos,2)
                %     % psd_ap(:,cc) = 10*(ap_cos(1,cc)-log10(params.freqs.^(ap_cos(2,cc))));
                %     psd_ap(:,cc) = 10*(ap_cos(1,cc)*offm-log10(km+params.freqs.^(ap_cos(2,cc)*expm)));
                % end
                cc=1;
                expm = 1;
                km = 0;
                offm = 1;
                psd_ap(:,cc) = 10*(ap_cos(1,1)*offm-log10(km+params.freqs.^(ap_cos(2,1)*expm)));
                cc=2;
                expm = 1.3;
                km = 0;
                offm = 1;
                psd_ap(:,cc) = 10*(ap_cos(1,1)*offm-log10(km+params.freqs.^(ap_cos(2,1)*expm)));
                cc=3;
                expm = 1;
                km = 0;
                offm = 1.3;
                psd_ap(:,cc) = 10*(ap_cos(1,1)*offm-log10(km+params.freqs.^(ap_cos(2,1)*expm)));
                %--
                % psd_in_c = psd_in_c(:,1:36);
                % psd_ap = psd_ap(:,1:36)
                %--
                % mui = psd_in_c; 
                % % mui = mean(psd_in_c,2);
                % stdi = std(psd_in_c,[],2);
                % tmp_line_struct.do_err_shading = false;
                % tmp_line_struct.err_bnd_vec = [mui+stdi/sqrt(size(stdi,2)), ...
                %     mui-stdi/sqrt(size(stdi,2))];
                % %--
                % % tmp_line_struct.do_err_shading = false;
                % tmp_line_struct.line_props = { ...
                %     'LineWidth',2, ...
                %     'LineStyle','-', ...
                %     'DisplayName',params.xtick_label{c_i}, ...
                %     'Color',[params.cmaps(c_i,:),0.6], ...
                %     };
                % % tmp_line_struct.line_color = [cmaps(c_i,:),0.65];
                % tmp_line_struct.err_props = { ...
                %     'LineStyle',':', ...
                %     'LineWidth',3, ...
                %     'FaceAlpha',0.6, ...
                %     'EdgeColor','none', ...
                %     'FaceColor',params.cmaps(c_i,:)};
                % %## PLOT PSD
                % [~,~,Li] = plot_psd(ax,mui,params.freqs, ...
                %     'LINE_STRUCT',tmp_line_struct, ...
                %     'PLOT_STRUCT',tmp_plot_struct);
                %--
                % mui = psd_ap;
                % % mui = mean(psd_ap,2);
                % stdi = std(psd_ap,[],2);
                % tmp_line_struct.do_err_shading = false;
                % tmp_line_struct.err_bnd_vec = [mui+stdi/sqrt(size(stdi,2)), ...
                %     mui-stdi/sqrt(size(stdi,2))];
                % tmp_line_struct.line_props = { ...
                %     'LineWidth',2, ...
                %     'LineStyle','--', ...
                %     'DisplayName',params.xtick_label{c_i}, ...
                %     'Color',[params.cmaps(c_i,:),0.6], ...
                %     };
                % [~,~,Li] = plot_psd(ax,mui,params.freqs, ...
                %     'LINE_STRUCT',tmp_line_struct, ...
                %     'PLOT_STRUCT',tmp_plot_struct);
                % hold on;
                % ylim([-30,-5]);
                
                %##
                cmaps = linspecer(3);
                for cc= 1:3
                    % psd_ap(:,cc) = 10*(ap_cos(1,1)*offm-log10(km+params.freqs.^(ap_cos(2,1)*expm)));
                    mui = psd_ap(:,cc);
                    tmp_line_struct.do_err_shading = false;
                    tmp_line_struct.line_props = { ...
                        'LineWidth',2, ...
                        'LineStyle','--', ...
                        'DisplayName',sprintf('%i',cc), ...
                        'Color',[cmaps(cc,:),0.6], ...
                        };
                    [~,~,Li] = plot_psd(ax,mui,params.freqs, ...
                        'LINE_STRUCT',tmp_line_struct, ...
                        'PLOT_STRUCT',tmp_plot_struct);
                end
                legend();
                hold on;
                ylim([-30,-5]);
                % title(sprintf('offm=%0.1f, km=%0.1f, expm=%0.1f',offm,km,expm));
            end
            %--
            exportgraphics(fig,[EEG.filepath filesep strjoin(tmp_fext,'_') filesep sprintf('%i_ap_fit_1stride_ex4-5.pdf',cl_i)], ...
                'ContentType','Vector');
            exportgraphics(fig,[EEG.filepath filesep strjoin(tmp_fext,'_') filesep sprintf('%i_ap_fit_1stride_ex4-5.tiff',cl_i)], ...
                'Resolution',900);
            %}
            
            %%
            %## EXTRACT IMU DATA
            fprintf('Performing gait kinematic calculations...\n');
            tband = (freqs_orig > theta_band_lims(1) & freqs_orig < theta_band_lims(2));
            aband = (freqs_orig > alpha_band_lims(1) & freqs_orig < alpha_band_lims(2));
            bband = (freqs_orig > beta_band_lims(1) & freqs_orig < beta_band_lims(2));
            body_posx = find(contains({EEG.chanlocs.labels},'body_xpos','IgnoreCase',true));
            body_posy = find(contains({EEG.chanlocs.labels},'body_ypos','IgnoreCase',true));
            body_posz = find(contains({EEG.chanlocs.labels},'body_zpos','IgnoreCase',true));
            EEG = eeg_checkset(EEG,'loaddata');
            %##
            steps_struct = repmat(def_step_structs,[1,length(EEG.epoch)*size(comp_arr,2)]);
            cnt = 1;
            %- loop through each condition
            for cond_i = 1:length(conds)
                %##
                % inds = find(cellfun(@(x) any(strcmp(x,conds{cond_i})),{EEG.epoch.eventcond}));
                inds_cond = find(cellfun(@(x) any(strcmp(x,conds{cond_i})),{trialinfo.cond}));
                %- imu data
                datx = squeeze(EEG.data(body_posx,:,inds_cond)); % AP (anteroposterior)
                daty = squeeze(EEG.data(body_posy,:,inds_cond)); % ML (mediolateral)
                datz = squeeze(EEG.data(body_posz,:,inds_cond)); % UD (Up-Down)
                %-- validation figures
                % fh = imu_valid_plots(EEG, ...
                %     EEG.filepath, ...
                %     sprintf('kindatgen_%s_',conds{cond_i}), ...
                %     'EXPORT_RES',100);
                % close(fh)
                %## ADVANCED EPOCH'ING?
                %- before & after inds
                do_log_s = zeros(1,length(inds_cond));
                tmp_step_struct = repmat(def_step_structs,[1,length(inds_cond)]);        
                % for s_i = 2:length(inds_cond)-1
                for s_i = 1:length(inds_cond)
                    %## DETERMINE LS FOOT EVENTS for the event and surrounding           
                    tmp = EEG.epoch(inds_cond(s_i));
                    rhs_1 = nan();
                    rhs_2 = nan();
                    lhs_1 = nan();
                    lto_1 = nan();
                    rto_1 = nan();
                    %## 
                    t_rhi = find(strcmp(tmp.eventtype,'RHS'));
                    lhi = find(strcmp(tmp.eventtype,'LHS'));
                    rti = find(strcmp(tmp.eventtype,'RTO'));
                    lti = find(strcmp(tmp.eventtype,'LTO'));
                    %##
                    good_s = true;
                    ii = 0;
                    while good_s
                        if length(t_rhi) > 2
                            rhi = [t_rhi(1+ii),t_rhi(2+ii)];
                        else
                            rhi = t_rhi;
                        end
                        rht = [tmp.eventlatency{rhi}]*TIME_FACTOR;
                        lht = [tmp.eventlatency{lhi}]*TIME_FACTOR;
                        rtt = [tmp.eventlatency{rti}]*TIME_FACTOR;
                        ltt = [tmp.eventlatency{lti}]*TIME_FACTOR;
                        %-
                        rhs_1 = rht(1);
                        rhs_2 = rht(2);
                        %- latency criteria
                        ind = lht > rht(1) & lht < rht(2);
                        lhs_1 = lht(ind);
                        ind = ltt > rht(1) & ltt < rht(2);
                        lto_1 = ltt(ind);
                        ind = rtt > rht(1) & rtt < rht(2);
                        rto_1 = rtt(ind);
                        %- flexible solve if weird duplicate events
                        if length(lhs_1) > 1
                            lhs_1 = lhs_1(1);
                        end
                        if length(lto_1) > 1
                            lto_1 = lto_1(1);
                        end
                        if length(rto_1) > 1
                            rto_1 = rto_1(1);
                        end
                        %(01/29/2025) JS, adding this in the case of weird events
                        %that have duplicate sequences. Seems reasonable given the
                        %timings seem to match better if you pick the firsts in the
                        %sequence.
                        %-
                        if all(~cellfun(@isempty,{rhs_1,lhs_1,lto_1,rto_1,rhs_2})) && ...
                                length([rhs_1,lhs_1,lto_1,rto_1,rhs_2]) == 5
                            %-
                            log_s = true;
                            good_s = false;
                            fprintf('x');
                        else
                            % fprintf('Error. Trying next consecutive strikes...\n');
                            fprintf('.');
                            ii = ii + 1;
                            if ii+2 > length(t_rhi)
                                fprintf('o');
                                good_s = false;
                                log_s = false;
                            end
                        end
                    end
                    %## CALCULATE IMU MEASURES
                    %- rhs_1, lto_1, lhs_1, rto_1, rhs_2
                    EVENT_ERR = 2.001; % potentially a millisecond discrepency for some reason
                    if log_s
                        %- rhs_1 ind
                        diff = EEG.times - rhs_1/TIME_FACTOR;
                        rhsi_1 = find(diff > -EVENT_ERR & diff < EVENT_ERR,1,'first');
                        %- rhs_1 ind
                        diff = EEG.times - lto_1/TIME_FACTOR;
                        ltoi_1 = find(diff > -EVENT_ERR & diff < EVENT_ERR,1,'first');
                        %- rhs_1 ind
                        diff = EEG.times - lhs_1/TIME_FACTOR;
                        lhsi_1 = find(diff > -EVENT_ERR & diff < EVENT_ERR,1,'first');
                        %- rhs_1 ind
                        diff = EEG.times - rto_1/TIME_FACTOR;
                        rtoi_1 = find(diff > -EVENT_ERR & diff < EVENT_ERR,1,'first');
                        %- rhs_2 ind
                        diff = EEG.times - rhs_2/TIME_FACTOR;
                        rhsi_2 = find(diff > -EVENT_ERR & diff < EVENT_ERR,1,'first');
                        %- calculate pelvis movement
                        %* step width (ml exc)
                        vec = [daty(rhsi_1,s_i),daty(lhsi_1,s_i)];
                        m1 = min(vec);
                        m2 = max(vec);
                        ml_exc = sqrt((m2-m1)^2);
                        %* ap exc
                        vec = [datx(rhsi_1,s_i),datx(lhsi_1,s_i)];
                        m1 = min(vec);
                        m2 = max(vec);
                        ap_exc = sqrt((m2-m1)^2);
                        %* ud exc
                        vec = [datz(rhsi_1,s_i),datz(lhsi_1,s_i)];
                        m1 = min(vec);
                        m2 = max(vec);
                        ud_exc = sqrt((m2-m1)^2);
                        %* min max exc (ml)
                        vec1 = [daty(rhsi_1:lhsi_1,s_i)];
                        vec2 = [daty(lhsi_1:rhsi_2,s_i)];
                        m1 = mean([min(vec1),min(vec2)]);
                        m2 = mean([max(vec1),max(vec2)]);
                        ml_mm_exc = sqrt((m2-m1)^2);
                        %* min max exc (ap)
                        vec1 = [datx(rhsi_1:lhsi_1,s_i)];
                        vec2 = [datx(lhsi_1:rhsi_2,s_i)];
                        m1 = mean([min(vec1),min(vec2)]);
                        m2 = mean([max(vec1),max(vec2)]);
                        ap_mm_exc = sqrt((m2-m1)^2);
                        %* min max exc (ud)
                        vec1 = [datz(rhsi_1:lhsi_1,s_i)];
                        vec2 = [datz(lhsi_1:rhsi_2,s_i)];
                        m1 = mean([min(vec1),min(vec2)]);
                        m2 = mean([max(vec1),max(vec2)]);
                        ud_mm_exc = sqrt((m2-m1)^2);
                        %* min max exc gc (ml)
                        vec = [daty(rhsi_1:rhsi_2,s_i)];
                        m1 = min(vec);
                        m2 = max(vec);
                        ml_mm_gc_exc = sqrt((m2-m1)^2);
                        %* min max exc gc (ap)
                        vec = [datx(rhsi_1:rhsi_2,s_i)];
                        m1 = min(vec);
                        m2 = max(vec);
                        ap_mm_gc_exc = sqrt((m2-m1)^2);
                        %* min max exc gc (ud)
                        vec = [datz(rhsi_1:rhsi_2,s_i)];
                        m1 = min(vec);
                        m2 = max(vec);
                        ud_mm_gc_exc = sqrt((m2-m1)^2);
                        %-
                        tmp_step_struct(s_i).stride_n = inds_cond(s_i);
                        tmp_step_struct(s_i).rhs_1 = rhs_1;
                        tmp_step_struct(s_i).rhs_2 = rhs_2;
                        tmp_step_struct(s_i).lhs_1 = lhs_1;
                        tmp_step_struct(s_i).lto_1 = lto_1;
                        tmp_step_struct(s_i).rto_1 = rto_1;
                        tmp_step_struct(s_i).gait_cycle_dur = rhs_2 - rhs_1;
                        tmp_step_struct(s_i).stance_dur = rto_1 - rhs_1;
                        tmp_step_struct(s_i).swing_dur = rhs_2 - rto_1;
                        tmp_step_struct(s_i).double_sup_dur = rto_1 - lhs_1;
                        tmp_step_struct(s_i).single_sup_dur = lhs_1 - lto_1;
                        tmp_step_struct(s_i).step_dur = rhs_2 - lhs_1;
                        tmp_step_struct(s_i).ml_exc = ml_exc;
                        tmp_step_struct(s_i).ap_exc = ap_exc;
                        tmp_step_struct(s_i).ud_exc = ud_exc;
                        % fprintf('ml_exc_mm: %0.2g\n',ml_mm_exc);
                        tmp_step_struct(s_i).ml_exc_mm = ml_mm_exc;
                        tmp_step_struct(s_i).ap_exc_mm = ap_mm_exc;
                        tmp_step_struct(s_i).ud_exc_mm = ud_mm_exc;
                        %-
                        tmp_step_struct(s_i).ml_exc_mm_gc = ml_mm_gc_exc;
                        tmp_step_struct(s_i).ap_exc_mm_gc = ap_mm_gc_exc;
                        tmp_step_struct(s_i).ud_exc_mm_gc = ud_mm_gc_exc;
                        %-
                        do_log_s(s_i) = true;
                    end            
                end
                fprintf('\n');
        
                %## ADD EEG MEASURES
                cs_inds = cellfun(@(x) strcmp(x,conds{cond_i}),{cond_struct.cond}); 
                tmp_cs = cond_struct(cs_inds);
                for iii = 1:length(tmp_cs)
                    %- struct
                    for i = 1:length(CLUSTER_PICS)
                        cl_i = CLUSTER_PICS(i);
                        comp_n = comp_arr(1,comp_arr(3,:) == cl_i);
                        if ~isempty(comp_n)
                            s_ii = tmp_cs(iii).indices;
                            stride_inds = [tmp_step_struct.stride_n];
                            [~,ss_inds] = intersect(stride_inds,s_ii);
                            %##
                            steps_struct(cnt).subj_char = EEG.subject;
                            steps_struct(cnt).cond_char = conds{cond_i}; %speed_n_chars{cond_i}; %conds{cond_i};
                            steps_struct(cnt).group_char = EEG.group;
                            steps_struct(cnt).speed_char = speed_n_chars{cond_i}; %conds{cond_i};
                            steps_struct(cnt).speed_n = double(string(speed_n_chars{cond_i})); %conds{cond_i};
                            steps_struct(cnt).model_char = 'speed';
                            steps_struct(cnt).model_n = des_i;
                            steps_struct(cnt).cluster_n = cl_i;
                            steps_struct(cnt).comp_n = comp_n; % new comp number
                            steps_struct(cnt).comp_n_old = comp_arr(2,comp_arr(3,:) == cl_i); % old comp number
                            tmpfe = fname_ext(~cellfun(@isempty,fname_ext));
                            steps_struct(cnt).epoch_type = strjoin(tmpfe,'_');
                            
                            %## MEAN MEAS.
                            steps_struct(cnt).mu_stride_n = mean([tmp_step_struct(ss_inds).stride_n]);
                            steps_struct(cnt).mu_rhs_1 = mean([tmp_step_struct(ss_inds).rhs_1]);
                            steps_struct(cnt).mu_rhs_2 = mean([tmp_step_struct(ss_inds).rhs_2]);
                            steps_struct(cnt).mu_lhs_1 = mean([tmp_step_struct(ss_inds).lhs_1]);
                            steps_struct(cnt).mu_lto_1 = mean([tmp_step_struct(ss_inds).lto_1]);
                            steps_struct(cnt).mu_rto_1 = mean([tmp_step_struct(ss_inds).rto_1]);
                            steps_struct(cnt).mu_gait_cycle_dur = mean([tmp_step_struct(ss_inds).gait_cycle_dur]);
                            steps_struct(cnt).mu_stance_dur = mean([tmp_step_struct(ss_inds).stance_dur]);
                            steps_struct(cnt).mu_swing_dur = mean([tmp_step_struct(ss_inds).swing_dur]);
                            steps_struct(cnt).mu_double_sup_dur = mean([tmp_step_struct(ss_inds).double_sup_dur]);
                            steps_struct(cnt).mu_single_sup_dur = mean([tmp_step_struct(ss_inds).single_sup_dur]);
                            steps_struct(cnt).mu_step_dur = mean([tmp_step_struct(ss_inds).step_dur]);
                            %-
                            steps_struct(cnt).mu_ml_exc_mm_gc = mean([tmp_step_struct(ss_inds).ml_exc_mm_gc]);
                            steps_struct(cnt).mu_ap_exc_mm_gc = mean([tmp_step_struct(ss_inds).ap_exc_mm_gc]);
                            steps_struct(cnt).mu_ud_exc_mm_gc = mean([tmp_step_struct(ss_inds).ud_exc_mm_gc]);
                            %-                    
                            steps_struct(cnt).mu_avg_theta = mean(tmp_psd(tband,tmp_cs(iii).ind,comp_n),1);
                            steps_struct(cnt).mu_avg_alpha = mean(tmp_psd(aband,tmp_cs(iii).ind,comp_n),1);
                            steps_struct(cnt).mu_avg_beta = mean(tmp_psd(bband,tmp_cs(iii).ind,comp_n),1);
                            steps_struct(cnt).med_avg_theta = median(tmp_psd(tband,tmp_cs(iii).ind,comp_n),1);
                            steps_struct(cnt).med_avg_alpha = median(tmp_psd(aband,tmp_cs(iii).ind,comp_n),1);
                            steps_struct(cnt).med_avg_beta = median(tmp_psd(bband,tmp_cs(iii).ind,comp_n),1);
                            
                            %## STANDARD DEV MEAS.
                            steps_struct(cnt).std_stride_n = std([tmp_step_struct(ss_inds).stride_n]);
                            steps_struct(cnt).std_rhs_1 = std([tmp_step_struct(ss_inds).rhs_1]);
                            steps_struct(cnt).std_rhs_2 = std([tmp_step_struct(ss_inds).rhs_2]);
                            steps_struct(cnt).std_lhs_1 = std([tmp_step_struct(ss_inds).lhs_1]);
                            steps_struct(cnt).std_lto_1 = std([tmp_step_struct(ss_inds).lto_1]);
                            steps_struct(cnt).std_rto_1 = std([tmp_step_struct(ss_inds).rto_1]);
                            steps_struct(cnt).std_gait_cycle_dur = std([tmp_step_struct(ss_inds).gait_cycle_dur]);
                            steps_struct(cnt).std_stance_dur = std([tmp_step_struct(ss_inds).stance_dur]);
                            steps_struct(cnt).std_swing_dur = std([tmp_step_struct(ss_inds).swing_dur]);
                            steps_struct(cnt).std_double_sup_dur = std([tmp_step_struct(ss_inds).double_sup_dur]);
                            steps_struct(cnt).std_single_sup_dur = std([tmp_step_struct(ss_inds).single_sup_dur]);
                            steps_struct(cnt).std_step_dur = std([tmp_step_struct(ss_inds).step_dur]);
                            %-
                            steps_struct(cnt).std_ml_exc_mm_gc = std([tmp_step_struct(ss_inds).ml_exc_mm_gc]);
                            steps_struct(cnt).std_ap_exc_mm_gc = std([tmp_step_struct(ss_inds).ap_exc_mm_gc]);
                            steps_struct(cnt).std_ud_exc_mm_gc = std([tmp_step_struct(ss_inds).ud_exc_mm_gc]);
                            %--
                            steps_struct(cnt).std_avg_theta = mean(tmp_psd_std(tband,tmp_cs(iii).ind,comp_n),1);
                            steps_struct(cnt).std_avg_alpha = mean(tmp_psd_std(aband,tmp_cs(iii).ind,comp_n),1);
                            steps_struct(cnt).std_avg_beta = mean(tmp_psd_std(bband,tmp_cs(iii).ind,comp_n),1);
        
                            %## COV MEAS.
                            steps_struct(cnt).cov_stride_n = std([tmp_step_struct(ss_inds).stride_n])/mean([tmp_step_struct(ss_inds).stride_n]);
                            steps_struct(cnt).cov_rhs_1 = std([tmp_step_struct(ss_inds).rhs_1])/mean([tmp_step_struct(ss_inds).rhs_1]);
                            steps_struct(cnt).cov_rhs_2 = std([tmp_step_struct(ss_inds).rhs_2])/mean([tmp_step_struct(ss_inds).rhs_2]);
                            steps_struct(cnt).cov_lhs_1 = std([tmp_step_struct(ss_inds).lhs_1])/mean([tmp_step_struct(ss_inds).lhs_1]);
                            steps_struct(cnt).cov_lto_1 = std([tmp_step_struct(ss_inds).lto_1])/mean([tmp_step_struct(ss_inds).lto_1]);
                            steps_struct(cnt).cov_rto_1 = std([tmp_step_struct(ss_inds).rto_1])/mean([tmp_step_struct(ss_inds).rto_1]);
                            steps_struct(cnt).cov_gait_cycle_dur = std([tmp_step_struct(ss_inds).gait_cycle_dur])/mean([tmp_step_struct(ss_inds).gait_cycle_dur]);
                            steps_struct(cnt).cov_stance_dur = std([tmp_step_struct(ss_inds).stance_dur])/mean([tmp_step_struct(ss_inds).stance_dur]);
                            steps_struct(cnt).cov_swing_dur = std([tmp_step_struct(ss_inds).swing_dur])/mean([tmp_step_struct(ss_inds).swing_dur]);
                            steps_struct(cnt).cov_double_sup_dur = std([tmp_step_struct(ss_inds).double_sup_dur])/mean([tmp_step_struct(ss_inds).double_sup_dur]);
                            steps_struct(cnt).cov_single_sup_dur = std([tmp_step_struct(ss_inds).single_sup_dur])/mean([tmp_step_struct(ss_inds).single_sup_dur]);
                            steps_struct(cnt).cov_step_dur = std([tmp_step_struct(ss_inds).step_dur])/mean([tmp_step_struct(ss_inds).step_dur]);
                            %--
                            steps_struct(cnt).cov_ml_exc_mm_gc = std([tmp_step_struct(ss_inds).ml_exc_mm_gc])/mean([tmp_step_struct(ss_inds).ml_exc_mm_gc]);
                            steps_struct(cnt).cov_ap_exc_mm_gc = std([tmp_step_struct(ss_inds).ap_exc_mm_gc])/mean([tmp_step_struct(ss_inds).ap_exc_mm_gc]);
                            steps_struct(cnt).cov_ud_exc_mm_gc = std([tmp_step_struct(ss_inds).ud_exc_mm_gc])/mean([tmp_step_struct(ss_inds).ud_exc_mm_gc]);
                            %--
                            % steps_struct(cnt).cov_avg_theta = 100*(mean(tmp_psd_std(tband,tmp_cs(iii).ind,comp_n),1)/abs(mean(tmp_psd(tband,tmp_cs(iii).ind,comp_n),1)));
                            % steps_struct(cnt).cov_avg_alpha = 100*(mean(tmp_psd_std(aband,tmp_cs(iii).ind,comp_n),1)/abs(mean(tmp_psd(aband,tmp_cs(iii).ind,comp_n),1)));
                            % steps_struct(cnt).cov_avg_beta = 100*(mean(tmp_psd_std(bband,tmp_cs(iii).ind,comp_n),1)/abs(mean(tmp_psd(bband,tmp_cs(iii).ind,comp_n),1)));
                            
                            steps_struct(cnt).cov_avg_theta = 100*(median(tmp_psd_std(tband,tmp_cs(iii).ind,comp_n),1)/abs(median(tmp_psd(tband,tmp_cs(iii).ind,comp_n),1)));
                            steps_struct(cnt).cov_avg_alpha = 100*(median(tmp_psd_std(aband,tmp_cs(iii).ind,comp_n),1)/abs(median(tmp_psd(aband,tmp_cs(iii).ind,comp_n),1)));
                            steps_struct(cnt).cov_avg_beta = 100*(median(tmp_psd_std(bband,tmp_cs(iii).ind,comp_n),1)/abs(median(tmp_psd(bband,tmp_cs(iii).ind,comp_n),1)));                            
                            % --
                            steps_struct(cnt).cov_i_avg_theta = squeeze(tmp_covt(1,tmp_cs(iii).ind,comp_n));
                            steps_struct(cnt).cov_i_avg_alpha = squeeze(tmp_cova(1,tmp_cs(iii).ind,comp_n));
                            steps_struct(cnt).cov_i_avg_beta = squeeze(tmp_covb(1,tmp_cs(iii).ind,comp_n));
                            %--
                            steps_struct(cnt).qcv_i_avg_theta = squeeze(tmp_qcvt(1,tmp_cs(iii).ind,comp_n));
                            steps_struct(cnt).qcv_i_avg_alpha = squeeze(tmp_qcva(1,tmp_cs(iii).ind,comp_n));
                            steps_struct(cnt).qcv_i_avg_beta = squeeze(tmp_qcvb(1,tmp_cs(iii).ind,comp_n));
                            %-- sanity chk
                            % fprintf('%s) theta COV1: mu = %0.2f, std = %0.2f\n',conds{cond_i},mean([steps_struct.cov_avg_theta]),std([steps_struct.cov_avg_theta]))
                            % fprintf('%s) alpha COV1: mu = %0.2f, std = %0.2f\n',conds{cond_i},mean([steps_struct.cov_avg_alpha]),std([steps_struct.cov_avg_alpha]))
                            % fprintf('%s) beta COV1: mu = %0.2f, std = %0.2f\n',conds{cond_i},mean([steps_struct.cov_avg_beta]),std([steps_struct.cov_avg_beta]))
                            % fprintf('%s) theta COVi: mu = %0.2f, std = %0.2f\n',conds{cond_i},mean([steps_struct.cov_i_avg_theta]),std([steps_struct.cov_i_avg_theta]))
                            % fprintf('%s) alpha COVi: mu = %0.2f, std = %0.2f\n',conds{cond_i},mean([steps_struct.cov_i_avg_alpha]),std([steps_struct.cov_i_avg_alpha]))
                            % fprintf('%s) beta COVi: mu = %0.2f, std = %0.2f\n',conds{cond_i},mean([steps_struct.cov_i_avg_beta]),std([steps_struct.cov_i_avg_beta]))
                            % fprintf('%s) theta QCV: mu = %0.2f, std = %0.2f\n',conds{cond_i},mean([steps_struct.qcv_i_avg_theta]),std([steps_struct.qcv_i_avg_theta]))
                            % fprintf('%s) alpha QCV: mu = %0.2f, std = %0.2f\n',conds{cond_i},mean([steps_struct.qcv_i_avg_alpha]),std([steps_struct.qcv_i_avg_alpha]))
                            % fprintf('%s) beta QCV: mu = %0.2f, std = %0.2f\n',conds{cond_i},mean([steps_struct.qcv_i_avg_beta]),std([steps_struct.qcv_i_avg_beta]))
                            
                            %## AP PARAMS
                            steps_struct(cnt).mu_ap_exponent = squeeze(mean(tmp_ap_mu(2,tmp_cs(iii).ind,comp_n),2));
                            steps_struct(cnt).std_ap_exponent = squeeze(mean(tmp_ap_std(2,tmp_cs(iii).ind,comp_n),2));
                            steps_struct(cnt).mu_ap_offset = squeeze(mean(tmp_ap_mu(1,tmp_cs(iii).ind,comp_n),2));
                            steps_struct(cnt).std_ap_offset = squeeze(mean(tmp_ap_std(1,tmp_cs(iii).ind,comp_n),2));
                            steps_struct(cnt).cov_ap_exponent = squeeze(mean(tmp_ap_cov(2,tmp_cs(iii).ind,comp_n),2));
                            steps_struct(cnt).cov_ap_offset = squeeze(mean(tmp_ap_cov(1,tmp_cs(iii).ind,comp_n),2));
                            %--
                            cnt = cnt + 1;
                        end
                    end
                end
                steps_struct = steps_struct(~cellfun(@isempty,{steps_struct.mu_avg_theta}));
                steps_struct = steps_struct(~cellfun(@isempty,{steps_struct.mu_ml_exc_mm_gc}));
                steps_struct = steps_struct(~cellfun(@isempty,{steps_struct.comp_n}));
                
            end
    
            %## FURTHER REMOVE BAD STRIDES
            steps_struct = steps_struct(~cellfun(@(x) x==0,{steps_struct.mu_avg_theta}));
            steps_struct = steps_struct(~cellfun(@(x) x==0,{steps_struct.mu_ml_exc_mm_gc}));
            %- and empty columns
            fns = fieldnames(steps_struct);
            for ffn = 1:length(fns)        
                %-
                chk1 = cellfun(@isempty,{steps_struct.(fns{ffn})});
                try
                    chk2 = cellfun(@(x) x==0,{steps_struct.(fns{ffn})});
                catch e
                    chk2 = false;
                    if ~strcmp(e.identifier,'MATLAB:cellfun:NotAScalarOutput')
                        fprintf('%s\n\n',getReport(e));
                    end
                end
                if all(chk1)            
                    steps_struct = rmfield(steps_struct,fns{ffn});
                end
                if all(chk2)
                    steps_struct = rmfield(steps_struct,fns{ffn});
                end
            end
            
            %-- sanity chk
            fprintf('%s) theta COV1: mu = %0.2f, std = %0.2f\n',EEG.subject,mean([steps_struct.cov_avg_theta]),std([steps_struct.cov_avg_theta]))
            fprintf('%s) alpha COV1: mu = %0.2f, std = %0.2f\n',EEG.subject,mean([steps_struct.cov_avg_alpha]),std([steps_struct.cov_avg_alpha]))
            fprintf('%s) beta COV1: mu = %0.2f, std = %0.2f\n',EEG.subject,mean([steps_struct.cov_avg_beta]),std([steps_struct.cov_avg_beta]))
            fprintf('%s) theta COVi: mu = %0.2f, std = %0.2f\n',EEG.subject,mean([steps_struct.cov_i_avg_theta]),std([steps_struct.cov_i_avg_theta]))
            fprintf('%s) alpha COVi: mu = %0.2f, std = %0.2f\n',EEG.subject,mean([steps_struct.cov_i_avg_alpha]),std([steps_struct.cov_i_avg_alpha]))
            fprintf('%s) beta COVi: mu = %0.2f, std = %0.2f\n',EEG.subject,mean([steps_struct.cov_i_avg_beta]),std([steps_struct.cov_i_avg_beta]))
            fprintf('%s) theta QCV: mu = %0.2f, std = %0.2f\n',EEG.subject,mean([steps_struct.qcv_i_avg_theta]),std([steps_struct.qcv_i_avg_theta]))
            fprintf('%s) alpha QCV: mu = %0.2f, std = %0.2f\n',EEG.subject,mean([steps_struct.qcv_i_avg_alpha]),std([steps_struct.qcv_i_avg_alpha]))
            fprintf('%s) beta QCV: mu = %0.2f, std = %0.2f\n',EEG.subject,mean([steps_struct.qcv_i_avg_beta]),std([steps_struct.qcv_i_avg_beta]))
            %## SAVE
            steps_struct = struct2table(steps_struct);
            fname_ext = fname_ext(~cellfun(@isempty,fname_ext));
            par_save(steps_struct,EEG.filepath,sprintf('kin_eeg_psd_meas_allcond_%s.mat',strjoin(fname_ext,'_')));
            fname_ext = fname_ext(~cellfun(@isempty,fname_ext));
    
            tmp_fname_exts{nn} = strjoin(fname_ext,'_');
            %- turn off when doing parallel processing.
            EEG.data = 'eeg_imu_epochs.fdt';
            % tmp_alleeg{subj_i} = EEG;
        end
    catch e
        disp(getReport(e));
    end
    fname_ext_store{subj_i,1} = tmp_fname_exts;
end
%% (LOAD & CONGREGATE PSD DATA) ===================================== %%
f_range = [3, 40];
%## FNAMES
fnames = {SBS_STUDY.datasetinfo.filename};
fpaths = {SBS_STUDY.datasetinfo.filepath};
subj_chars = {SBS_STUDY.datasetinfo.subject};
%-- loop
% fname_ext_store = {''}
for nn = 1:length(STRIDES_AVG)
    fextr = fname_ext_store{1}{nn};
    %## GET FREQ DIMS
    dat = par_load(fpaths{1},sprintf('psd_output_%s.mat',fextr));
    psd_mean = dat.psd_mean;
    %--
    psd_dat_out = nan(length(subj_chars),size(psd_mean,1),250,length(CLUSTER_PICS));
    psd_std_dat_out = nan(length(subj_chars),size(psd_mean,1),250,length(CLUSTER_PICS));
    cond_dat_out = cell(length(subj_chars),250,length(CLUSTER_PICS));
    subj_dat_out = nan(length(subj_chars),250,length(CLUSTER_PICS));
    group_dat_out = cell(length(subj_chars),250,length(CLUSTER_PICS));
    % conds_extract = {'0p25','0p5','0p75','1p0'};    
    %##
    for subj_i = 1:length(subj_chars)
        %-- initiate params
        tmp_sbs_study = SBS_STUDY;
        tmp_cl_study = CL_STUDY;
        tmp_subjs = {tmp_cl_study.datasetinfo.subject};
        try
            %## LOAD SET FIlE
            EEG = load([fpaths{subj_i} filesep fnames{subj_i}],'-mat');
            
            %## LOAD PSD DATA
            dat = par_load(fpaths{subj_i},sprintf('psd_output_%s.mat',fextr));
            %--
            psd_mean = dat.psd_mean; %[frequency, epoch/splice, channel/component]
            psd_std = dat.psd_std;
            fooof_freqs = dat.freqs;
            cond_struct = dat.cond_struct;
            fprintf('conds: %s\n',strjoin(unique({cond_struct.cond}),','));
            comp_arr = dat.indx_to_comp;
            psd_mean = psd_mean(:,[cond_struct.ind],:);
            psd_std = psd_std(:,[cond_struct.ind],:);
        
            %## GET CLUSTER DATA
            fprintf('Getting cluster information...\n');
            [~,subs] = intersect(comp_arr(3,:),main_cl_inds);
            tmp_arr = comp_arr(:,subs);
            
            %## EXTRACT DATA TO CLUSTERED ARRAYS
            fprintf('Extract data and assigning to big array...\n');    
            psd_dat_out(subj_i,:,1:size(psd_mean,2),tmp_arr(3,:)) = psd_mean(:,:,tmp_arr(1,:));
            psd_std_dat_out(subj_i,:,1:size(psd_std,2),tmp_arr(3,:)) = psd_std(:,:,tmp_arr(1,:));
            subj_dat_out(subj_i,1:size(psd_mean,2),tmp_arr(3,:)) = repmat(subj_i,[1,size(psd_mean,2),length(tmp_arr)]);
            tmp = {cond_struct.cond};
            cond_dat_out(subj_i,1:length(cond_struct),tmp_arr(3,:)) = repmat(tmp',[1,length(tmp_arr)]);
            group_dat_out(subj_i,1:length(cond_struct),tmp_arr(3,:)) = repmat({EEG.group},[1,size(psd_mean,2),length(tmp_arr)]);
            fprintf('psd len: %i\ncond len: %i\n',size(psd_mean,2),length(cond_struct));
            % psd_dat_out(subj_i) = psd_mean(:,:,tmp_arr(1,:));
        catch e
            fprintf('Couldn''t load %s...\n%s\n',EEG.subject,getReport(e));
        end
    end
    inds = ~all(cellfun(@isempty,group_dat_out),[4,3,2]);
    psd_dat_out = psd_dat_out(inds,:,:,:);
    psd_std_dat_out = psd_std_dat_out(inds,:,:,:);
    cond_dat_out = cond_dat_out(inds,:,:,:);
    subj_dat_out = subj_dat_out(inds,:,:,:);
    group_dat_out = group_dat_out(inds,:,:,:);
    dat_out_struct = struct('psd_dat',psd_dat_out, ...
        'psd_std_dat',psd_std_dat_out, ...
        'cond_dat',{cond_dat_out}, ...
        'subj_dat',subj_dat_out, ...
        'group_dat',{group_dat_out});
    par_save(dat_out_struct,[cluster_k_dir filesep 'kin_eeg_step_to_step' filesep sprintf('raw_psd_dat_%s.mat',fextr)]);
    dat_out_struct = struct.empty;
    %## LOAD
    %{
    dat_out_struct = par_load([cluster_k_dir filesep 'kin_eeg_step_to_step' filesep sprintf('raw_psd_dat_%s.mat',fextr)]);
    psd_dat_out = dat_out_struct.psd_dat;
    cond_dat_out = dat_out_struct.cond_dat;
    subj_dat_out = dat_out_struct.subj_dat;
    %}
end
%% (CREATE BIG TABLE ACROSS SPEEDS) ==================================== %%
% data_store = table.empty;
groups = {'H1000','H2000','H3000'};
group_name = {'younger_adults','older_high_function','older_low_function'};
% fname_ext_store = {{'perstridefb_nfslidingb36'}};
for nn = 1:length(STRIDES_AVG)
    data_store = cell(1,length(SBS_STUDY.datasetinfo));
    fextr = fname_ext_store{1}{nn};
    for i = 1:length(SBS_STUDY.datasetinfo)
        try
            tmp = par_load(SBS_STUDY.datasetinfo(i).filepath,sprintf('kin_eeg_psd_meas_allcond_%s.mat',fextr));
            vals = regexp(tmp.subj_char{1},'[nN]?[hH](\d)\d*','tokens');
            vals = double(string(vals{1}{1}));
            tt = repmat(vals,[1,size(tmp,1)]);
            tmp.group_n = tt';
            tt = repmat(groups(vals),[1,size(tmp,1)]);
            [tmp.group_char] = tt';
            tt = repmat(group_name(vals),[1,size(tmp,1)]);
            [tmp.group_name] = tt';
            data_store{i} = tmp;
        catch e
            fprintf('Couldn''t load %s...\n%s\n',SBS_STUDY.datasetinfo(i).subject,getReport(e));
        end
    end
    data_store = util_resolve_table(data_store);
    writetable(data_store,[save_dir filesep sprintf('sbs_eeg_psd_kin_allcond_%s.xlsx',fextr)]);
    par_save(data_store,save_dir,sprintf('sbs_eeg_psd_kin_allcond_%s.mat',fextr));
end

