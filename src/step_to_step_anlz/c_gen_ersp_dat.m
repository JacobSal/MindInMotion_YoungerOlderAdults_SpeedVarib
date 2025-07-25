%   Project Title: MIM OA & YA SPEED & KINETICS ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/step_to_step_anlz/run_sts_c_gen_itc_dat.sh

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
FORCE_RECALC_ERSP = true;
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
NEWTIMEF_STRUCT = struct(...
    'timewindow',[], ...
    'timewarp',[], ...
    'timewarpms',[], ...
    'ntimesout',[], ...
    'pointrange',[], ...
    'baseline',nan(), ...
    'timelimits',[], ...
    'freqs',[3,250], ...
    'freqscale','log', ...
    'do_plot_ersp',false, ...
    'do_plot_itc',false, ...
    'do_plot_phase',false, ...
    'cycles',[3,0.8], ...
    'padratio',1, ...
    'nfreqs',200, ...
    'alpha',nan(), ...
    'srate',[] ...
    );
%% (PATHS)
%- datset name
DATA_SET = 'MIM_dataset';
%- study name
STUDY_DNAME = '02202025_mim_yaoa_powpow0p3_crit_speed';
STUDY_FNAME = 'kin_eeg_epoch_study';
STUDY_EPOCH_FNAME = 'epoch_study';
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
%## LOAD EPOCH STUDY
% need this for the gait timewarp information
if ~ispc
    tmp = load('-mat',[studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep sprintf('%s_UNIX.study',STUDY_EPOCH_FNAME)]);
    E_STUDY = tmp.STUDY;
else
    tmp = load('-mat',[studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep sprintf('%s.study',STUDY_EPOCH_FNAME)]);
    E_STUDY = tmp.STUDY;
end
%% (STUDY DESIGN) ====================================================== %%
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

%## TIME WARP
SRATE = 500;
TIMEWARP_NTIMES = floor(SRATE/pi); % conservative nyquist frequency. making this too big can cause overlap between gait cyles
averaged_warpto_events = E_STUDY.etc.a_epoch_process.avg_warpto_events;
ERSP_CROP_TIMES = [averaged_warpto_events(1), averaged_warpto_events(end)+1];
disp(['Grand average (across all subj) warp to: ',num2str(averaged_warpto_events)]);
%--
CL_STUDY.cluster = cl_struct;
%% ERSP
SRATE = 500;
WINDOW_LEN = 500; %ms (for cycle=[3,0.8], the largest window len is usually around 500
WINDOW_SLIDE = 25; %ms
spin = SRATE*(sp/1000);
nt = length(4500)/spin;
lf = floor(pi()*(SRATE/(WINDOW_LEN)));   
MAX_CL_NUM = 14;
MOD_CHARS = {{'flat','low','med','high'},{'0p25','0p5','0p75','1p0'}};
TF_DAT_STRUCT = struct(...
        'subj_char',{''}, ...
        'group_char',{''}, ...
        'cond_char',{''}, ...
        'mod_n',[], ...
        'comp_n',[], ...
        'cluster_n',[], ...
        'itc_dat',[], ...
        'itc_times',[], ...
        'itc_freqs',[] ...
        );
%-- newtimef func params
NTF_STRUCT = struct(...
    'timewindow',[], ...
    'timewarp',[], ...
    'timewarpms',[], ...
    'ntimesout',floor(SRATE/2), ...
    'pointrange',[], ...
    'baseline',nan(), ...
    'timelimits',[], ...
    'freqs',[4,SRATE/2], ...
    'itctype','phasecoher', ...  % this is the measure AS used
    'freqscale','log', ...
    'do_plot_ersp',false, ...
    'do_plot_itc',false, ...
    'do_plot_phase',false, ...
    'do_single_dat',true, ...
    'do_timewarp',true, ...
    'cycles',[3,0.8], ...
    'padratio',2^2, ...
    'nfreqs',floor(SRATE/pi()), ...
    'alpha',nan(), ...
    'srate',SRATE ...
);
DAT_FEXT = 'allersp';

%## GET ITERS
%-- iters
iters = zeros(length(CL_STUDY.datasetinfo)*30,3);
%-- data store
sdinfo = {CL_STUDY.datasetinfo.subject};
scinfo = CL_STUDY.cluster;
%-- iter assign
cnt = 1;
for s_i = 1:length(CL_STUDY.datasetinfo)
    sc = CL_STUDY.datasetinfo(s_i).subject;
    tmp = scinfo(1);
    tc = tmp.comps;
    ts = tmp.sets;
    tcn = tc(ts == s_i);
    for k_i = 1:length(tcn)
        cc = false;
        cl_i = 2;
        while ~any(cc) && cl_i < length(scinfo)
            cl_i = cl_i + 1;
            ind = find(strcmp(sc,sdinfo));
            ind = scinfo(cl_i).sets == ind;
            cc = scinfo(cl_i).comps(ind) == k_i;                
        end
        iters(cnt,1) = cl_i;
        iters(cnt,2) = k_i;
        iters(cnt,3) = s_i;
        cnt = cnt + 1;
    end
end
inds = all(iters == 0,2) | iters(:,1) > MAX_CL_NUM;
iters = iters(~inds,:);
% parfor s_i = 1:length(SBS_STUDY.datasetinfo)
parfor i = 1:size(iters,1)
    %## PARAM COPIES
    tmpi = iters;
    cl_i = tmpi(i,1);
    k_i = tmpi(i,2);
    s_i = tmpi(i,3);
    tmp_cl_study = CL_STUDY;
    tmp_ersp_params = ERSP_PARAMS;
    tmp_ntf_struct = NTF_STRUCT;
    tmp_tfd_struct = TF_DAT_STRUCT;
    %--
    tt = tic();
    fprintf('%s) Running ITC calculation:\n',tmp_cl_study.datasetinfo(s_i).subject);

    %-- storage
    try
        EEG = pop_loadset('filepath',tmp_cl_study.datasetinfo(s_i).filepath, ...
            'filename',tmp_cl_study.datasetinfo(s_i).filename);
        fprintf('Running Subject %s\n',EEG.subject);
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        end
        EEG.icaact = EEG.icaact(k_i,:,:);
        ics_calc = EEG.etc.urreject.ic_keep;

        %## PARAMETERS IN
        tlimits = [EEG(1).xmin EEG(1).xmax]*1000;
        pointrange1 = round(max((tlimits(1)/1000-EEG(1).xmin)*EEG(1).srate, 1));
        pointrange2 = round(min(((tlimits(2)+1000/EEG(1).srate)/1000-EEG(1).xmin)*EEG(1).srate, EEG(1).pnts));
        pointrange = (pointrange1:pointrange2);
        
        %-- newtimef func params
        tmp_ntf_struct.timewarp =  EEG.timewarp.latencies;
        tmp_ntf_struct.timewarpms = averaged_warpto_events;
        tmp_ntf_struct.timelimits = [EEG(1).xmin,EEG(1).xmax]*1000;
        tmp_ntf_struct.pointrange = pointrange;

        %-- get conditiion info
        ind = find(strcmp(EEG.subject,{tmp_cl_study.datasetinfo.subject}));
        trialinfo = std_combtrialinfo(tmp_cl_study.datasetinfo,ind);      
        %--
        fprintf('%s) Running ITC calculation:\n',EEG.subject);
        tt = tic();
        %--           
        EEG.icaact = EEG.icaact(k_i,:,:);
        %-- run newtimef
        [~,~,~,ttimes,lfreqs,~,~,all_e_dat] = ...
            eeglab_newtimef_iter(EEG,1, ...
            'NEWTIMEF_STRUCT',tmp_ntf_struct);
        %--
        %(04/29/2025) JS, newtimefitc seems to take advantage of the
        %precomputed time-freq measures from std_precompute, yet still
        %generates the same values as the newtimef "itc" output;
        %however, by default both use "phasecoher" as itctype.
        
        %## MI INDEX (Tort et al. 2010)
        tmpd = EEG.icaact(:,tmp_ntf_struct.pointrange,:);
        sigd  = reshape(tmpd,[1,length(tmp_ntf_struct.pointrange)*size(tmpd,3)]);
        [comod,p_vec,f_vec] = mod_index_calc(sigd);            
       
        %## STORE IN STRUCT  
        twmed = mean(EEG.timewarp.latencies(inds_cond,:),1);
        %--
        tmp_tfd_struct.subj_char = EEG.subject;
        tmp_tfd_struct.group_char = EEG.group;
        tmp_tfd_struct.subj_char = EEG.subject;
        tmp_tfd_struct.comp_n = k_i;
        tmp_tfd_struct.cluster_n = cl_i;
        tmp_tfd_struct.all_e_dat = single(all_e_dat);
        tmp_tfd_struct.twc_times = twmed;
        tmp_tfd_struct.twg_times = averaged_warpto_events; 
        %-- store
        % itc_dato{i} = tmp_itc_dato;

        % fprintf('%s) ERSP generation took %0.2f min\n',EEG.subject,toc(tt)/60);
        %## RESHAPE DATA

        %## ASSIGN DATA & ?SAVE?
        timef_data_struct = struct(...
            'itc_dat_struct',tmp_tfd_struct, ...
            'datatype','TIMEF', ...
            'datafiles',[EEG.filepath filesep EEG.filename], ...
            'datatrials',ics_calc, ...
            'parameters',tmp_ntf_struct, ...
            'trialinfo',trialinfo);
        par_save(timef_data_struct,[EEG.filepath filesep sprintf('%s_cl%i_custom_%s.mat',tmp_cl_study.datasetinfo(s_i).subject,cl_i,DAT_FEXT)]);
    catch e
        fprintf('%s) Something happened during ERSP generation...\n%s\n',tmp_cl_study.datasetinfo(s_i).subject,getReport(e));
    end
end

