%   Project Title: MIM OA & YA SPEED & KINETICS ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_speed_kin/step_to_step_anlz/run_dddd_eeg_ersp_imuls_sts_anl.sh

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
        STUDY_DIR = fileparts(SCRIPT_DIR); % change this if in sub folder
        SRC_DIR = fileparts(fileparts(STUDY_DIR));
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        STUDY_DIR = getenv('STUDY_DIR');
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
    STUDY_DIR = fileparts(SCRIPT_DIR); % change this if in sub folder
    SRC_DIR = fileparts(fileparts(STUDY_DIR));
end
%## Add Study, Src, & Script Paths
addpath(SCRIPT_DIR)
addpath(SRC_DIR);
addpath(STUDY_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
% [SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca_speed');
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
%- study name
STUDY_DNAME = '10172024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej_speed';
STUDY_FNAME = 'kin_eeg_epoch_study';
ANALYSIS_DNAME = 'kin_eeg_ersp_step_to_step';
%- statistics & conditions
speed_trials = {'0p25','0p5','0p75','1p0'};
terrain_trials = {'flat','low','med','high'};
COND_DESIGNS = {speed_trials,terrain_trials};
%* ERSP PARAMS
ERSP_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','cluster',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
SPEC_PARAMS = struct('freqrange',[1,200],...
    'subject','',...
    'specmode','psd',...
    'freqfac',4,...
    'logtrials','on',...
    'comps','all',...
    'plot_freqrange',[4,60],...
    'plot_ylim',[-35,-8],...
    'subtractsubjectmean','on',...
    'plotmode','normal');
ERSP_PARAMS = struct('subbaseline','off',...
    'timerange',[],...
    'ersplim',[-2,2],...
    'freqfac',4,...
    'cycles',[3,0.8],...
    'freqrange',[1,200]);
%% (PATHS)
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
%- load cluster
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb10'];
% cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb3'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
cluster_k_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
%% ===================================================================== %%
%## LOAD ALLEEG & STUDY
%{
if ~ispc
    [STUDY,ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME '_UNIX.study'],'filepath',[studies_fpath filesep sprintf('%s',STUDY_DNAME)]);
else
    [STUDY,ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME '.study'],'filepath',[studies_fpath filesep sprintf('%s',STUDY_DNAME)]);
end
%% CALCULATE GRANDAVERAGE WARPTOs
for subj_i = 1:length(ALLEEG)
    %- assign percondition timewarping
    ALLEEG(subj_i).timewarp.warpto = nanmedian(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto));
end
allWarpTo = nan(length(ALLEEG),size(ALLEEG(1).timewarp.warpto,2));
for subj_i = 1:length(ALLEEG)
    allWarpTo(subj_i,:) = ALLEEG(subj_i).timewarp.warpto; 
end
STUDY.etc.timewarp.avg_warpto_times = floor(nanmean(allWarpTo));
STUDY.etc.timewarp.ntimes = floor(ALLEEG(1).srate/pi);
STUDY.etc.timewarp.crop_times = [STUDY.etc.timewarp.avg_warpto_times(1), STUDY.etc.timewarp.avg_warpto_times(end)+1];
[STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                                        STUDY.filename,STUDY.filepath,...
                                        'RESAVE_DATASETS','off');
%}
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
    STUDY = tmp.STUDY;
else
    tmp = load('-mat',[studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep sprintf('%s.study',STUDY_FNAME)]);
    STUDY = tmp.STUDY;
end
%-
%-
cl_struct = par_load(cluster_k_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
[comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
CLUSTER_PICS = main_cl_inds;
%-
save_dir = [cluster_k_dir filesep ANALYSIS_DNAME];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% (STUDY DESIGN) ====================================================== %%
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
STUDY = pop_statparams(STUDY,args{:});
args = eeglab_struct2args(ERSP_PARAMS);
STUDY = pop_erspparams(STUDY,args{:});
STUDY.cache = [];
for des_i = 1:length(STUDY_DESI_PARAMS)
    [STUDY] = std_makedesign(STUDY,[],des_i,STUDY_DESI_PARAMS{des_i}{:});
end
%% (SUBJECT LOOP)
%## PARAMETERS
des_i = 2;
speed_n_chars = {'0.25','0.50','0.75','1.0'};
conds = STUDY.design(des_i).variable.value; %event.cond});
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
%## VALIDATION PLOT STRUCT
%-
freq_bound = [4,60];
%-
epoched_fPath = strsplit(CL_STUDY.datasetinfo(1).filepath,filesep);
icatimef_f = load('-mat',[strjoin(epoched_fPath,filesep) filesep sprintf('%s.icatimef',STUDY.datasetinfo(1).subject)]);
ersp_times = icatimef_f.times;
ersp_freqs = icatimef_f.freqs;
timewarpms = icatimef_f.parameters{find(strcmp(icatimef_f.parameters,'timewarpms'))+1}; %STUDY.etc.timewarp.avg_warpto_times;
% ntimes = STUDY.etc.timewarp.ntimes;
% crop_times = STUDY.etc.timewarp.crop_times;
time_bound = [timewarpms(1),timewarpms(end)];
freq_crop = find(ersp_freqs>=freq_bound(1) & ersp_freqs<=freq_bound(2));
time_crop = find(ersp_times>=time_bound(1) & ersp_times<=time_bound(2));
%-
cmap_terrain = linspecer(4);
custom_yellow = [254,223,0]/255;
cmap_terrain = [cmap_terrain(3,:);custom_yellow;cmap_terrain(4,:);cmap_terrain(2,:)];
cmap_speed = linspecer(4*3);
cmap_speed = [cmap_speed(1,:);cmap_speed(2,:);cmap_speed(3,:);cmap_speed(4,:)];
%-
PLOT_STRUCT = struct('figure_position_inch',[0.5,0.5,6.5,9],...
            'alltitles',{{}},...
            'xlabel','Gait Cycle Time (ms)',...
            'ylabel','Frequency (Hz)',...
            'xticklabel_times',timewarpms,...
            'xticklabel_chars',{{'RHS','LTO','LHS','RTO','RHS'}},...
            'xticklabel_angle',45,...
            'clim',[],...
            'font_size',8,...
            'font_name','Arial',...
            'freq_lims',freq_bound,...
            'time_lims',time_bound,...
            'stats_title','F Stat (p<0.05)',...
            'figure_title','',...
            'contourf_grain',ceil((500/pi())),...
            'group_titles',{{}},...
            'group_titles_shift_x',0.0,...
            'group_titles_shift_y',0.65,...
            'subplot_width',0.13,...
            'subplot_height',0.16,... %(02/17/2024) was 0.2
            'subplot_shift_x',0.035,...
            'subplot_shift_y',0.05,...
            'subplot_init_y',0.7,...
            'subplot_init_x',0.06,...
            'colorbar_shift_x',0.145,...
            'colorbar_shift_y',0,...
            'colorbar_label_shift_x',0,...
            'colorbar_label_shift_y',0.005,...
            'colorbar_label','\Delta Power (dB)',...
            'colorbar_fontsize',8,...
            'colorbar_fontweight','bold',...
            'display_bandmarks',true,...
            'bandmarks',{{'\theta','\alpha','\beta','\gamma'}},...
            'bandmarks_shift_y',[-0.42,-0.2,0.05,0.325],...
            'bandmarks_shift_x',-0.03,...
            'bandmarks_fontsize',8,...
            'bandmarks_fontweight','bold');
%##
tmp_alleeg = cell(length(STUDY.datasetinfo));
% data_store = [];
%%
% for subj_i = 1
parfor subj_i = 1:length(STUDY.datasetinfo)
    tmp_study = STUDY;
    tmp_cl_study = CL_STUDY;
    fooof_freqs = [];
    %## LOAD EEG DATA
    try
        EEG = pop_loadset('filepath',tmp_study.datasetinfo(subj_i).filepath,'filename',tmp_study.datasetinfo(subj_i).filename);
        subj_ind = strcmp(EEG.subject,{tmp_cl_study.datasetinfo.subject});
        fprintf('Running subject %s...\n',EEG.subject);
    catch e
        fprintf('%s',getReport(e));
        continue
    end
    if ~any(subj_ind)
        continue
    end
    %## LOAD PSD DATA
    % epoched_fPath = strsplit(EEG.filepath,filesep);
    epoched_fPath = strsplit(tmp_cl_study.datasetinfo(subj_ind).filepath,filesep);
    icatimef_f = [strjoin(epoched_fPath,filesep) filesep sprintf('%s.icatimef',EEG.subject)];
    %- load .icatimef load-in parameters
    fprintf('Loading Time-Frequency Data...\n');
    tmp = load(icatimef_f,'-mat');
    %- extract ersp data & info
    ersp_times = tmp.times;
    ersp_freqs = tmp.freqs;
    timewarpms = tmp.parameters{find(strcmp(tmp.parameters,'timewarpms'))+1}; %STUDY.etc.timewarp.avg_warpto_times;
    timewarp_trial = tmp.parameters{find(strcmp(tmp.parameters,'timewarp'))+1};
    time_bound = [timewarpms(1),timewarpms(end)+1]; %STUDY.etc.timewarp.crop_times; %[timewarpms(1),timewarpms(end)];
    freq_bound = [4,60];
    freq_cropd = ersp_freqs(freq_crop);
    time_cropd = ersp_times(time_crop);
    trialinfo = tmp.trialinfo;
    %- 
    fn = fieldnames(tmp);
    inds = find(contains(fieldnames(tmp),'comp'));
    test = tmp.(fn{inds(1)});
    %- 
    % eeg_ersp = zeros(size(test,1),size(test,2),size(test,3),length(inds),'single');
    eeg_ersp = zeros(length(freq_crop),length(time_crop),size(test,3),length(inds),'single');
    for i = 1:length(inds)
        tt = tmp.(fn{inds(i)});
        tt = tt(freq_crop,time_crop,:,:);
        eeg_ersp(:,:,:,i) = tt; % freq x time x epoch x chan
    end
    eeg_ersp = 20*log10(abs(eeg_ersp));
    %## SANITY CHECK
    %{
    i = 4;
    c_i = 4;
    inds = find(cellfun(@(x) any(strcmp(x,tmp_conds{c_i})),{trialinfo.cond}));
    tmp_ersp = {mean(squeeze(eeg_ersp(:,:,inds,i)),3)};
    tmp_ersp = {squeeze(eeg_ersp(:,:,inds(1),i))};
    fig = plot_txf_conds_tftopo(tmp_ersp,ersp_times,ersp_freqs,[],[],...
                'PLOT_STRUCT',PLOT_STRUCT);
    %}
    %- 
    fprintf('baselining ERSP data...\n');
    % tmp_conds = unique({trialinfo.trial_num_code});
    % tmp_conds = tmp_conds(~cellfun(@isempty,tmp_conds));
    tmp_conds = unique({trialinfo.cond});
    tmp_eegersp = zeros(size(eeg_ersp));
    %## LOOP THROUGH COMPONENTS
    for i = 1:size(tmp_eegersp,4)
        %## REMOVE MEAN OF EACH STEP
        % allersp = squeeze(eeg_ersp(:,:,:,i));
        % [allersp,~,~,~] = eeglab_baseln({allersp},ersp_times,ersp_freqs,...
        %     time_bound,freq_bound,...
        %     'DO_COMMON_BASE',false,...
        %     'DO_SUBJ_BASE',true);
        % allersp = allersp{1};

        %## REMOVE MEAN OF EACH CONDITION | TRIAL
        fprintf('Removing mean of each condition...\n');
        for c_i = 1:length(tmp_conds)
            %-
            inds = cellfun(@(x) any(strcmp(x,sprintf('%s',tmp_conds{c_i}))),{trialinfo.cond});
            %-
            % inds = cellfun(@(x) any(strcmp(x,sprintf('%s',tmp_conds{c_i}))),{trialinfo.trial_num_code});
            %(01/07/2025) JS, this could be potentially buggy given some
            %incosistencies in my trial determination algorithm.
            %(01/08/2025) JS, it is, in fact, buggy.
            if any(inds)
                allersp = squeeze(eeg_ersp(:,:,inds,i)); % gather data
                cond_mean = mean(squeeze(allersp),3); % mean across strides
                cond_mean = squeeze(mean(cond_mean,2)); % mean across time
                cond_mean = repmat(cond_mean,[1,size(eeg_ersp,2)]); % repeat across time
                cond_mean = repmat(cond_mean,[1,1,length(find(inds))]); % repeat across strides            
                allersp = allersp - cond_mean; % baseline
                tmp_eegersp(:,:,inds,i) = allersp;
            end
        end
        %## SANITY CHECK
        %{
        disp(c_i)
        c_i = 4;
        disp(tmp_conds{c_i})
        inds = cellfun(@(x) any(strcmp(x,tmp_conds{c_i})),{trialinfo.cond});
        tmpersp = tmp_eegersp(:,:,inds);
        tmpersp = mean(tmpersp,3);
        fig = plot_txf_conds_tftopo({tmpersp},time_cropd,freq_cropd,[],[],...
            'PLOT_STRUCT',PLOT_STRUCT);
        %}

        %## REMOVE MEAN OF WHOLE EXP. (TERRAIN | SPEED BASE)
        fprintf('Removing mean of each design...\n');
        %- whole design mean
        inds = cellfun(@(x) any(strcmp(x,{'0p25','0p5','0p75','1p0'})),{trialinfo.cond});
        %-
        % inds = cellfun(@(x) any(contains(x,{'0p25','0p5','0p75','1p0'})),{trialinfo.trial_num_code});
        %-
        tmp = tmp_eegersp(:,:,inds,i);
        all_mean = mean(squeeze(tmp(:,:,:)),3);
        all_mean = squeeze(mean(all_mean,2));
        all_mean = repmat(all_mean,[1,size(tmp,2)]);
        all_mean = repmat(all_mean,[1,1,size(tmp,3)]);
        allersp = tmp - all_mean;
        %- per stride baseline        
        [allersp,~,~,~] = eeglab_baseln({allersp},time_cropd,freq_cropd,...
            [time_cropd(1),time_cropd(end)],[freq_cropd(1),freq_cropd(end)],...
            'DO_COMMON_BASE',false,...
            'DO_SUBJ_BASE',true); 
        allersp = allersp{1};
        %## SANITY CHECK
        %{
        c_i = 4;
        disp(tmp_conds{c_i})
        inds = cellfun(@(x) any(strcmp(x,tmp_conds{c_i})),{trialinfo.cond});
        tmpersp = allersp(:,:,inds);
        tmpersp = mean(tmpersp,3);
        fig = plot_txf_conds_tftopo({tmpersp},time_cropd,freq_cropd,[],[],...
            'PLOT_STRUCT',PLOT_STRUCT);
        %}
        %- assign to array
        tmp_eegersp(:,:,inds,i) = allersp;
    end
    inds_keep = squeeze(~all(squeeze(tmp_eegersp(:,:,:,1)) == 0,[1,2]));
    eeg_ersp = tmp_eegersp(:,:,inds_keep,:);
    trialinfo = trialinfo(inds_keep);
    tmp_eegersp = double.empty;
    tmp_ersp = double.empty;
    allersp = double.empty;
    %-
    % eeg_ersp = tmp;
    % freqs = fooof_freqs;
    %## GET CLUSTER DATA
    fprintf('Getting cluster information...\n');
    tmp_subjs = {tmp_cl_study.datasetinfo.subject};
    comp_arr = zeros(3,length(CLUSTER_PICS)); %(1, comps; 2, old_comps; 3, cluster
    for i = 1:length(CLUSTER_PICS)
        cl_i = CLUSTER_PICS(i);
        ind = find(strcmp(EEG.subject,tmp_subjs));
        ind = tmp_study.cluster(cl_i).sets == ind;
        if any(ind)
            comp_arr(1,i) = tmp_study.cluster(cl_i).comps(ind);
            comp_arr(3,i) = cl_i;
            comp_arr(2,i) = EEG.etc.urreject.ic_keep(comp_arr(1,i));
        end
    end
    comp_arr = comp_arr(:,all(comp_arr,1));
    %##
    fprintf('Performing gait kinematic calculations...\n');
    tband = (ersp_freqs > theta_band_lims(1) & ersp_freqs < theta_band_lims(2));
    aband = (ersp_freqs > alpha_band_lims(1) & ersp_freqs < alpha_band_lims(2));
    bband = (ersp_freqs > beta_band_lims(1) & ersp_freqs < beta_band_lims(2));
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
        %- custom alg from epoch info
        % inds = find(cellfun(@(x) any(strcmp(x,conds{cond_i})),{EEG.epoch.eventcond}));
        %- use info from .icatimef
        inds = find(cellfun(@(x) strcmp(x,conds{cond_i}),{trialinfo.cond}));
        %- imu data
        datx = squeeze(EEG.data(body_posx,:,inds)); % AP (anteroposterior)
        daty = squeeze(EEG.data(body_posy,:,inds)); % ML (mediolateral)
        datz = squeeze(EEG.data(body_posz,:,inds)); % UD (Up-Down)
        %## SANITY CHECK        
        %{
        figure;
        title(conds{cond_i})
        hold on;
        inds = randi(size(datx,2),10);
        for s_i = inds
            plot3(datx(:,s_i),daty(:,s_i),datz(:,s_i));
        end
        hold off;
        %-
        figure;
        hold on;
        for s_i = inds
            plot(EEG.times,datx(:,s_i),'DisplayName','x');
            plot(EEG.times,daty(:,s_i),'DisplayName','y');
            plot(EEG.times,datz(:,s_i),'DisplayName','z');
        end
        hold off;
        %}
        %## ADVANCED EPOCH'ING?        
        for s_i = 1:length(inds)-1
            %## USING .ICATIMEF INFORMATION
            log_s = true;
            rhs_1 = timewarp_trial(inds(s_i),1)*TIME_FACTOR;
            rhs_2 = timewarp_trial(inds(s_i),5)*TIME_FACTOR;
            lhs_1 = timewarp_trial(inds(s_i),2)*TIME_FACTOR;
            lto_1 = timewarp_trial(inds(s_i),3)*TIME_FACTOR;
            rto_1 = timewarp_trial(inds(s_i),4)*TIME_FACTOR;
            %- calculate pelvis movement
            EVENT_ERR = 1.001;
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
            vec = [daty(rhsi_1:lhsi_1,s_i)];
            m1 = min(vec);
            m2 = max(vec);
            ml_mm_exc = sqrt((m2-m1)^2);
            %* min max exc (ml_contra)
            vec = [daty(lhsi_1:rhsi_2,s_i)];
            m1 = min(vec);
            m2 = max(vec);
            ml_contra_mm_exc = sqrt((m2-m1)^2);
            %* min max exc (ap)
            vec = [datx(rhsi_1:lhsi_1,s_i)];
            m1 = min(vec);
            m2 = max(vec);
            ap_mm_exc = sqrt((m2-m1)^2);
            %* min max exc (ap_contra)
            vec = [datx(lhsi_1:rhsi_2,s_i)];
            m1 = min(vec);
            m2 = max(vec);
            ap_contra_mm_exc = sqrt((m2-m1)^2);
            %* min max exc (ud)
            vec = [datz(rhsi_1:lhsi_1,s_i)];
            m1 = min(vec);
            m2 = max(vec);
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
            % if s_i == 1
            %     fprintf('step width: %0.1g',ml_mm_exc);
            % elseif s_i < length(inds)-2
            %     fprintf(', %0.1g',ml_mm_exc);
            % else
            %     fprintf('\n');
            % end
            %## LOG STRUCT           
            %- struct
            if log_s %&& all(~isempty([rhs_1,rhs_2,lhs_1,lto_1,rto_1]))
                for i = 1:length(CLUSTER_PICS)
                    cl_i = CLUSTER_PICS(i);
                    comp_n = comp_arr(1,comp_arr(3,:) == cl_i);
                    if ~isempty(comp_n)
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
                        steps_struct(cnt).stride_n = inds(s_i);
                        steps_struct(cnt).rhs_1 = rhs_1;
                        steps_struct(cnt).rhs_2 = rhs_2;
                        steps_struct(cnt).lhs_1 = lhs_1;
                        steps_struct(cnt).lto_1 = lto_1;
                        steps_struct(cnt).rto_1 = rto_1;
                        steps_struct(cnt).gait_cycle_dur = rhs_2 - rhs_1;
                        steps_struct(cnt).stance_dur = rto_1 - rhs_1;
                        steps_struct(cnt).swing_dur = rhs_2 - rto_1;
                        steps_struct(cnt).double_sup_dur = rto_1 - lhs_1;
                        steps_struct(cnt).single_sup_dur = lhs_1 - lto_1;
                        steps_struct(cnt).step_dur = rhs_2 - lhs_1;
                        %-
                        steps_struct(cnt).step_width = ml_exc;
                        steps_struct(cnt).ap_exc = ap_exc;
                        steps_struct(cnt).ud_exc = ud_exc;
                        %-
                        steps_struct(cnt).step_width_mm = ml_mm_exc;
                        steps_struct(cnt).step_width_contra_mm = ml_contra_mm_exc;
                        steps_struct(cnt).ap_contra_mm_exc = ap_contra_mm_exc;
                        steps_struct(cnt).ap_exc_mm = ap_mm_exc;
                        steps_struct(cnt).ud_exc_mm = ud_mm_exc;
                        %-
                        steps_struct(cnt).step_width_gc_mm = ml_mm_gc_exc;
                        steps_struct(cnt).ap_exc_gc_mm = ap_mm_gc_exc;
                        steps_struct(cnt).ud_exc_gc_mm = ud_mm_gc_exc;
                        %- ADD GC PARAMS!!!!!!!

                        % step_structs(cnt).eeg_psd = squeeze(eeg_psd(:,cnt,comp_n));
                        % step_structs(cnt).ersp_freqs = ersp_freqs;
                        %- RHS-LHS of the next step >>> timewarpms::{'RHS','LTO','LHS','RTO','RHS'}
                        stance_1 = ersp_times >= timewarpms(1) & ersp_times <= timewarpms(2);
                        swing_1 = ersp_times >= timewarpms(2) & ersp_times <= timewarpms(3);
                        stance_2 = ersp_times >= timewarpms(3) & ersp_times <= timewarpms(4);
                        swing_2 = ersp_times >= timewarpms(4) & ersp_times <= timewarpms(5);
                        %- calculate mean across time THEN frequency.
                        % steps_struct(cnt).avg_theta = mean(eeg_ersp(tband,step_time_inds,inds(s_i),comp_n),[2,1]);
                        % steps_struct(cnt).avg_alpha = mean(eeg_ersp(aband,step_time_inds,inds(s_i),comp_n),[2,1]);
                        % steps_struct(cnt).avg_beta = mean(eeg_ersp(bband,step_time_inds,inds(s_i),comp_n),[2,1]);
                        % int_order = [1,2]; % mean across time then freq
                        int_order = [2,1]; % mean acros freq then time
                        %## IN STRIDE
                        steps_struct(cnt).avg_theta_stance1 = mean(eeg_ersp(tband,stance_1,inds(s_i),comp_n),int_order);
                        steps_struct(cnt).avg_alpha_stance1 = mean(eeg_ersp(aband,stance_1,inds(s_i),comp_n),int_order);
                        steps_struct(cnt).avg_beta_stance1 = mean(eeg_ersp(bband,stance_1,inds(s_i),comp_n),int_order);
                        steps_struct(cnt).avg_theta_stance2 = mean(eeg_ersp(tband,stance_2,inds(s_i),comp_n),int_order);
                        steps_struct(cnt).avg_alpha_stance2 = mean(eeg_ersp(aband,stance_2,inds(s_i),comp_n),int_order);
                        steps_struct(cnt).avg_beta_stance2 = mean(eeg_ersp(bband,stance_2,inds(s_i),comp_n),int_order);
                        %-
                        steps_struct(cnt).avg_theta_swing1 = mean(eeg_ersp(tband,swing_1,inds(s_i),comp_n),int_order);
                        steps_struct(cnt).avg_alpha_swing1 = mean(eeg_ersp(aband,swing_1,inds(s_i),comp_n),int_order);
                        steps_struct(cnt).avg_beta_swing1 = mean(eeg_ersp(bband,swing_1,inds(s_i),comp_n),int_order);
                        steps_struct(cnt).avg_theta_swing2 = mean(eeg_ersp(tband,swing_2,inds(s_i),comp_n),int_order);
                        steps_struct(cnt).avg_alpha_swing2 = mean(eeg_ersp(aband,swing_2,inds(s_i),comp_n),int_order);
                        steps_struct(cnt).avg_beta_swing2 = mean(eeg_ersp(bband,swing_2,inds(s_i),comp_n),int_order);
                        %## POST
                        steps_struct(cnt).avg_theta_post_stance1 = mean(eeg_ersp(tband,stance_1,inds(s_i+1),comp_n),int_order);
                        steps_struct(cnt).avg_alpha_post_stance1 = mean(eeg_ersp(aband,stance_1,inds(s_i+1),comp_n),int_order);
                        steps_struct(cnt).avg_beta_post_stance1 = mean(eeg_ersp(bband,stance_1,inds(s_i+1),comp_n),int_order);
                        steps_struct(cnt).avg_theta_post_stance2 = mean(eeg_ersp(tband,stance_2,inds(s_i+1),comp_n),int_order);
                        steps_struct(cnt).avg_alpha_post_stance2 = mean(eeg_ersp(aband,stance_2,inds(s_i+1),comp_n),int_order);
                        steps_struct(cnt).avg_beta_post_stance2 = mean(eeg_ersp(bband,stance_2,inds(s_i+1),comp_n),int_order);
                        %-
                        steps_struct(cnt).avg_theta_post_swing1 = mean(eeg_ersp(tband,swing_1,inds(s_i+1),comp_n),int_order);
                        steps_struct(cnt).avg_alpha_post_swing1 = mean(eeg_ersp(aband,swing_1,inds(s_i+1),comp_n),int_order);
                        steps_struct(cnt).avg_beta_post_swing1 = mean(eeg_ersp(bband,swing_1,inds(s_i+1),comp_n),int_order);
                        steps_struct(cnt).avg_theta_post_swing2 = mean(eeg_ersp(tband,swing_2,inds(s_i+1),comp_n),int_order);
                        steps_struct(cnt).avg_alpha_post_swing2 = mean(eeg_ersp(aband,swing_2,inds(s_i+1),comp_n),int_order);
                        steps_struct(cnt).avg_beta_post_swing2 = mean(eeg_ersp(bband,swing_2,inds(s_i+1),comp_n),int_order);
                        cnt = cnt + 1;
                    end
                end
            end
        end
        
        %- grab condition
        inds_cond = strcmp({steps_struct.cond_char},conds{cond_i});
        tmp_ss = steps_struct(inds_cond);
        %- per condition measures (step dur)
        musd = mean([tmp_ss.step_dur]);
        stdsd = std([tmp_ss.step_dur]);
        tmp = [tmp_ss.step_dur];
        %-
        vals = sqrt((tmp-musd).^2);
        vals = num2cell(vals);
        [steps_struct(inds_cond).var_step_dur_1] = deal(vals{:});
        %- COV-like
        vals = stdsd./sqrt((tmp-musd).^2);
        vals = num2cell(vals);
        [steps_struct(inds_cond).var_step_dur_2] = deal(vals{:});
        %-
        % vals =sqrt((tmp-musd).^2)./stdsd;
        % vals = num2cell(vals);
        % [steps_struct(inds_cond).var_step_dur_3] = deal(vals{:});
        %- per condition measures (step width)
        musd = mean([tmp_ss.step_width]);
        stdsd = std([tmp_ss.step_width]);
        tmp = [tmp_ss.step_width];
        %-
        vals = sqrt((tmp-musd).^2);
        vals = num2cell(vals);
        [steps_struct(inds_cond).var_step_width_1] = deal(vals{:});
        %- COV-like
        vals = stdsd./sqrt((tmp-musd).^2);
        vals = num2cell(vals);
        [steps_struct(inds_cond).var_step_width_2] = deal(vals{:});
        fprintf('\n');
    end
    %- sanity check
    % p2 = steps_struct.
    % fig = plot_txf_conds_tftopo(p2,ersp_times,ersp_freqs,[],[],...
    %     'PLOT_STRUCT',PLOT_STRUCT);
    %- assign to struct
    % vals = num2cell(vals)';
    % [steps_struct(inds_cond)] = deal(vals{:});
    %- per design measures
    %## SAVE    
    steps_struct = steps_struct(~cellfun(@isempty,{steps_struct.comp_n}));
    steps_struct = steps_struct(~cellfun(@isempty,{steps_struct.avg_theta_post_stance1}));
    steps_struct = steps_struct(~cellfun(@isempty,{steps_struct.step_width}));
    steps_struct = struct2table(steps_struct);
    par_save(steps_struct,EEG.filepath,'ls_eeg_ersp_struct_trialbase_fix.mat');
    % par_save(steps_struct,EEG.filepath,'ls_eeg_struct_indvfooof.mat');
    %- turn off when doing parallel processing.
    % data_store = cat(1,data_store,step_structs);
    EEG.data = 'eeg_imu_epochs.fdt';
    tmp_alleeg{subj_i} = EEG;
end
%% (RESAVE STUDY) ====================================================== %%
%- remove bugged out subjects
tmp_alleeg = tmp_alleeg(~cellfun(@isempty,tmp_alleeg));
%- clear data for memory
% for i = 1:length(tmp_alleeg)
%     tmp_alleeg{i}.data = 'eeg_imu_epochs.fdt';
% end
%## BOOKKEEPING (i.e., ADD fields not similar across EEG structures)
tmp_alleeg = util_resolve_struct(tmp_alleeg);
%##
[STUDY, ALLEEG] = std_editset([],tmp_alleeg,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name','d_eeg_ersp_imuls_sts_anl',...
                                'filename','d_eeg_ersp_imuls_sts_anl',...
                                'filepath',save_dir);
[STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
%## POP PARAMS
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
STUDY = pop_statparams(STUDY,args{:});
args = eeglab_struct2args(ERSP_PARAMS);
STUDY = pop_erspparams(STUDY,args{:});
STUDY.cache = [];
for des_i = 1:length(STUDY_DESI_PARAMS)
    [STUDY] = std_makedesign(STUDY,[],des_i,STUDY_DESI_PARAMS{des_i}{:});
end
%## SAVE
[STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                                        STUDY.filename,STUDY.filepath,...
                                        'RESAVE_DATASETS','off');
%% (TEST CHANGES ACROSS SPEEDS) ======================================== %%
data_store = table.empty;
groups = {'H1000','H2000','H3000'};
group_name = {'younger_adults','older_high_function','older_low_function'};
for i = 1:length(STUDY.datasetinfo)
    try
        tmp = par_load(STUDY.datasetinfo(i).filepath,'ls_eeg_ersp_struct_trialbase_fix.mat');
        vals = regexp(tmp.subj_char{1},'[nN]?[hH](\d)\d*','tokens');
        vals = double(string(vals{1}{1}));
        tt = repmat(vals,[1,size(tmp,1)]);
        tmp.group_n = tt';
        tt = repmat(groups(vals),[1,size(tmp,1)]);
        [tmp.group_char] = tt'; %deal(tt{:});
        tt = repmat(group_name(vals),[1,size(tmp,1)]);
        [tmp.group_name] = tt'; %deal(tt{:});
        data_store = [data_store; tmp];
    catch e
        fprintf('Couldn''t load %s...\n%s',STUDY.datasetinfo(i).subject,getReport(e));
    end
    % data_store = cat(1,data_store,tmp);
end
% steps_struct = data_store;
% data_store = struct.empty;
writetable(data_store,[save_dir filesep 'step_by_step_eeg_ersp_trialbase_table_fix.xlsx']);
par_save(data_store,save_dir,'step_by_step_eeg_ersp_trialbase_table_fix.mat');
% writetable(data_store,[save_dir filesep 'step_by_step_eeg_table_indvfooof.xlsx']);
% par_save(data_store,save_dir,'step_by_step_eeg_table_indvfooof.mat');



