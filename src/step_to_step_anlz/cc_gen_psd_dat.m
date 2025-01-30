%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/_bash_sh_files/run_sts_cc_gen_psd_dat.sh

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
    SRC_DIR = fileparts(SCRIPT_DIR); % change this if in sub folder
end
%## Add Study, Src, & Script Paths
addpath(SCRIPT_DIR);
addpath(SRC_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
%% (PARAMETERS) ======================================================== %%
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
%% (PATHS) ========================================================== %%
%- datset name
DATA_SET = 'MIM_dataset';

%## PROCESSED STUDY
STUDY_DNAME = '01192025_mim_yaoa_nopowpow_crit_speed';

%## SAVE INFO
STUDY_FNAME_EPOCH = 'kin_eeg_epoch_study';
studies_dir = [PATHS.data_dir filesep DATA_SET filesep '_studies'];
save_dir = [studies_dir filesep sprintf('%s',STUDY_DNAME)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
% if ~ispc
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME '_UNIX.study'],'filepath',save_dir);
% else
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME '.study'],'filepath',save_dir);
% end
%## LOAD STUDY
if ~ispc
    tmp = load('-mat',[save_dir filesep sprintf('%s_UNIX.study',STUDY_FNAME_EPOCH)]);
    STUDY = tmp.STUDY;
else
    tmp = load('-mat',[save_dir filesep sprintf('%s.study',STUDY_FNAME_EPOCH)]);
    STUDY = tmp.STUDY;
end
%% (ERSP PLOT PREP) PREPARE STUDYFILE FOR EXTRACTION (BLACK-HAWK DOWN!)
%## ersp plot per cluster per condition
STUDY = pop_statparams(STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);

SPEC_PARAMS.subtractsubjectmean = 'on';
STUDY = pop_specparams(STUDY,'subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
    'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed',...
    'plotconditions','together','ylim',SPEC_PARAMS.plot_ylim,'plotgroups','together');
%% (PRECOMPUTE MEASURES) COMPUTE ERSPs
subject_chars = {STUDY.datasetinfo.subject};
% for subj_i = 1:length(ALLEEG)
parfor subj_i = 1:length(STUDY.datasetinfo)
    %- TEMP PARAMS
    tmp_study = STUDY;
    tmp_spec_params = SPEC_PARAMS;
    % tmp_epoch_params = EPOCH_PARAMS;
    %-
    tt = tic;
    EEG = pop_loadset('filepath',tmp_study.datasetinfo(subj_i).filepath,...
        'filename',tmp_study.datasetinfo(subj_i).filename);
    fprintf('Running Subject %s\n',EEG.subject);
    EEG = eeg_checkset(EEG,'loaddata');
    if isempty(EEG.icaact)
        fprintf('%s) Recalculating ICA activations\n',EEG.subject);
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
    end
    %## GENERATE PSDS
    tmp_study_calc = tmp_study;
    %- overrride datasetinfo to trick std_precomp to run.
    tmp_study_calc.datasetinfo = tmp_study_calc.datasetinfo(subj_i);
    tmp_study_calc.datasetinfo(1).index = 1;
    [~, ~] = std_precomp(tmp_study_calc, EEG,...
                'components',...
                'recompute','on',...
                'spec','on',...
                'scalp','on',...
                'savetrials','on',...
                'specparams', ...
                {'specmode',tmp_spec_params.specmode, ...
                'freqfac',tmp_spec_params.freqfac,...
                'freqrange',tmp_spec_params.freqrange, ...
                'logtrials',tmp_spec_params.logtrials});
    fprintf('Done calculating PSDs for %s\n',EEG.subject);
end
% %% (SPCA IMPLEMENTATION?) ============================================== %%
% DEF_EPOCH_PARAMS = struct('epoch_method','timewarp',...
%     'percent_overlap',0,...
%     'epoch_event_char','RHS',...
%     'epoch_time_lims',[-0.5,4.5],...
%     'baseline_time_lims',[-0.5,4.5-2],...
%     'tw_stdev',3,...
%     'tw_events',{{'RHS','LTO','LHS','RTO','RHS'}},...
%     'path_ext','gait_epoched',...
%     'cond_field','cond',...
%     'appx_cond_len',3*60,...
%     'slide_cond_chars',{{}},...
%     'gait_trial_chars',{{'0p25','0p5','0p75','1p0','flat','low','med','high'}},...
%     'rest_trial_char',{{}},...
%     'do_recalc_epoch',true);
% %% (PATHS)
% %- datset name & 
% DATA_SET = 'MIM_dataset';
% ICA_DIR_FNAME = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
% STUDY_DNAME = '01102025_mim_yaoa_spca_calcs';
% STUDY_FNAME_GAIT = 'spca_gait_epoch_study';
% STUDY_FNAME_REST = 'spca_rest_slide_study';
% %-
% studies_dir = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
% save_dir = [studies_dir filesep sprintf('%s',STUDY_DNAME)];
% ica_data_dir = [studies_dir filesep ICA_DIR_FNAME];
% %- create new study directory
% if ~exist(save_dir,'dir')
%     mkdir(save_dir);
% end
% %% LOAD STUDY
% if ~ispc
%     tmp = load('-mat',[save_dir filesep sprintf('%s_UNIX.study',STUDY_FNAME_GAIT)]);
%     STUDY_GAIT = tmp.STUDY;
% else
%     tmp = load('-mat',[save_dir filesep sprintf('%s.study',STUDY_FNAME_GAIT)]);
%     STUDY_GAIT = tmp.STUDY;
% end
% %%
% subj_ind = 1;
% epoched_fPath = strsplit(STUDY.datasetinfo(subj_ind).filepath,filesep);
% icatimef_f = [strjoin(epoched_fPath,filesep) filesep sprintf('%s.icaspec',STUDY.datasetinfo(subj_ind).subject)];
% %-
% EEG = pop_loadset('filepath',STUDY.datasetinfo(subj_ind).filepath,'filename',STUDY.datasetinfo(subj_ind).filename);
% %- load .icatimef load-in parameters
% fprintf('Loading Time-Frequency Data...\n');
% tmp = load(icatimef_f,'-mat');
% %-
% fn = fieldnames(tmp);
% inds = find(contains(fieldnames(tmp),'comp'));
% test = tmp.(fn{inds(1)});
% %-
% eeg_psd = zeros(size(test,1),size(test,2),length(inds),'double');
% for i = 1:length(inds)
%     eeg_psd(:,:,i) = tmp.(fn{inds(i)}); % freq x epoch x chan
% end
% %-
% % nolog_eeg_psd = 10.^(eeg_psd/10);
% freqs = tmp.freqs;
% trialinfo = tmp.trialinfo;
% %-
% figure;
% chan_i = 1;
% plot(freqs,squeeze(eeg_psd(:,:,chan_i)));
% %%
% tmp_epoch_params = DEF_EPOCH_PARAMS;
% fprintf('Loading Gait Frequency Data...\n');
% epoched_fPath = STUDY_GAIT.datasetinfo(subj_ind).filepath;
% icaspec_f = [epoched_fPath filesep sprintf('%s.icaspec',STUDY.datasetinfo(subj_ind).subject)];
% %-
% EEG = pop_loadset('filepath',STUDY_GAIT.datasetinfo(subj_ind).filepath,'filename',STUDY_GAIT.datasetinfo(subj_ind).filename);
% %- load .icaspec load-in parameters
% tmp = load(icaspec_f,'-mat');
% %-
% fn = fieldnames(tmp);
% inds = find(contains(fieldnames(tmp),'comp'));
% test = tmp.(fn{inds(1)});
% %-
% spca_eeg_psd = zeros(size(test,1),size(test,2),length(inds),'double');
% for i = 1:length(inds)
%     spca_eeg_psd(:,:,i) = tmp.(fn{inds(i)}); % freq x epoch x chan
% end
% spca_freqs = tmp.freqs;
% spca_trialinfo = tmp.trialinfo;
% %-
% figure;
% chan_i = 1;
% plot(freqs,squeeze(spca_eeg_psd(:,:,chan_i)));