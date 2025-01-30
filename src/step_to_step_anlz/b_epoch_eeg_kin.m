%   Project Title: MIM OA & YA SPEED & KINETICS ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/_bash_sh_files/run_sts_b_epoch_eeg_kin.sh

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
ADD_ALL_SUBMODS = true;
%## Determine Working Directories
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = fileparts(SCRIPT_DIR); % change this if in sub folder
        SRC_DIR = STUDY_DIR;
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
    SRC_DIR = STUDY_DIR;
end
%## Add Study, Src, & Script Paths
addpath(SRC_DIR);
addpath(STUDY_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
% [SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('test');
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca_speed');
%## UNIT TEST
% SUBJ_PICS = {{'H1004','H1007'},{'H2013','H2018_FU'},{'H3120','NH3129'}};
% GROUP_NAMES = {'H1000''s','H2000''s','H3000''s'};
% SUBJ_ITERS = {[1,2],[1,2],[1,2]};
%% (PARAMETERS) ======================================================== %%
%- epoching params
DO_STANDARD_TRIALS = false;
MIN_STANDARD_TRIALS = 100;
%## EPOCH PARAMS
EPOCH_PARAMS = struct('epoch_method','timewarp',...
    'percent_overlap',0,...
    'epoch_event_char','RHS',...
    'epoch_time_lims',[-0.5,4.5],...
    'baseline_time_lims',[-0.5,4.5-2],...
    'tw_stdev',3,...
    'tw_events',{{'RHS','LTO','LHS','RTO','RHS'}},...
    'path_ext','gait_epoched',...
    'cond_field','cond',...
    'appx_cond_len',3*60,...
    'slide_cond_chars',{{}},...
    'gait_trial_chars',{{'0p25','0p5','0p75','1p0','flat','low','med','high'}},...
    'rest_trial_char',{{}},...
    'do_recalc_epoch',true,...
    'do_combine_trials',true);
%% (PATHS) ========================================================== %%
%- datset name
DATA_SET = 'MIM_dataset';
%## PREPROCESSED ICA 
ICA_DIR_FNAME = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
SUBJ_FNAME_REGEX = 'cleanEEG_EMG_HP3std_iCC0p65_iCCEMG0p4_ChanRej0p7_TimeRej0p4_winTol10';

%## PROCESSED STUDY
STUDY_DNAME = '01192025_mim_yaoa_nopowpow_crit_speed';
STUDY_FNAME = 'all_comps_study';

%## SAVE INFO
KIN_DNAME_EXT = 'kin_eeg_anl';
STUDY_FNAME_EPOCH = 'kin_eeg_epoch_study';
studies_dir = [PATHS.data_dir filesep DATA_SET filesep '_studies'];
ica_data_dir = [studies_dir filesep ICA_DIR_FNAME]; % JACOB,SAL(02/23/2023)
save_dir = [studies_dir filesep sprintf('%s',STUDY_DNAME)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% (LOAD PROCESSED STUDY)
%## LOAD STUDY
if ~ispc
    tmp = load('-mat',[save_dir filesep sprintf('%s_UNIX.study',STUDY_FNAME)]);
    STUDY_PROC = tmp.STUDY;
else
    tmp = load('-mat',[save_dir filesep sprintf('%s.study',STUDY_FNAME)]);
    STUDY_PROC = tmp.STUDY;
end
%% Store fNames and fPaths
subj_chars          = [SUBJ_PICS{:}];
fNames              = cell(1,length([SUBJ_PICS{:}]));
fPaths              = cell(1,length([SUBJ_PICS{:}]));
chanlocs_fPaths     = cell(1,length([SUBJ_PICS{:}]));
dipfit_norm_fPaths  = cell(1,length([SUBJ_PICS{:}]));
%-
for subj_i = 1:length(subj_chars)
    fPaths{subj_i} = [ica_data_dir filesep subj_chars{subj_i} filesep 'clean'];
    fNames{subj_i} = sprintf('%s_%s.set',subj_chars{subj_i},SUBJ_FNAME_REGEX);
    fprintf('Loading Subject %s\n',subj_chars{subj_i})
    if ~exist([fPaths{subj_i} filesep fNames{subj_i}],'file')
        fprintf('No .set file found...\n')
        dipfit_norm_fPaths{subj_i} = [];
    else
        chanlocs_fPaths{subj_i} = [PATHS.data_dir filesep DATA_SET filesep subj_chars{subj_i} filesep 'MRI' filesep 'CustomElectrodeLocations.mat'];
        dipfit_norm_fPaths{subj_i} = [fPaths{subj_i} filesep 'dipfit_fem_norm_ants.mat'];
        fprintf('ICA Exists: %i\n',(exist([fPaths{subj_i} filesep fNames{subj_i}],'file') && exist([fPaths{subj_i} filesep 'W'],'file')))
        fprintf('Normalized DIPFIT Exists: %i\n',exist(dipfit_norm_fPaths{subj_i},'file'));
    end
end
%- remove subjects without a dipole fit
inds = logical(cellfun(@(x) exist(x,'file'),dipfit_norm_fPaths));
chanlocs_fPaths = chanlocs_fPaths(inds);
dipfit_norm_fPaths = dipfit_norm_fPaths(inds);
fPaths = fPaths(inds);
fNames = fNames(inds);
subj_chars = subj_chars(inds);
%% (PROC 1)
tmp_alleeg = cell(length(fPaths),1);
%## PATHING UPDATES
% path(unix_genpath([PATHS.submods_dir filesep 'Gait Tracking With x-IMU']),path);
path(unix_genpath([PATHS.submods_dir filesep 'gait_tracking_w_imu']),path);
%##
parfor subj_i = 1:length(subj_chars)
    tmp_epoch_params = EPOCH_PARAMS;
    EEG = []; 
    tmp_biom = {};
    tmp_study_proc = STUDY_PROC;
    cont_proc = true;
    try
        trial_fpath = [PATHS.data_dir filesep DATA_SET filesep subj_chars{subj_i} filesep 'EEG' filesep 'Trials'];
        %## Find All Trials of Interest
        fileList_TM = dir([trial_fpath filesep 'TM*.set']);
        fileList_SP = dir([trial_fpath filesep 'SP*.set']);
        fileList_Rest = dir([trial_fpath filesep 'Rest.set']);
        fileList = [fileList_Rest; fileList_TM; fileList_SP];
        tmp_savedir = [save_dir filesep subj_chars{subj_i} filesep 'kin_eeg_valid_plots'];
        if ~exist(tmp_savedir,'dir')
            mkdir(tmp_savedir);
        end
        tmp_biom = cell(size(fileList,1),1);
        for trial_i = 1:size(fileList,1)
            %- Load trial
            trial_i = 10; trial_i = 15;
            EEG = pop_loadset('filename',fileList(trial_i).name,'filepath',trial_fpath);
            %(11/22/2024) JS, trying to fix a bug where the IMU channels don't get
            %loaded if I load using the pop_loadset function. Actually they get
            %loaded but they are full of nan's not sure where the entire trace gets
            %removed
            %(12/01/2024) JS, see eeg_getdatact.m (lines 203-254'ish) for
            %loading procedure on channel level data
            %{
            EEG        = load('-mat', [trial_fpath filesep fileList(trial_i).name]);
            EEG        = EEG.EEG;
            fid         = fopen([trial_fpath filesep tmpt.data], 'r', 'ieee-le');
            data        = fread(fid, [tmpt.nbchan Inf], 'float32');
            data        = reshape(data, size(data,1), tmpt.pnts, tmpt.trials);
            data        = single(data);
            EEG.data   = data;
            % data = tmp.data;
            fprintf('percent data points nan: %0.2g\n',(sum(logical(isnan(data)),'all')/(size(data,1)*size(data,2))*100));
            
            back_imu_chans = find(contains({tmpt.chanlocs.labels},'Back','IgnoreCase',true));
            ls_chans = find(contains({tmpt.chanlocs.labels},'GRF','IgnoreCase',true));
            EEG = pop_select(tmpt,'channel',sort([back_imu_chans,ls_chans]));
            % data = data(sort([back_imu_chans,ls_chans]),:);
            %}
            fprintf('percent data points nan: %0.2g\n',(sum(logical(isnan(EEG.data)),'all')/(size(EEG.data,1)*size(EEG.data,2))*100));
            back_imu_chans = find(contains({EEG.chanlocs.labels},'Back','IgnoreCase',true));
            ls_chans = find(contains({EEG.chanlocs.labels},'GRF','IgnoreCase',true));
            EEG.data(isnan(EEG.data)) = 0;
            EEG = pop_select(EEG,'channel',sort([back_imu_chans,ls_chans]));
            %- calculate body and world positions
            %(12/11/2024) JS, changing the position grabbing alg. to on a
            %per trial basis instead of at the total experiment level.
            if any(subj_i==[1,40,70])
                export_res = 900;
            else
                export_res = 300;
            end
            EEG = imu_get_pos_coords(EEG);
            nn = strsplit(fileList(trial_i).name,'.')
            fh = imu_valid_plots(EEG,tmp_savedir,sprintf('%s_%s_',subj_chars{subj_i},nn{1}),...
                'EXPORT_RES',export_res);
            close(fh);
            % tmp_biom.nbchan = length(tmp_biom.chanlocs);
            tmp_biom{trial_i} = EEG;
            % pop_eegplot( EEG, 1, 1, 1);
            fprintf('percent data points nan: %0.2g\n',sum(logical(isnan(tmp_biom{trial_i}.data)),'all')/(size(tmp_biom{trial_i}.data,1)*size(tmp_biom{trial_i}.data,2)));
        end        
        %##
        tmp_biom = tmp_biom(~cellfun(@isempty,tmp_biom));
        %## BOOKKEEPING (i.e., ADD fields not similar across EEG structures)
        tmp_biom = util_resolve_struct(tmp_biom);
        tmp_biom = pop_mergeset(tmp_biom,1:length(tmp_biom), 0);
        fh = imu_valid_plots(tmp_biom,trial_fpath,sprintf('%s_all_',subj_chars{subj_i}));
        close(fh);
        %- check data
        fprintf('percent data points nan: %0.2g\n',(sum(logical(isnan(tmp_biom.data)),'all')/(size(tmp_biom.data,1)*size(tmp_biom.data,2))*100));
        %## CALCULATE POSITION COORDSf
        % tmp_biom = imu_get_pos_coords(tmp_biom);
        %(12/01/2024) JS, fixed bug where NaN values were causing the filter to
        %NaN entire trace.
        %(12/11/2024) JS, no longer doing this for now.
        body_chans = find(contains({tmp_biom.chanlocs.labels},'body','IgnoreCase',true));
        world_chans = find(contains({tmp_biom.chanlocs.labels},'world','IgnoreCase',true));
        tmp_biom.nbchan = length(tmp_biom.chanlocs);
        tmp_biom = pop_select(tmp_biom,'channel',sort([body_chans,world_chans]));
    catch e
        fprintf('%s) Error. in BioM struct construction...\n\n%s',subj_chars{subj_i},getReport(e))
        cont_proc = false;
        % tmp_alleeg{subj_i} = [];
    end
    %## EPOCHING
    if cont_proc
        try
            %## LOAD ICA INFORMATION
            [~,EEG,~] = eeglab_loadICA(fNames{subj_i},fPaths{subj_i});        
            %-
            eeg_chans = find(strcmpi('EEG',{EEG.chanlocs.type}));
            biom_chans = find(strcmpi('BioM',{tmp_biom.chanlocs.type}));
            %- check data
            fprintf('percent data points nan: %0.2g\n',sum(logical(isnan(tmp_biom.data)),'all')/(size(tmp_biom.data,1)*size(tmp_biom.data,2)));
            
            %## MERGE BIOM & EEG DATA
            tmp = EEG;
            tmp.data = [tmp.data(:,:);
                tmp_biom.data(:,tmp.etc.clean_artifacts.clean_sample_mask)];
            %- add in chanlocs information
            n_eeg_chans = length(eeg_chans);
            for chan = 1:length(sort(biom_chans))
                %- fix not consistent chan info
                tmp.chanlocs(n_eeg_chans+chan).urchan = [];
                tmp.chanlocs(n_eeg_chans+chan).ref = [];
                tmp.chanlocs(n_eeg_chans+chan).sph_theta_besa = [];
                tmp.chanlocs(n_eeg_chans+chan).sph_phi_besa = [];
                %-
                if ~isfield(tmp_biom.chanlocs(chan),'sph_theta_besa')
                    tmp_biom.chanlocs(chan).sph_theta_besa = [];
                    tmp_biom.chanlocs(chan).sph_phi_besa = [];
                end
                tmp.chanlocs(n_eeg_chans+chan) = tmp_biom.chanlocs(chan);
            end
            tmp.nbchan = size(tmp.data,1);
            %## REJECT NON BRAIN COMPONENTS
            % THRESH_BRAIN_SCORE = 8;
            % tmp_icrej = load([save_dir filesep subj_chars{subj_i} filesep 'ICA' filesep sprintf('%s_allcond_ICA_TMPEEG.set',subj_chars{subj_i})],'-mat');
            ind = strcmp({tmp_study_proc.datasetinfo.subject},subj_chars{subj_i})
            tmp_icrej = load([tmp_study_proc.datasetinfo(ind).filepath filesep tmp_study_proc.datasetinfo(ind).filename],'-mat');
            biom_chans = find(strcmpi('BioM',{tmp.chanlocs.type}));
            ics_keep = tmp_icrej.etc.urreject.ic_keep;
            urreject = tmp_icrej.etc.urreject;
            %- robust recreation of rejection crit
            %{
            icrej_crit = load([save_dir filesep subj_chars{subj_i} filesep 'ICA' filesep sprintf('%s_ICRej.mat',subj_chars{subj_i})],'-mat');
            reject_struct = icrej_crit.Output_ICRejection;
            chk = (reject_struct.IC_all_brain >= 8 & reject_struct.IC_all_brain ~= 9);
            chk_w_powpow = unique(cat(1,find(chk),reject_struct.IC_powpow_rej));
            tmp_bad = setdiff(find((1:size(ALLEEG.icaweights,1))),chk_w_powpow);
            tmp_good = chk_w_powpow';
            urreject = [];
            urreject.crit = reject_struct;
            urreject.ic_keep = tmp_good;
            urreject.ic_rej = tmp_bad;
            urreject.dipfit = [];
            %}
            tmp = pop_select(tmp,'channel',sort([ics_keep, biom_chans]));
            tmp.etc.urreject = urreject;
            %## REASSIGN
            EEG = tmp;
            tmp = struct.empty; % memory clean-up
            EEG.subject = subj_chars{subj_i};
            
            %## DIPFIT INFO
            fprintf('Trying to load .mat file...\n');
            out = load(dipfit_norm_fPaths{subj_i});
            EEG.dipfit = out.SAVEVAR;
            
            %## EPOCH INFORMATION
            if isfield(EEG.event,'trialName')
                EEG.event = rmfield(EEG.event,'trialName');
            end
            if isfield(EEG.event,'channel')
                EEG.event = rmfield(EEG.event,'channel');
            end
            if isfield(EEG.event,'code')
                EEG.event = rmfield(EEG.event,'code');
            end
            if isfield(EEG.event,'bvtime')
                EEG.event = rmfield(EEG.event,'bvtime');
            end
            if isfield(EEG.event,'bvmknum')
                EEG.event = rmfield(EEG.event,'bvmknum');
            end
            if isfield(EEG.event,'datetime')
                EEG.event = rmfield(EEG.event,'datetime');
            end
            EEG.etc.chanlocs_biom_eeg = EEG.chanlocs;
            [ALLEEG,timewarp_struct] = mim_imu_ls_trial_parser(EEG, ...
                'EPOCH_PARAMS',tmp_epoch_params);
            if DO_STANDARD_TRIALS
                for i = 1:length(ALLEEG)
                    %##
                    inds = randi(length(ALLEEG(i).epoch),MIN_STANDARD_TRIALS,1);
                    %##
                    tmp_all = ALLEEG(i);
                    tmp_all = pop_selectevent(tmp_all,'event',inds);
                    while length(tmp_all.epoch)~=MIN_STANDARD_TRIALS
                        tmp_all = ALLEEG(i);
                        inds = randi(length(tmp_all.epoch),MIN_STANDARD_TRIALS,1)
                        tmp_all = pop_selectevent(tmp_all,'event',inds);
                    end
                    ALLEEG(i) = tmp_all;
                end
            end
            ALLEEG = pop_mergeset(ALLEEG,1:length(ALLEEG),1);
            
            %## TIMEWARP ACROSS CONDITIONS
            timewarp = make_timewarp(ALLEEG,tmp_epoch_params.tw_events,'baselineLatency',0, ...
                    'maxSTDForAbsolute',inf,...
                    'maxSTDForRelative',inf);
            %- subject specific warpto (later use to help calc grand avg warpto across subjects)
            timewarp.warpto = nanmedian(timewarp.latencies);        
            goodepochs  = sort([timewarp.epochs]);
            %- probably not needed? 
            sedi = setdiff(1:length(ALLEEG.epoch),goodepochs);
            %- reject outlier strides
            ALLEEG = pop_select(ALLEEG,'notrial',sedi);
            %- store timewarp structure in EEG
            ALLEEG.timewarp = timewarp;
            %- store condition-by-conditino timewarpings
            ALLEEG.etc.timewarp_by_cond = timewarp_struct;
            %## FINAL DATA CHECK
            % pop_eegplot( ALLEEG, 1, 1, 1);
            %## STRUCT EDITS
            groups = {'H1000','H2000','H3000'};
            group_name = {'younger_adults','older_high_function','older_low_function'};
            vals = regexp(ALLEEG.subject,'[nN]?[hH](\d)\d*','tokens');
            vals = double(string(vals{1}{1}));
            ALLEEG.urevent = []; % might be needed
            ALLEEG.etc.epoch_params = tmp_epoch_params;
            ALLEEG.group = groups{vals};
            ALLEEG.session = 1;
            ALLEEG.run = 1;
            %## SAVE
            if ~exist([save_dir filesep ALLEEG.subject filesep KIN_DNAME_EXT],'dir')
                mkdir([save_dir filesep ALLEEG.subject filesep KIN_DNAME_EXT]);
            end
            ALLEEG = pop_saveset(ALLEEG,'filepath',[save_dir filesep ALLEEG.subject filesep KIN_DNAME_EXT],'filename','eeg_imu_epochs.set');
            tmp_alleeg{subj_i} = ALLEEG;
        catch e
            fprintf('%s) Error. in EEG/BioM epoching process...\n\n%s',subj_chars{subj_i},getReport(e))
            % tmp_alleeg{subj_i} = [];
        end
    else
        fprintf('%s) Error. no epoching performed...\n\n',subj_chars{subj_i});
        % tmp_alleeg{subj_i} = [];
    end
end
%%
%- remove bugged out subjects
tmp_alleeg = tmp_alleeg(~cellfun(@isempty,tmp_alleeg));
%## BOOKKEEPING (i.e., ADD fields not similar across EEG structures)
tmp_alleeg = util_resolve_struct(tmp_alleeg);
disp(tmp_alleeg);

%##
[STUDY, ALLEEG] = std_editset([],tmp_alleeg,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name',STUDY_FNAME_EPOCH,...
                                'filename',STUDY_FNAME_EPOCH,...
                                'filepath',save_dir);
[STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
STUDY.etc.a_epoch_process.epoch_params = EPOCH_PARAMS;
[STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                                        STUDY.filename,STUDY.filepath,...
                                        'RESAVE_DATASETS','on');