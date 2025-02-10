% Create study for cluster ICs on HPG. This code only works for cluster without
% using ERSP. Precompute ERSP needed to be done on Hipergator
% Chang Liu - 2021-11-23 - V1
% Chang Liu - 2022-07-23 - use on hipergator

%Run after DIPFIT and epoching. This puts all good dipoles into a study for
%clustering and ERSP plotting.

%   NJacobsen notes
%   When timewarping data, save values as EEG.timewarp = timewarp;
%   EEG.timewarp.medianlatency = median(timewarp.latencies(:,:));%Warping to the median latency of my 5 events
%   By default, std_ersp will use the median of all subject's
%   timewarp.latencies(:,:) as 'timewarpms' unless individual subject
%   warpto is indiciated using 'timewarpms', 'subject tw matrix'

% 2023-12-30 TH trying to update this script to cycle through different
% event epoched datasets to make an individual STUDY

%{
2024-09-09: TH creating study with subjects that have enough epochs with
the correct leading events prior to the window. subject = 11, folder_n=8; 

create_study_PS is for use with separated baseline and trial sets being put into the same
study. 

11-26-24: TH on the side check results of study with and without flagged
subjects. folder n=2

12/14/24: n=3 re running epoching in rnTW folder to create a study with
reduced sample size to compare ERSPs


CURRENT VERSION: create_study_PS.m creates study for just pong trials. uses
_allCond to precompute measures
 
%}

%%   Setup
close all; clearvars
if isunix
    PATH_ROOT = ['/blue' filesep 'dferris' filesep 'theresa.hauge'];
elseif ispc
    PATH_ROOT =['M:' filesep 'theresa.hauge'];
end
addpath([PATH_ROOT filesep 'eeglab2021.0'],[PATH_ROOT filesep 'HY' filesep 'scripts' filesep 'code_Preprocess'])
P3_HY_config_params;

eeglab;
pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0, ...
    'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1);

Process_S = 0;
Process_PS = 1;
Process_TS = 0;
Process_BLandPS = 0;
Process_D = 0; %condStr = '_D';
n = 4; % subjPt n=3
%{
if Process_S == 1
    conds = {'Pong','TableTennis'};
    condStr = '_S';
    condList = {'_S'};
elseif Process_D == 1
    conds = {'PongDual','TennisDual','PongSingle','TennisSingle'};
elseif Process_PS == 1
    conds = {'Pong'};
    condStr = '_PS';
    condList = {'PS'};
elseif Process_TS == 1
    conds = {'TableTennis'};
    condStr = '_TS';
    condList = {'_BL','_TS'};
elseif Process_BLandPS ==1
    conds = {'Baseline','Pong'};
    condStr = '_PS';
    condList = {'BL','PS'};
end
%}
%-DEFINE FILEPATHS AND FOLDERS
condStr = '_PS';
STUDY_date = '20241113';
folder_name =  ['BATCH-3-ICA_',STUDY_date,condStr];
filepath = ['STUDY-',STUDY_date,'_PS'];
Study_folder = folder_name;
% studyname = ['STUDY-',STUDY_date,condStr,'.study'];

load_folder = fullfile(PATH_ROOT,'HY','Data','P3_preprocess',filepath, folder_name); brain_score = 8;
save_study_folder = fullfile(load_folder,'BATCH-5-Epoch'); %, num2str(n));
save_epoch_folder = fullfile(load_folder,'BATCH-5-Epoch');
load_IC_rej = fullfile(load_folder,'BATCH-4-ICRejection');
if ~exist(save_study_folder,'dir'), mkdir(save_study_folder), end
%{
suffix_list = {'_SubEp','_ResEp',...
    '_SubEpWarp_GrAvg_m1p40_rnTW','_ResEpWarp_GrAvg_m1p40_rn',...
    '_SubServeEp','_ResServeEp',...
    '_SubServeEpWarp_GrAvg','_ResServeEpWarp_GrAvg'};
% updating suffix_list to include other events

suffix_list = {'_SubEp','_ResEp',...
    '_SubEpWarp_GrAvg_m1p40_rn','_ResEpWarp_GrAvg_m1p40_rn'...
    '_SubServeEp','_ResServeEp',...
    '_SubServeEpWarp_GrAvg','_ResServeEpWarp_GrAvg',...
    '_Wall_T','_Wall_B',...
    '_Wall_TWarp','_Wall_BWarp',...
    '_SubjPtEp','_ResPtEp',...
    '_SubjPtWarp','_ResPtWarp',...
    '_BL_PS_epoched_m1p40',... %17 
    '_SubPooledEp','_ResPooledEp',... %18 19
    '_SubPoolWarp','_ResPoolWarp',...
    '_ResWallSub_Warp','_SubWallRes_Warp'};

numSTUDY = length(suffix_list);
%}
% last_suff = '_numTrials';

%% new information for multiple designs in one study, also assuming there you are not making individual studies for each design
CondGroup = {'ResPt','SubjPt'...
    'RallyRSR','RallySRS',...
    'ResServeS','SubServeR',...
    'RallyRWallS','RallySWallR',...
    'ResWallPt','SubWallPt'};
set_suffix = '_MergeWarp_GrAvg';
design_flags = [1 2; 3 4;5 6; 7 8; 9 10;...
    5 1; 6 1; 7 9; 8 10];
condGen = {'Points','Rally','Serve','WallRally','WallPoint',...
    'RServePt','SServePt','RWallPt','SWallPt'};
n=4;
condList = {CondGroup{design_flags(n,:)}};
CondDesign = condGen{n};
f_design = strcat('design',num2str(n));
studyname = ['STUDY-',STUDY_date,'_',CondDesign,'.study'];
save_study_folder = fullfile(save_epoch_folder,f_design);

SUBJ_RM = {'S21_0214','S24_0225','S36_0408'}; SUBJ_PICS = all_subjStr(~contains(all_subjStr,SUBJ_RM));
for i = 1:length(SUBJ_PICS)
    subjStr = SUBJ_PICS{i};
    %     for des_i = 1:length(design_flags)

    EEG = pop_loadset('filename',strcat(subjStr,'_',f_design,set_suffix,'.set'),...
        'filepath', fullfile(save_epoch_folder,f_design,subjStr));
    % _PS and _TS only data
    for b = 1:length(EEG.event)
        EEG.event(b).trialCond = string(EEG.event(b).trialCond);
        EEG.event(b).trialName = string(EEG.event(b).trialName);
    end

    [ALLEEG,EEG,~] = eeg_store(ALLEEG,EEG,0);
    %     end
end

%{
for study_i = 23
    set_suffix = '_MergeWarp_GrAvg'; %char(suffix_list(study_i));
%     set_nm = set_suffix(2:end);
    %     studyname = ['STUDY-',STUDY_date,condStr,set_suffix,last_suff,'.study'];
    save_study = fullfile(save_study_folder);
    if ~exist(save_study,'dir'), mkdir(fullfile(save_study)), end
    TWon = contains(set_suffix,'Warp'); disp(TWon)
    if TWon
        TW = 1; groupmedian_timewarpms = 1;
    else
        TW = 0; groupmedian_timewarpms = 0;
    end

    SUBJ_RM = {'S21_0214','S24_0225','S36_0408'}; SUBJ_PICS = all_subjStr(~contains(all_subjStr,SUBJ_RM));
%     SUBJ_RM = {'S01_0220','S17_0206','S18_0208','S19_0210','S25_0301','S28_0307','S29_0309','S24_0225','S36_0408'};
%     SUBJ_PICS = all_subjStr(~contains(all_subjStr,SUBJ_RM));
%     load(fullfile(save_epoch_folder,set_nm,'viable_subs.mat')); keep_subj = KEEP_SUBJ;
    for i = 1:length(SUBJ_PICS)
        subjStr = SUBJ_PICS{i};
        if ~exist(fullfile(save_epoch_folder,set_nm,num2str(n),subjStr),'dir')==7
            mkdir(fullfile(fullfile(save_epoch_folder,set_nm,num2str(n),subjStr)))
        end
        %         for a = 1:length(condList)
        %             condStr =  % condList{a};
        EEG = pop_loadset('filename',strcat(subjStr,f_design,set_suffix,'.set'),...
            'filepath', fullfile(save_epoch_folder,f_design,subjStr));
        % _PS and _TS only data
        for b = 1:length(EEG.event)
            EEG.event(b).trialCond = string(EEG.event(b).trialCond);
            EEG.event(b).trialName = string(EEG.event(b).trialName);
        end
        mkdir(fullfile(load_folder,'BATCH-5-Epoch',set_nm,num2str(n),subjStr))
%         EEG.cond = 'Pong';
        EEG = pop_saveset(EEG,'filename',strcat(subjStr,condStr,set_suffix, '.set'),'filepath',...
                fullfile(load_folder,'BATCH-5-Epoch',set_nm,num2str(n),subjStr),'savemode','twofiles','version','7.3');
       
        % EEG.cond = EEG.condition;
        % EEG.condition = [];
        %             EEG = pop_saveset(EEG,'filename',strcat(subjStr,condStr,set_suffix, '.set'),'filepath',fullfile(load_folder,'BATCH-5-Epoch',set_nm,subjStr));
        [ALLEEG,EEG,~] = eeg_store(ALLEEG,EEG,0);
%{
        % % adding in loadset for the baseline epoched data
        if Process_BLandPS
            EEG = pop_loadset('filename',strcat(subjStr,'_PS','_BL_epoched.set'),...
                'filepath', fullfile(save_epoch_folder,set_nm,'numTrials',subjStr));
            EEG.setname = [subjStr,'_BL_PS_post_epochs_rmBL'];
            for b = 1:length(EEG.event)
                EEG.event(b).trialCond = string(EEG.event(b).trialCond);
                EEG.event(b).trialName = string(EEG.event(b).trialName);
            end
            EEG = pop_saveset(EEG,'filename',strcat(subjStr,'_PS','_BL_epoched.set'),...
                'filepath', fullfile(save_epoch_folder,set_nm,'numTrials',num2str(n),subjStr));
            [ALLEEG,EEG,~] = eeg_store(ALLEEG,EEG,0);
        end
%}
        disp('next subject')
        %         end
    end
%}

DO_FEM_DIPFIT = 1;
win_epoch = [-1 4]; %ALLEEG(1).etc.epochinfo.win_epoch;
STUDY_PICS = SUBJ_PICS;
%% from Jacob Salminen, setting paths to BEM even if FEM was computed
for ii = 1:length(ALLEEG)
    EEG = ALLEEG(ii);
    tmp = strsplit(path,';');
    % tmp = strsplit(path,';');
    b1 = regexp(tmp,'eeglab','end');
    b2 = tmp(~cellfun(@isempty,b1));
    PATH_EEGLAB = b2{1}(1:b1{1});
    fprintf('EEGLAB path: %s\n',PATH_EEGLAB);
    %- set default paths for boundary element head model
    PATH_EEGLAB_BEM  = [PATH_EEGLAB filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM' filesep];
    MNI_MRI = [PATH_EEGLAB_BEM filesep 'standard_mri.mat'];
    MNI_VOL = [PATH_EEGLAB_BEM filesep 'standard_vol.mat'];
    %     MNI_CHAN_1005 = [PATH_EEGLAB_BEM filesep 'elec' filesep 'standard_1005.elc'];
    MNI_CHAN_1005 = [PATH_EEGLAB_BEM 'elec' filesep 'standard_1005.elc'];
    COORD_TRANSFORM_MNI = [0 0 0 0 0 -1.5708 1 1 1];
    tmp = [];
    if ~isfield(EEG.dipfit,'coord_transform')
        EEG.dipfit.coord_transform = [0 0 0 0 0 0 1 1 1]; %COORD_TRANSFORM_MNI;
        tmp = [tmp, 'added default coord_transform; '];
    end
    if ~isfield(EEG.dipfit,'mrifile')
        EEG.dipfit.mrifile = MNI_MRI;
        tmp = [tmp, 'added default mrifile; '];
    end
    if ~isfield(EEG.dipfit,'hdmfile')
        EEG.dipfit.hdmfile = MNI_VOL;
        tmp = [tmp, 'added default hdmfile; '];
    end
    if ~isfield(EEG.dipfit,'coordformat')
        EEG.dipfit.coordformat = 'MNI';
        tmp = [tmp, 'added default coordformat; '];
    end
    if ~isfield(EEG.dipfit,'chanfile')
        EEG.dipfit.chanfile = MNI_CHAN_1005;
        tmp = [tmp, 'added default chanfile; '];
    end
    if ~isfield(EEG.dipfit,'chansel')
        EEG.dipfit.chansel = (1:EEG.nbchan);
        tmp = [tmp, 'added default chansel; '];
    end

    for i = 1:length(STUDY_PICS) %(SUBJ_PICS)
        if DO_FEM_DIPFIT
            fprintf('Trying Finite Element Model Dipfit...\n');
            try
                EEG.dipfit.model;
            catch e
                fprintf('A valid dipfit model needs to be calculated for FEM dipfit analysis...\n');
                fprintf(['error. identifier: %s\n',...
                    'error. %s\n',...
                    'error. on subject %s\n'],e.identifier,e.message,EEG.subject);
                exit();
            end
            [EEG] = fem_eeglab_dipfit(EEG,COORD_TRANSFORM_MNI,MNI_MRI,MNI_VOL,MNI_CHAN_1005);
        end
    end
    ALLEEG(ii) = EEG;
    %     [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG,EEG,i);
end


%% Change the folder path
% ------------------------------------
disp('FEM complete')

subDirStudy = 1;
Process_Imagined = 0;

% -----------------------------------
% mkdir(fullfile(save_study_folder,num2str(subDirStudy)));

load_study = 0;
if load_study
    tic
    [STUDY, ALLEEG] = pop_loadstudy('filename', studyname, 'filepath', save_study);
    toc
end
%     % THERESA this needs to be updated but for now comment it out
%     temp = dir(fullfile(load_folder,'BATCH-5-Epoch','*S*'));
%     temp = temp(~contains({temp.name},'BATCH'));
%     all_subjStr = {temp(:).name};

%% SET PARAMS - if you aren't time warping around 2+ events, you set TW = 0 and groupmedian_timewarpms =0.
%     TW = 0; %Time warping option %THERESA -- this is set to 0 while you are still working on response locked and stimulus locked ERSP and ERP analysis
preclust = 1; %1 - precluster study; 0 - no preclustering
clust = 1; %1 - precluster study; 0 - no preclustering
precomp_nonERSPs = 0; %1 - pulls up precompute gui; 0 - no gui
precompute_ERSP = 1; %1 - precompute ersps; 0 - don't precompute ersps
precompute_spec_scalp = 1;

erspComp='full'; %'light' - quicker computation; 'full' - with usual parameters (takes longer)
showClusterPlotGUI=1; %1 - show cluster plot gui at end; 0 - don't show it (and clear study)
%     groupmedian_timewarpms = 0;%NJacobsen; warp each subject's
%                           tw matrix to the entire group's median event
%                           latencies [1=ON], or use individual subject's
%                           median event latencies [0=OFF]. TW must be ON
%                           for this setting to do anything
clustering_weights.dipoles = 1; %other options: 5, 1; % this value seems to be different across researchers. not sure why but will determine what is best.
clustering_weights.scalp = 0; %other options = 5, 0; trying 1 and 1
cluster_alg = 'kmeans';
do_multivariate_data = 1;

%% this is the original script that I copied from CL, might be more straightforward to understand the logic.
% make the study
disp('starting study creation')
% alleeg2 = ALLEEG;
EEG=[]; ALLEEG=[]; CURRENTSET=[]; ALLCOM=[]; CURRENTSTUDY = []; STUDY = [];

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;  %Theresa commented out to try
% to fix datasetinfo issues.
subjToAnalyze = {};
indx = 1;
disp('Check do not use eeglab on R drive');
save_study = fullfile(save_study_folder);
% Create Study - type of study should pull conditions you want
for i = 1:length(STUDY_PICS)    %(SUBJ_PICS)               % THERESA check this loop start -- should it be 2 or something different?
    subjStr = STUDY_PICS{i}; %SUBJ_PICS
    subjStr = subjStr(1:8);
    %for ii = 1:length(condList)
    % condStr = condList{i};
    if  exist(fullfile(save_study,subjStr,[subjStr,'_',f_design,set_suffix,'.set']),'file')==2 % exist(fullfile(save_epoch_folder,set_nm,'numTrials',subjStr,[subjStr,condStr,set_suffix,'_numTrials','.set']),'file')==2
        temp_EEG = importdata(fullfile(save_study,subjStr,[subjStr,'_',f_design,set_suffix,'.set'])); % temp_EEG = importdata(fullfile(save_epoch_folder,set_nm,'numTrials',subjStr,[subjStr,condStr,set_suffix,'_numTrials','.set']));
        disp('LET''S GET IT STARTED IN HAH')
        dipole_validxyz = vertcat(temp_EEG.dipfit.model(:).posxyz);
        dipole_valid = find(~isnan(dipole_validxyz(:,1)));

        % Remove bad ICs first
        try
            load(fullfile(load_IC_rej,subjStr,[subjStr,condStr,'_ICRej.mat']))
        catch
            load(fullfile(load_IC_rej,subjStr,[subjStr,'_ICRej.mat']))
        end
        badcompInds = setdiff(1:size(Output_ICRejection.IC_all_brain,1),find((Output_ICRejection.IC_all_brain >= brain_score & Output_ICRejection.IC_all_brain ~= 9)));
        brainIC = find(Output_ICRejection.IC_all_brain >= brain_score & Output_ICRejection.IC_all_brain ~= 9);
        if isempty(brainIC)
            continue
        end

        dipole_valid_all = intersect(dipole_valid,brainIC);
        % ------ I used eeglab 2021 on hipergator
        %             [STUDY ALLEEG] = std_editset( STUDY, ALLEEG,...
        %                 'name',studyname,...
        %                 'commands',...
        %                 {{'index',indx,'load',fullfile(save_study,subjStr,[subjStr,condStr,set_suffix,last_suff,'.set']),... %min50 numTrials_50
        %                 'comps',dipole_valid_all',...
        %                 'subject',[subjStr],...
        %                 'condition',condStr}},'updatedat','on'); % THERESA CHANGED SUBJECT NAME SO IT'S CONSISTENT

        [STUDY, ALLEEG] = std_editset( STUDY, ALLEEG,...
            'name',studyname,...
            'commands',...
            {{'index',indx,'load',...
            fullfile(save_study,subjStr,[subjStr,'_',f_design,set_suffix,'.set']),... %min50 numTrials_50
            'comps',dipole_valid_all',...
            'subject',subjStr,...
            'session',1,...
            'run',1,...
            'condition',condGen{n}}},...   % 'cond',condGen{n}
            'updatedat','on');
        %                 {'index',indx+1,'load',...
        %                 fullfile(save_study,subjStr,[subjStr,'_PS_BL_epoched','.set']),... %min50 numTrials_50
        %                 'comps',dipole_valid_all',...
        %                 'subject',subjStr,...
        %                 'session',1,...
        %                 'run',1,...
        %                 'condition','Baseline',...
        %                 'conds','BL'}},...
        %                 'updatedat','on');
        %                 {'index',indx,'run',2},...
        %                 {'index',indx+1,'run',1}},...
        % THERESA CHANGED SUBJECT NAME SO IT'S CONSISTENT, added conds to see if it helps
        indx = indx + 1;
        subjToAnalyze{i} = subjStr; disp(subjToAnalyze)
        %             indx = indx + 2;
    end

    disp('LET''S GET IT STARTED IN HERE')
    %end
end

disp(subjToAnalyze)
[STUDY ALLEEG] = std_checkset(STUDY, ALLEEG);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
% eeglab redraw

%     %Make STUDY design - idk why it is erroring out, but i'm commenting out
%     %because I can change my study design later, i'm pretty sure
[STUDY] = std_makedesign(STUDY, ALLEEG, 1, 'subjselect', subjToAnalyze,...
    'variable1','cond','values1',condList); %changed from condList to be consistent down the line
% Save STUDY
[STUDY, ALLEEG] = pop_savestudy( STUDY, EEG, 'filename',studyname,'filepath',save_study);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);

% find best warp for all subjects
TW=1;
if TW == 1
    allWarpTo = []; SubjListEpoched=STUDY.subject;
    for set_i = 1:length(SubjListEpoched)
        if Process_PS || Process_TS
            % Assuming ALLEEG is a structure array and SubjListEpoched is a cell array
            %                 subjectMatches = find(contains({ALLEEG.subject}, SubjListEpoched{set_i}));
            %                 setnameMatches = find(~cellfun('isempty', regexp({ALLEEG(subjectMatches).filename}, 'Warp')));
            %                 matchingIndices = subjectMatches(setnameMatches);
            %
            %                 allWarpTo(set_i,:) = ALLEEG(matchingIndices).timewarp.warpto; %stack subject specific median event latencies - using this for single condition sets _PS and _TS

            allWarpTo(set_i,:) = ALLEEG(set_i).timewarp.grandAvgWarpTo; %temporary fix to deal with multiple conditions with different timewarp lengths... _PS and _TS
        elseif Process_S
            allWarpTo(set_i,:) = ALLEEG(set_i).timewarp.grandAvgWarpTo; % THERESA update to deal with multiple conditions with different timewarp lengths... 20240630 use for _S
        end
    end
    grandAvgWarpTo = median(allWarpTo); disp(grandAvgWarpTo) % trying this 20240630
    % grandAvgWarpTo = mean(allWarpTo);
    disp(['Grand average (across all subj) warp to: ', num2str(grandAvgWarpTo)]);
else
    disp('you are not timewarping your data, hopefully this will resolve some of your issues')
    allWarpTo = []; SubjListEpoched=STUDY.subject;
    for set_i = 1:length(SubjListEpoched)
        allWarpTo = []; % changed from allWarpTo(set_i,:) = [];
    end
    grandAvgWarpTo = [];
end

%%add in a block of code to get printout of all epoch numbers for study
%%split by conditions THERESA 2024-04-23


%% params
%     if isunix
%         MainDirectory = fullfile('/blue/','dferris/','theresa.hauge/','HY/','Data/','P3_preprocess/',filepath,'/', folder_name,'/BATCH-5-Epoch/',set_nm,'/','numTrials'); %decide which study to use
%     elseif ispc
%         MainDirectory = fullfile('Z:','theresa.hauge','HY','Data','P3_preprocess',filepath, folder_name,'BATCH-5-Epoch/',set_nm,'/','numTrials');
%     end
MainDirectory = save_study; %fullfile(save_epoch_folder,set_nm,'numTrials',num2str(n));
%choose which params to calculate and send into the precompute function
%20231215 Theresa turning off
params.CalcERP = true;
params.CalcSpectra = true;
params.CalcERSP = true;
%Add ability to calc scalp topo!
params.chanOrComp = 'comp'; %default = 'comp'. Use 'chan' to examine subset of raw sensor data (experimental)
params.overwriteBool.Spectra = 'on';
params.overwriteBool.ERSP = 'on';
params.baselineCorrectBool = false;                                     %if baseline correction is set to true, will remove average spectral power in ersps.
%- additional params
params.condStr = 'Pong'; %condStr
params.win_epoch = win_epoch;
%     params.padratio = 4;
if TW==1
    params.spec_timerange = [(win_epoch(1)*1000) (win_epoch(2)*1000)]; %[-500 4500]; %
    %         params.spec_timerange = [-1000 2500]; %[-500 4500]; %
    params.timewarpBool = true; %
else
    params.spec_timerange = [-1000 1500];
    params.timewarpBool = false; % turned off unless timewarping
end

%THERESA add spec_params? JSAL makes a SPEC_PARAMS and ERSP_PARAMS
%structures, which is very neat imo...
SPEC_PARAMS = struct('freqrange',[1,200],...
    'subject','',...
    'specmode','psd',...
    'freqfac',4,...
    'logtrials','on',...
    'comps','all',...
    'plot_freqrange',[4,100],...
    'plot_ylim',[-35,-8],...
    'subtractsubjectmean','on',...
    'plotmode','normal');
ERSP_PARAMS = struct('subbaseline','off',...
    'timerange',[],...
    'ersplim',[win_epoch(1), win_epoch(2)],...
    'freqfac',4,...
    'cycles',[3,0.5],... % 3,0.8
    'padratio',4,...
    'freqrange',[1,200]);
params.SPEC_PARAMS = SPEC_PARAMS; params.ERSP_PARAMS = ERSP_PARAMS;
%% do precompute
addpath([PATH_ROOT filesep 'HY' filesep 'scripts' filesep 'code_Preprocess' filesep '7_precomp']) % if
% running on HPG
% addpath 'Z:\\theresa.hauge\HY\scripts\code_Preprocess\7_precomp' % if running locally in Matlab
PrecomputeMeasuresGivenStudy_function_allCond( MainDirectory,studyname, params);

% end
%% ======================== SUBFUNCTIONS START HERE ========================================= %%
function [EEG] = fem_eeglab_dipfit(EEG,COORD_TRANSFORM_MNI,MNI_MRI,MNI_VOL,MNI_CHAN_1005)
tmp = [];
if ~isfield(EEG.dipfit,'coord_transform')
    EEG.dipfit.coord_transform = [0 0 0 0 0 0 1 1 1]; %COORD_TRANSFORM_MNI;
    tmp = [tmp, 'added default coord_transform; '];
end
if ~isfield(EEG.dipfit,'mrifile')
    EEG.dipfit.mrifile = MNI_MRI;
    tmp = [tmp, 'added default mrifile; '];
end
if ~isfield(EEG.dipfit,'hdmfile')
    EEG.dipfit.hdmfile = MNI_VOL;
    tmp = [tmp, 'added default hdmfile; '];
end
if ~isfield(EEG.dipfit,'coordformat')
    EEG.dipfit.coordformat = 'MNI';
    tmp = [tmp, 'added default coordformat; '];
end
if ~isfield(EEG.dipfit,'chanfile')
    EEG.dipfit.chanfile = MNI_CHAN_1005;
    tmp = [tmp, 'added default chanfile; '];
end
if ~isfield(EEG.dipfit,'chansel')
    EEG.dipfit.chansel = (1:EEG.nbchan);
    tmp = [tmp, 'added default chansel; '];
end
end
function [EEG] = custom_update_chanlocs(EEG,chanlocs_fPath)
tmp = load(chanlocs_fPath);
chanlocs_new = tmp.chanlocs_new;
nodatchans_new = tmp.nodatchans_new;
%- update the EEG electrode locations
% Be cautious that not all electrodes are EEG
% Sanity check: if we have 120 electrodes digitized
fprintf('Found total of %i electrodes',length(chanlocs_new));
for p = 1:length(chanlocs_new)
    elec_idx = find(strcmpi(chanlocs_new(p).labels,{EEG.chanlocs(:).labels}));
    if ~isempty(elec_idx)
        % update all available fields
        EEG.chanlocs(elec_idx).X = chanlocs_new(p).X;
        EEG.chanlocs(elec_idx).Y = chanlocs_new(p).Y;
        EEG.chanlocs(elec_idx).Z = chanlocs_new(p).Z;
        EEG.chanlocs(elec_idx).theta = chanlocs_new(p).theta;
        EEG.chanlocs(elec_idx).radius = chanlocs_new(p).radius;
        EEG.chanlocs(elec_idx).sph_theta = chanlocs_new(p).sph_theta;
        EEG.chanlocs(elec_idx).sph_phi = chanlocs_new(p).sph_phi;
    end
end
% Add fiducials location
if isempty(EEG.chaninfo.nodatchans)
    EEG.chaninfo.nodatchans = nodatchans_new;
end
EEG = eeg_checkchanlocs(EEG); % check the consistency of the chanloc structure
end

function [EEG] = bem_eeglab_dipfit(EEG,COORD_TRANSFORM_MNI,MNI_MRI,MNI_VOL,MNI_CHAN_1005,DIP_NUM,DIP_PLOT)
%## FIT DIPOLES TO HEADMODEL IF NEEDED
fprintf('MNI pop_dipfit_settings...\n');
%## Check Files
tmp = [];
if ~isfield(EEG.dipfit,'coord_transform')
    EEG.dipfit.coord_transform = COORD_TRANSFORM_MNI;
    tmp = [tmp, 'added default coord_transform; '];
end
if ~isfield(EEG.dipfit,'mrifile')
    EEG.dipfit.mrifile = MNI_MRI;
    tmp = [tmp, 'added default mrifile; '];
end
if ~isfield(EEG.dipfit,'hdmfile')
    EEG.dipfit.hdmfile = MNI_VOL;
    tmp = [tmp, 'added default hdmfile; '];
end
if ~isfield(EEG.dipfit,'coordformat')
    EEG.dipfit.coordformat = 'MNI';
    tmp = [tmp, 'added default coordformat; '];
end
if ~isfield(EEG.dipfit,'chanfile')
    EEG.dipfit.chanfile = MNI_CHAN_1005;
    tmp = [tmp, 'added default chanfile; '];
end
if ~isfield(EEG.dipfit,'chansel')
    EEG.dipfit.chansel = (1:EEG.nbchan);
    tmp = [tmp, 'added default chansel; '];
end
EEG.dipfit.comment = tmp;

%- pop_multifit.m
%- DIPFIT (see. ft_dipolefitting())
fprintf('pop_multifit...\n');
EEG = pop_multifit(EEG,[],'dipoles',DIP_NUM,'dipplot',DIP_PLOT);
%     dipfit = EEG.dipfit;
%     save([fPath filesep sprintf('%s_dipfit_mni.mat',EEG.subject)],'dipfit');
end



%{
    % find numTrials needed to reduce down the number of trials for a given
    % subject to be the same
        gameList = {'Pong','TableTennis'};
        epochNums = zeros(length(ALLEEG),2);
        for g = 1:length(ALLEEG)
            eeg2=ALLEEG(g);
            for cond_i = 1:length(gameList)
                CurrentCond = gameList{cond_i};
                if any(strcmpi(CurrentCond,{eeg2.event.cond}))
                    subsetEEG = pop_selectevent(eeg2, 'cond',CurrentCond,'deleteevents','off','deleteepochs','on','invertepochs','off');
                    disp(['selected ' CurrentCond ' trials']);
                end
                fprintf('There are %u epochs for %s within condition %s\n',length(subsetEEG.epoch),subsetEEG.subject,CurrentCond )
                epochNums(g,cond_i) = length(subsetEEG.epoch);
            end
            numTrials(g) = min(epochNums(g,:))';
        end
        disp(epochNums)
        min(epochNums(:,1)); min(epochNums(:,2));
    
    % only select number of trials that will mean the same number between
    % conditions for each participant
    for f = 1:length(SUBJ_PICS)
        subjStr = SUBJ_PICS{f};
        EEG2 = ALLEEG(f);
        for cond_i = 1:length(gameList)                 % go through each condition and select random epochs
            CurrentCond = gameList{cond_i};
            if any(strcmpi(CurrentCond,{EEG2.event.cond}))
                subsetEEG = pop_selectevent(EEG2, 'cond',CurrentCond,'deleteevents','off','deleteepochs','on','invertepochs','off');
                disp(['selected ' CurrentCond ' trials']);
                subsetEEG = pop_select(subsetEEG,'trial',randsample(1:size(subsetEEG.data,3),numTrials(f)));
            end
            [ALLEEG,subsetEEG,CURRENTSET] = eeg_store(ALLEEG,subsetEEG,0);
        end
        RemergedEEG = pop_mergeset( ALLEEG(end-1:end), 1:size(ALLEEG(end-1:end),2),1); %if running multiple conditions together, probably comment out for _PS /_TS only
        RemergedEEG.filename= [subjStr,condStr,set_suffix,'_numTrials'];
        mkdir(fullfile(load_folder,'BATCH-5-Epoch',set_nm,'numTrials',subjStr)) % SAVE Epoched dataset
        RemergedEEG = pop_saveset( RemergedEEG, 'filepath', fullfile(load_folder,'BATCH-5-Epoch',set_nm,'numTrials',subjStr), 'filename', RemergedEEG.filename);
        [ALLEEG,RemergedEEG,CURRENTSET] = eeg_store(ALLEEG,RemergedEEG,0);
    end
%}

