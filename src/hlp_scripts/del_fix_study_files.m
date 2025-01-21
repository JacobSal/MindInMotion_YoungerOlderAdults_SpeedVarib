%   Project Title: MIM OA & YA SPEED & KINETICS ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/_bash_sh_files/run_a_epoch_process.sh

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
% [SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca_speed');
% [SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('ya');
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
%- study name
STUDY_DNAME = '10172024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej_speed';
% STUDY_FNAME = 'kin_eeg_epoch_study';
% ANALYSIS_DNAME = 'kin_eeg_step_to_step';
STUDY_FNAME_CONT = 'all_comps_study';
STUDY_FNAME_EPOCH = 'epoch_study';
%% (PATHS)
studies_fpath = [PATHS.data_dir filesep DATA_SET filesep '_studies'];
%- load cluster
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
% cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb10'];
% cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb3'];
cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb5'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
cluster_k_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
%% (EPOCH STUDY RECREATION) ============================================ %%
if ~ispc
    tmp = load('-mat',[studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep sprintf('%s_UNIX.study',STUDY_FNAME_EPOCH)]);
    tmp_study = tmp.STUDY;
else
    tmp = load('-mat',[studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep sprintf('%s.study',STUDY_FNAME_EPOCH)]);
    tmp_study = tmp.STUDY;
end
%## LOAD SUBJECTS 
dirext = ['0p250p50p751p0flatlowmedhigh'];
ALLEEG = cell(length(tmp_study.datasetinfo),1);
parfor subj_i = 1:length(tmp_study.datasetinfo)
    subj_char = tmp_study.datasetinfo(subj_i).subject;
    fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep subj_char filesep dirext];
    EEG = pop_loadset('filepath',fpath, ...
        'filename',tmp_study.datasetinfo(subj_i).filename);
    ALLEEG{subj_i} = pop_saveset(EEG, ...
        'filepath',fpath, ...
        'filename',tmp_study.datasetinfo(subj_i).filename);

end
%## REMAKE AND SAVE STUDY
ALLEEG = ALLEEG(~cellfun(@isempty,ALLEEG));
ALLEEG = util_resolve_struct(ALLEEG);
[STUDY, ALLEEG] = std_editset([],ALLEEG,...
                                    'updatedat','off',...
                                    'savedat','off',...
                                    'name',STUDY_FNAME_EPOCH,...
                                    'filename',STUDY_FNAME_EPOCH,...
                                    'filepath',[studies_fpath filesep sprintf('%s',STUDY_DNAME)]);
STUDY.etc = tmp_study.etc;
[STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
[STUDY,~] = parfunc_save_study(STUDY,ALLEEG,...
                                STUDY.filename,STUDY.filepath,...
                                'RESAVE_DATASETS','on');  
%% (CLUSTER STUDY RECREATION) ========================================== %%
if ~ispc
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_NAME)]);
    tmp_study = tmp.STUDY;
else
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s.study',CLUSTER_STUDY_NAME)]);
    tmp_study = tmp.STUDY;
end
%## LOAD SUBJECTS 
dirext = ['0p250p50p751p0flatlowmedhigh'];
ALLEEG = cell(length(tmp_study.datasetinfo),1);
for subj_i = 1:length(tmp_study.datasetinfo)
    subj_char = tmp_study.datasetinfo(subj_i).subject;
    fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep subj_char filesep dirext];
    ALLEEG{subj_i} = pop_loadset('filepath',fpath, ...
                    'filename',tmp_study.datasetinfo(subj_i).filename);
end
%## REMAKE AND SAVE STUDY
ALLEEG = ALLEEG(~cellfun(@isempty,ALLEEG));
ALLEEG = util_resolve_struct(ALLEEG);
[STUDY, ALLEEG] = std_editset([],ALLEEG,...
                                    'updatedat','off',...
                                    'savedat','off',...
                                    'name',STUDY_FNAME_EPOCH,...
                                    'filename',STUDY_FNAME_EPOCH,...
                                    'filepath',[studies_fpath filesep sprintf('%s',STUDY_DNAME)]);
STUDY.etc = tmp_study.etc;
[STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
[STUDY,~] = parfunc_save_study(STUDY,ALLEEG,...
                                STUDY.filename,STUDY.filepath,...
                                'RESAVE_DATASETS','on');  
%% (DELETE UNUSED FILES) =============================================== %%
% for subj_i = 1:length(subj_chars)
%     tt = tic();
%     tmp_save_dir = [save_dir filesep subj_chars{subj_i}];
%     %-
%     fpaths_del = {[tmp_save_dir filesep 'ICA' filesep sprintf('%s_allcond_ICA_TMPEEG.set',subj_chars{subj_i})], ...
%         [tmp_save_dir filesep 'ICA' filesep sprintf('%s_allcond_ICA_TMPEEG.fdt',subj_chars{subj_i})], ...
%         [tmp_save_dir filesep [DEF_EPOCH_PARAMS.gait_trial_chars{:}] filesep sprintf('%s.icatimef',subj_chars{subj_i})]};
%     dirs_del = {[tmp_save_dir filesep 'ICA' filesep 'conn_slide'], ...
%         [tmp_save_dir filesep 'slide_conn_study']};
%     %-
%     for ff = 1:length(fpaths_del)
%         try
%             delete(fpaths_del{ff});
%         catch e
%              fprintf(['error. identifier: %s\n',...
%                      'error. %s\n',...
%                      'error. on subject %s\n',...
%                      'stack. %s\n'],e.identifier,e.message,subj_chars{subj_i},getReport(e));
%         end
%     end
%     %-
%     for ff = 1:length(dirs_del)
%         try
%             tmpd = dir(dirs_del{ff});
%             for fi = 1:length(tmpd)
%                 delete([tmpd(fi).folder filesep tmpd(fi).name]);
%             end
%             rmdir(dirs_del{ff});
%         catch e
%              fprintf(['error. identifier: %s\n',...
%                      'error. %s\n',...
%                      'error. on subject %s\n',...
%                      'stack. %s\n'],e.identifier,e.message,subj_chars{subj_i},getReport(e));
%         end
%     end
%     %-
%     fprintf('%s) deleting data done: %0.2f s\n',subj_chars{subj_i},toc(tt))
% end