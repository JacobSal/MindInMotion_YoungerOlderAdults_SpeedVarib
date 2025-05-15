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
%## R-STATS LOADING
r_stats_dir = [PATHS.src_dir filesep 'r_scripts' filesep 'sbs_lme_mods'];
%--
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
%-
cl_struct = par_load(cluster_k_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
[comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
CLUSTER_PICS = main_cl_inds;
%% (SAVE TABLE) ======================================================== %%
% fext = 'phasec_notw_mw';
% fext = 'phasec_notw_mw_based';
fext = 'phasec_notw_mw_based';
%## LOOP
itc_so = cell(length(CL_STUDY.datasetinfo),1);
parfor subj_i = 1:length(CL_STUDY.datasetinfo)
    try
        tmp_sbs_study = SBS_STUDY;
        tmp_cl_study = CL_STUDY;
        subj_char = tmp_sbs_study.datasetinfo(subj_i).subject;
        subj_ind = strcmp(tmp_cl_study.datasetinfo(subj_i).subject,{tmp_sbs_study.datasetinfo.subject});
        %--
        if ~isempty(subj_ind) && any(subj_ind)
            itc_dato = par_load([tmp_cl_study.datasetinfo(subj_i).filepath filesep sprintf('%s_custom_itc_%s.mat',tmp_cl_study.datasetinfo(subj_i).subject,fext)]);
            %--
            fprintf('%s) Adding ITC structure\n',tmp_cl_study.datasetinfo(subj_i).subject);
            tmp = itc_dato.itc_dat_struct;            
            %## CREATE AVERAGE ERSP MEASURES
            tt = cellfun(@(x) mean(abs(x),3),{tmp.itc_dat_slide},'UniformOutput',false);
            [tmp(:).itc_dat_slide] = tt{:}; %abs(mean(tmp.itc_dat_slide))
            %--
            tt = cellfun(@(x) mean(abs(x),3),{tmp.ersp_dat_slide},'UniformOutput',false);
            [tmp(:).ersp_dat_slide] = tt{:}; %abs(mean(tmp.itc_dat_slide))
            %--
            tt = cellfun(@(x) mean(abs(x),3),{tmp.erspb_dat_slide},'UniformOutput',false);
            [tmp(:).erspb_dat_slide] = tt{:}; %abs(mean(tmp.itc_dat_slide))

            %## RMV FIELDS & STORE
            % tmp = rmfield(tmp,{'itc_freqs','itc_times','ersp_dat_slide','erspb_dat_slide','itc_dat'});
            tmp = rmfield(tmp,{'ersp_dat_slide','erspb_dat_slide','itc_dat'});
            tmp = struct2table(tmp);            
            itc_so{subj_i} = tmp;
            
            %## CREATE AVERAGE ERSP MEASURES
            % i = 1;
            % b_i = 1;
            % itc_times = tmp(i).itc_times;
            % itc_freqs = tmp(i).itc_freqs;
            % dat = tmp(i).ersp_dat_slide(:,:,b_i);
            % dat = tmp(i).erspb_dat_slide(:,:,b_i);
            % dat = mean(abs(tmp(i).itc_dat_slide(:,:,b_i)),3);
            % dat = std(abs(tmp(i).itc_dat_slide(:,:,:)),[],3);
            % dat = tmp(i).itc_dat;
            % figure;
            % contourf(itc_times,itc_freqs,abs(dat),30,'lines','none');
            % set(gca,'fontsize',14);
            % ylabel('Frequency (Hz)');
            % xlabel('Time (ms)');
            % colorbar;
        end
    catch e 
        fprintf('%s) %s\n',tmp_cl_study.datasetinfo(subj_i).subject ,getReport(e));
    end
end
%--
itc_so = util_resolve_table(itc_so);
%-- remove all entries of itc_freqs and itc_times except 1
tmp = num2cell(zeros(size(itc_so,1)-1,1));
itc_so.itc_freqs(2:size(itc_so,1)) = tmp(:);
itc_so.itc_times(2:size(itc_so,1)) = tmp(:);
itc_so.mi_freqs(2:size(itc_so,1),:) = zeros(size(itc_so,1)-1,size(itc_so.mi_freqs(1,:),2)); %tmp(:);
itc_so.mi_phase(2:size(itc_so,1),:) = zeros(size(itc_so,1)-1,size(itc_so.mi_phase(1,:),2));
%--
% itc_so = util_resolve_table(itc_so);
par_save(itc_so,save_dir,sprintf('itc_table_%s.mat',fext));
%(04/30/2025) JS, NH3128 gets an OOM error on the fext='phasec_notw_mw'
%run. Need to include them for the final analysis

%% (SAVE R TABLE) =================================================== %%
% fext = 'phasec_notw_mw';
% fext = 'phasec_notw_mw_based';
fext = 'phasec_notw_mw_based';
CL_NUM_CUTOFF = 13;
GROUP_CHARS = {'H1000','H2000','H3000'};
COND_CHARS = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
FREQ_BOUND = [3,80];
TIME_BOUND = [0,1500];

%## LOOP
itc_so = cell(length(CL_STUDY.datasetinfo),1);
parfor subj_i = 1:length(CL_STUDY.datasetinfo)
    tmp_sbs_study = SBS_STUDY;
    tmp_cl_study = CL_STUDY;
    subj_char = tmp_sbs_study.datasetinfo(subj_i).subject;
    subj_ind = strcmp(tmp_cl_study.datasetinfo(subj_i).subject,{tmp_sbs_study.datasetinfo.subject});
    tmp_tbnd = TIME_BOUND;
    tmp_fbnd = FREQ_BOUND;

    %## TRY
    try      
        if ~isempty(subj_ind) && any(subj_ind)
            itc_dato = par_load([tmp_cl_study.datasetinfo(subj_i).filepath filesep sprintf('%s_custom_itc_%s.mat',tmp_cl_study.datasetinfo(subj_i).subject,fext)]);
            %--
            fprintf('%s) Adding ITC structure\n',tmp_cl_study.datasetinfo(subj_i).subject);
            tmp = itc_dato.itc_dat_struct;
            twp = itc_dato.parameters;
            %--
            alltimes = tmp(1).itc_times';
            allfreqs = tmp(1).itc_freqs';
            tinds = alltimes > tmp_tbnd(1) & alltimes < tmp_tbnd(2);
            finds = allfreqs > tmp_fbnd(1) & allfreqs < tmp_fbnd(2);
            %--
            tclu = unique([tmp.cluster_n]);
            % tcu = unique({tmp.cond_char});
            dats_store = cell(length(tclu)*length(COND_CHARS),1);
            cnt = 1;
            tti = tic();
            for cl_i = 1:length(tclu)
                if tclu(cl_i) <= CL_NUM_CUTOFF
                    for c_i = 1:length(COND_CHARS)
                        inds = [tmp.cluster_n] == tclu(cl_i) & ...
                            strcmp({tmp.cond_char},COND_CHARS{c_i});
                        tt = tmp(inds);
                        %-- convert & take abs of data
                        % rtd = abs(tt.itc_dat_slide(finds,tinds,:));
                        % rtd = abs(tt.erspb_dat_slide(finds,tinds,:));
                        % rtd = abs(tt.ersp_dat_slide(finds,tinds,:));
                        % rtd = abs(tt.itc_dat(finds,tinds));
                        %--
                        rtd = mean(abs(tt.itc_dat_slide(finds,tinds,:)),3);
                        
                        %## STORE DATA
                        gg = strcmp(tt.group_char,GROUP_CHARS);
                        cc = strcmp(tt.cond_char,COND_CHARS);
                        tmp_dats = struct(...
                            'subj_char',tt.subj_char, ...
                            'subj_n',subj_i, ...
                            'group_char',tt.group_char, ...
                            'group_n',find(gg), ...
                            'cond_char',tt.cond_char, ...
                            'cond_n',find(cc), ...
                            'mod_n',tt.mod_n, ...
                            'cluster_n',tt.cluster_n, ...
                            'itc_dat',rtd, ...
                            'itc_freq',tt.itc_freqs(finds), ...
                            'itc_time',tt.itc_times(tinds) ...
                        );
                        dats_store{cnt} = tmp_dats;
                        cnt = cnt + 1;
                    end
                end
            end
            tmp = util_resolve_struct(dats_store);
            itc_so{subj_i} = tmp;
            fprintf('%s) loading and conversion done. %0.2f min',tmp_cl_study.datasetinfo(subj_i).subject,toc(tti)/60)
        end
    catch e 
        fprintf('%s) %s\n',tmp_cl_study.datasetinfo(subj_i).subject ,getReport(e));
    end
end
%--
% itc_so = util_resolve_table(itc_so);
itc_so = util_resolve_struct(itc_so);
%--
tt = num2cell(zeros(size(itc_so,1),1));
[itc_so(2:length(itc_so)).itc_freq] = tt{:}; %abs(mean(tmp.itc_dat_slide))
itc_so.itc_freq
%--
par_save(itc_so,save_dir,sprintf('itc_rdata_struct_%s.mat',fext));
save([save_dir filesep sprintf('itc_rdata_struct_%s.mat',fext)],'itc_so','-mat')

%% ===================================================================== %%
%## VALIDATION PLOT
%{
tmp_plot_struct = struct( ...
    'alltitles',{''},...
    'xlabel','Gait Cycle Time (ms)',...
    'ylabel','Frequency (Hz)',...
    'xticklabel_times',tmp_ntf_struct.timewarpms,...
    'xticklabel_chars',{{'RHS','RTO','LHS','LTO','RHS'}},...
    'xticklabel_angle',45,...
    'clim',double([min(real(itc_dat),[],'all'),max(real(itc_dat),[],'all')]),...
    'font_size',8,...
    'font_name','Arial',...
    'freq_lims',[3,60],...
    'time_lims',[tmp_ntf_struct.timewarpms(1),tmp_ntf_struct.timewarpms(end)],...    
    'contourf_grain',ceil((500/pi())),...
    'alpha_multiple',0.6,...
    'position',[0.1,0.1,0.5,0.5], ...
    'colorbar_shift',[1,1], ...
    'do_display_bandmarks',true);
%--
allersp_mask = ones(size(itc_dat));
allersp_pcond = zeros(size(itc_dat));
% allersp = real(itc_dat); 
allersp = abs(itc_dat);
alltimes = ttimes;
allfreqs = lfreqs;
tmp_plot_struct.clim = double([min(allersp,[],'all'),max(allersp,[],'all')]);
figure;
ax = axes();
[ax] = plot_contourf_cell(ax,allersp,alltimes,allfreqs,...
    allersp_mask,allersp_pcond, ...
    'PLOT_STRUCT',tmp_plot_struct)
%}
%% IMPORT R FLASSO ===================================================== %%
% fext = 'itc_rdata_table_phasec_notw_mw_based';
% fext = sprintf('%s_tables_figs',fext);
% fext = 'itc_rdata_table_phasec_notw_mw_based_flasso_results_bsz5';
fext = 'itc_rdata_table_phasec_notw_mw_based_fl_res_bsz5_nob';
COND_CHARS = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
GROUP_CHARS = {'H1000','H2000','H3000'};
clusters = [3,4,6,8];
% conds = [5,6,7,8];
conds = [1,2,3,4];
%## LOOP
stats_struct = struct('itc_dat',[], ...
    'itc_mod',[], ...
    'lamb_val',[], ...
    'itc_freqs',[], ...
    'itc_times',[], ...
    'subj_char',{''}, ...
    'cond_char',{''}, ...
    'cluster_n',[], ...
    'group_char',[]);
stats_struct = repmat(stats_struct,[1,length(COND_CHARS)*length(GROUP_CHARS)*90]);
cnt = 1;
for cl_i = 1:length(clusters)
    for c_i = 1:length(conds)
        flasso_fpath = [r_stats_dir filesep fext filesep sprintf('allmat_cl%i-c%i.mat',clusters(cl_i),conds(c_i))];
        tmp = load(flasso_fpath,'-mat');
        tmp = tmp.flasso_tbl;
        subjs = tmp.svec;
        for s_i = 1:length(subjs)
            tmp_cl_study = CL_STUDY;
            %--
            itc_freqs = tmp.freqs;
            itc_times = tmp.times;
            %--            
            stats_struct(cnt).subj_char = tmp_cl_study.datasetinfo(subjs(s_i)).subject;
            stats_struct(cnt).cond_char = COND_CHARS{conds(c_i)};
            stats_struct(cnt).group_char = tmp_cl_study.datasetinfo(subjs(s_i)).group;
            stats_struct(cnt).cluster_n = clusters(cl_i);
            stats_struct(cnt).itc_dat = reshape(tmp.odat(:,s_i),length(itc_freqs),length(itc_times));
            stats_struct(cnt).itc_mod = reshape(tmp.mmat(:,s_i),length(itc_freqs),length(itc_times));
            stats_struct(cnt).lamb_val = tmp.lvec(s_i);
            stats_struct(cnt).itc_freqs = itc_freqs;
            stats_struct(cnt).itc_times = itc_times;
            %--     
            % dat = reshape(tmp.mmat(:,s_i),length(itc_freqs),length(itc_times));
            % dat = reshape(tmp.odat(:,s_i),length(itc_freqs),length(itc_times));
            % figure;
            % contourf(itc_times,itc_freqs,dat,30,'lines','none');
            % set(gca,'fontsize',14);
            % ylabel('Frequency (Hz)');
            % xlabel('Time (ms)');
            % colorbar;
            cnt = cnt + 1;
        end
    end
end

%## CONVERT TO TABLE
stats_struct = stats_struct(~cellfun(@isempty,{stats_struct.subj_char}));
tmp = struct2table(stats_struct);
par_save(tmp,save_dir,sprintf('itc_rdata_table_%s_cl%s_c%s.mat',fext,strjoin(string(conds),''),strjoin(string(clusters),'')));
