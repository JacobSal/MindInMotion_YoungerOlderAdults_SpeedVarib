%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/spca_scripts/run_spca_dd_ersp_plots_group.sh

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
ADD_CLEANING_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        SRC_DIR = fileparts(SCRIPT_DIR);
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
    SRC_DIR = fileparts(SCRIPT_DIR);    
end
%## Add Study, Src, & Script Paths
addpath(SRC_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
%% (PARAMETERS) ======================================================== %%
%- file regexps
ICA_REGEXP = '%s_cleanEEG_EMG_HP3std_iCC0p65_iCCEMG0p4_ChanRej0p7_TimeRej0p4_winTol10.set';
%- datset name
DATA_SET = 'MIM_dataset';
%- studies paths
ICA_DNAME = '02212025_YAOAN117_iccR0p65_iccREMG0p4_chanrej_samprej';
% STUDY_DNAME = '10172024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej_speed';
STUDY_DNAME =  '02202025_mim_yaoa_powpow0p3_crit_speed';
SPCA_STUDY_DNAME = '02202025_mim_yaoa_spca_calcs';
STUDY_FNAME_GAIT = 'spca_gait_epoch_study_all';
% STUDY_FNAME_REST = 'spca_rest_slide_study';
ICLABEL_EYE_CUTOFF = 0.75;
%- study group and saving
studies_fpath = [PATHS.data_dir filesep DATA_SET filesep '_studies'];
spca_dir = [studies_fpath filesep sprintf('%s',SPCA_STUDY_DNAME)];
ica_data_dir = [studies_fpath filesep ICA_DNAME]; % JACOB,SAL(02/23/2023)
%- load cluster
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
% cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb5']; % rb10 & rb3 available
% cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb3'];
cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_allcond_rb3']; 
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
cluster_k_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
%-- 
save_dir = [cluster_k_dir filesep 'spca_ersp_plots'];
if ~exist(save_dir,'dir')
    mkdir(save_dir)
end
%% ===================================================================== %%
%## LOAD STUDY
if ~ispc
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
else
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
end
%## LOAD STUDY
%- gait
if ~ispc
    tmp = load('-mat',[spca_dir filesep sprintf('%s_UNIX.study',STUDY_FNAME_GAIT)]);
    SPCA_STUDY = tmp.STUDY;
else
    tmp = load('-mat',[spca_dir filesep sprintf('%s.study',STUDY_FNAME_GAIT)]);
    SPCA_STUDY = tmp.STUDY;
end
%% ADD DESIGNS ========================================================= %%
%-
ERSP_STAT_PARAMS_COND = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','cluster',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
%-
ERSP_STAT_PARAMS_GROUP = ERSP_STAT_PARAMS_COND;
ERSP_STAT_PARAMS_GROUP.groupstats = 'on';
ERSP_STAT_PARAMS_GROUP.condstats = 'off';
%- 
ERSP_STAT_PARAMS_GC = ERSP_STAT_PARAMS_COND;
ERSP_STAT_PARAMS_GC.groupstats = 'on';
%-
ERSP_PARAMS = struct('subbaseline','off',...
    'timerange',[],...
    'ersplim',[-2,2],...
    'freqfac',4,...
    'cycles',[3,0.8],...
    'freqrange',[1,200]);
%--
STUDY_DESI_PARAMS = {{'subjselect',{},...
            'variable2','cond','values2',{'flat','low','med','high'},...
            'variable1','group','values1',{'H1000','H2000','H3000'}},...
            {'subjselect',{},...
            'variable2','cond','values2',{'0p25','0p5','0p75','1p0'},...
            'variable1','group','values1',{'H1000','H2000','H3000'}}};

%## ersp plot per cluster per condition
args = eeglab_struct2args(ERSP_STAT_PARAMS_COND);
STUDY = pop_statparams(STUDY,args{:});
args = eeglab_struct2args(ERSP_PARAMS);
STUDY = pop_erspparams(STUDY,args{:});
STUDY.cache = [];
for des_i = 1:length(STUDY_DESI_PARAMS)
    [STUDY] = std_makedesign(STUDY,[],des_i,STUDY_DESI_PARAMS{des_i}{:});
end

%% ADD CLUSTER INFORMATION 
cl_struct = par_load([cluster_study_fpath filesep sprintf('%i',CLUSTER_K)],sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
[comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);

%% ===================================================================== %%
%## HARD SETS
BOOT_NITERS = 2000;
BOOT_ALPHA = 0.05;
BOOT_CLUST_THRESH = 300;
COLOR_PRCTILE= [15,95];
TERRAIN_REF_CHAR = 'flat';
SPEED_REF_CHAR = '1p0';
SPEED_OVERRIDE_CHARS = {'0.25m/s','0.5m/s','0.75m/s','1.0m/s'};
FREQ_BOUND = [4,60];
TW_CHARS = {'RHS','LTO','LHS','RTO','RHS'};

%##
groups_ind = [1,2,3];
groups = {'YA','OHFA','OLFA'};
group_chars = unique({STUDY.datasetinfo.group});
gait_conds = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
COND_ITERS = {1:4,5:8};
%-- colormaps (not needed)
% COLORS_MAPS_TERRAIN = linspecer(4);
% custom_yellow = [254,223,0]/255;
% COLORS_MAPS_TERRAIN = [COLORS_MAPS_TERRAIN(3,:);custom_yellow;COLORS_MAPS_TERRAIN(4,:);COLORS_MAPS_TERRAIN(2,:)];
% COLOR_MAPS_SPEED = linspecer(4*3);
% COLOR_MAPS_SPEED = [COLOR_MAPS_SPEED(1,:);COLOR_MAPS_SPEED(2,:);COLOR_MAPS_SPEED(3,:);COLOR_MAPS_SPEED(4,:)];

%## GET TIMEWARPING INFORMATION
tmpf = par_load(SPCA_STUDY.datasetinfo(1).filepath,'gait_ersp_spca.mat');
timef_params = tmpf.icatimefopts;
timef_params.timewarpms = tmpf.warptimes;
hardcode_times = tmpf.icatimefopts.times;
hardcode_freqs = tmpf.icatimefopts.freqs;
%-- bounds and crops
time_bound = [timef_params.timewarpms(1),timef_params.timewarpms(end)];
freq_crop = find(hardcode_freqs>=FREQ_BOUND(1) & hardcode_freqs<=FREQ_BOUND(2));
time_crop = find(hardcode_times>=time_bound(1) & hardcode_times<=time_bound(2));
tw_times = tmpf.warptimes;

%## PLOT STRUCT
PLOT_STRUCT = struct('figure_position_inch',[0.5,0.5,6.5,9],...
    'alltitles',{{}},...
    'xlabel','Gait Cycle Events',...
    'ylabel','Frequency (Hz)',...
    'xticklabel_times',tw_times,...
    'xticklabel_chars',{TW_CHARS},...
    'xticklabel_angle',45,...
    'clim',[],...
    'font_size',8,...
    'font_name','Arial',...
    'freq_lims',FREQ_BOUND,...
    'time_lims',time_bound,...
    'stats_title','F Stat (p<0.05)',...
    'figure_title','',...
    'contourf_grain',ceil((500/pi())),...
    'alpha_multiple',0.6,...
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
SAVE_STATS = false;
%%
SPCA_TABLE = par_load(cluster_k_dir,'spca_cluster_table_ersp.mat');
alltimes = hardcode_times(time_crop);
allfreqs = hardcode_freqs(freq_crop);
CLUSTER_PICKS = main_cl_inds;
%% LOOP ============================================================= %%
parfor ii = 1:length(CLUSTER_PICKS)
% for ii = 1:length(CLUSTER_PICKS)
    %## PARAMS
    cl_i = CLUSTER_PICKS(ii);
    tmp_hct = hardcode_times;
    tmp_hcf = hardcode_freqs;
    alltimes = hardcode_times(time_crop);
    allfreqs = tmp_hcf(freq_crop);
    tmp_plot_struct = PLOT_STRUCT;
    tmp_spca_table = SPCA_TABLE;
    tmp_study = STUDY;
    atlas_name = tmp_study.cluster(cl_i).agg_anat_all;
    tmp_plot_struct.figure_title = atlas_name;

    %## INITIATION
    tmp_study = STUDY;
    tmp_spca_study = SPCA_STUDY;

    %%
    STORAGE = cell(length(COND_ITERS),20);
    for des_i = 1:length(COND_ITERS)
        %## REPOP STATS
        %-
        args = eeglab_struct2args(ERSP_STAT_PARAMS_COND);
        % args = eeglab_struct2args(ERSP_STAT_PARAMS_GROUP);
        tmp_study = pop_statparams(tmp_study,args{:});
        %-
%         args = eeglab_struct2args(ERSP_STAT_PARAMS_GC);
%         TMP_STUDY = pop_statparams(TMP_STUDY,args{:});
        %##
        con_con = COND_ITERS{des_i};
        allersp = cell(length(con_con),length(groups_ind)); 
        allgpm = cell(length(con_con),length(groups_ind)); 
        subjs = zeros(length(con_con),length(groups_ind));
        trial_order = cell(length(con_con),length(groups_ind));
        allersp_orig = cell(length(con_con),length(groups_ind));
        allgpm_orig = cell(length(con_con),length(groups_ind));
        spec_ersp = cell(length(con_con),length(groups_ind));
        spec_gpm = cell(length(con_con),length(groups_ind));
        spec_ersp_orig = cell(length(con_con),length(groups_ind));
        spec_gpm_orig = cell(length(con_con),length(groups_ind));
        cnt = 1;
        for group_i = 1:length(groups_ind)
            for cond_i = 1:length(con_con)
                %- get cluster & condition indices
                inds = tmp_spca_table.cluster_c == cl_i & ...
                    strcmp(tmp_spca_table.cond_c,gait_conds{con_con(cond_i)}) & ...
                    groups_ind(group_i) == tmp_spca_table.group_n;
                trial_order{cond_i,group_i} = gait_conds{con_con(cond_i)};
                g_inds = cellfun(@(x) strcmp(x,group_chars{groups_ind(group_i)}),{tmp_study.datasetinfo(tmp_study.cluster(cl_i).sets).group});           
                fprintf('True subjects of group %s in cluster: %i, Alg subjects in cluster: %i\n',group_chars{groups_ind(group_i)},sum(g_inds),sum(inds))
                %- extract corrected ersp
                tmp = cat(3,tmp_spca_table.tf_erspcorr_c{inds});
                tmp = permute(tmp,[2,1,3]);
                allersp{cond_i,group_i} = tmp;
                spec_ersp{cond_i,group_i} = squeeze(mean(tmp(freq_crop,time_crop,:),1));
                %- extract corrected gpm
                tmp = cat(3,tmp_spca_table.tf_gpmcorr_c{inds});
                tmp = permute(tmp,[2,1,3]);
                allgpm{cond_i,group_i} = tmp;
                spec_gpm{cond_i,group_i} = squeeze(mean(tmp(freq_crop,time_crop,:),1));
            end
        end
        %## SUBJECT-SPECIFIC WITHIN CONDITION BASELINE
        [allersp_sb_f,allersp_sb,~,~] = eeglab_baseln(allersp,tmp_hct,tmp_hcf,...
            time_bound,FREQ_BOUND,...
            'DO_COMMON_BASE',false,...
            'DO_SUBJ_BASE',true);
        [allgpm_sb_f,allgpm_sb,~,~] = eeglab_baseln(allgpm,tmp_hct,tmp_hcf,...
            time_bound,FREQ_BOUND,...
            'DO_COMMON_BASE',false,...
            'DO_SUBJ_BASE',true);
        %## COMMON BASELINE
        [allersp_com_f,allersp_com,~,~] = eeglab_baseln(allersp,tmp_hct,tmp_hcf,...
            time_bound,FREQ_BOUND,...
            'DO_COMMON_BASE',true,...
            'DO_SUBJ_BASE',false);
        [allgpm_com_f,allgpm_com,~,~] = eeglab_baseln(allgpm,tmp_hct,tmp_hcf,...
            time_bound,FREQ_BOUND,...
            'DO_COMMON_BASE',true,...
            'DO_SUBJ_BASE',false);
        %## CROP NON-BASELINED FOR PLOTTINGS
        allersp_f = allersp;
        allgpm_f = allgpm;
        allersp = cellfun(@(x) x(freq_crop,time_crop,:),allersp,'uniformoutput',false);
        allgpm = cellfun(@(x) x(freq_crop,time_crop,:),allgpm,'uniformoutput',false);
        %## STORE
        STORAGE{des_i,1} = con_con;
        STORAGE{des_i,2} = subjs;
        STORAGE{des_i,3} = trial_order;
        %-
        STORAGE{des_i,4} = allersp;
        STORAGE{des_i,5} = allgpm;
        %-
        STORAGE{des_i,8} = allersp_f;
        STORAGE{des_i,9} = allgpm_f;
        %-
        STORAGE{des_i,12} = allersp_com;
        STORAGE{des_i,13} = allgpm_com;
        STORAGE{des_i,14} = allersp_com_f;
        STORAGE{des_i,15} = allgpm_com_f;
        %-
        STORAGE{des_i,16} = allersp_sb;
        STORAGE{des_i,17} = allgpm_sb;
        STORAGE{des_i,18} = allersp_sb_f;
        STORAGE{des_i,19} = allgpm_sb_f;
    end
    %% CLIM
    HIGH = 99.52;
    LOW = 100-HIGH;
    %## GPM PLOTS
    INT=5;
    data = [];
    for des_i = 1:length(COND_ITERS)
        tmp = STORAGE{des_i,INT};
        tmp = cellfun(@(x) mean(x,3),tmp,'UniformOutput',false);
        data = cat(3,data,tmp{:});
    end
    bound = max([abs(prctile([data],LOW,'all')),abs(prctile([data],HIGH,'all'))]);
    GPM_CLIM = sort([-round(bound,1),round(bound,1)]);
    %## ERSP PLOTS
    INT=12;
    data = [];
    for des_i = 1:length(COND_ITERS)
        tmp = STORAGE{des_i,INT};
        tmp = cellfun(@(x) mean(x,3),tmp,'UniformOutput',false);
        data = cat(3,data,tmp{:});
    end
    bound = max([abs(prctile([data],LOW,'all')),abs(prctile([data],HIGH,'all'))]);
    ERSP_CLIM = sort([-round(bound,1),round(bound,1)]);
    %%
    for des_i = 1:length(COND_ITERS)
        %## VARS
        con_con = STORAGE{des_i,1};
        subjs = STORAGE{des_i,2};
        trial_order = STORAGE{des_i,3};
        %-
        allersp = STORAGE{des_i,4};
        allgpm = STORAGE{des_i,5};
        % allersp_orig = STORAGE{des_i,6};
        % allgpm_orig = STORAGE{des_i,7};
        %-
        allersp_f  = STORAGE{des_i,8};
        allgpm_f = STORAGE{des_i,9};
        % allersp_orig_f = STORAGE{des_i,10};
        % allgpm_orig_f = STORAGE{des_i,11};
        %-
        allersp_com = STORAGE{des_i,12};
        allgpm_com = STORAGE{des_i,13};
        allersp_com_f = STORAGE{des_i,14};
        allgpm_com_f = STORAGE{des_i,15};
        %-
        allersp_sb = STORAGE{des_i,16};
        allgpm_sb = STORAGE{des_i,17};
        allersp_sb_f = STORAGE{des_i,18};
        allgpm_sb_f = STORAGE{des_i,19};
        %-
        p1 = cell(length(con_con),1);
        p2 = cell(length(con_con),1);
        p3 = cell(length(con_con),1);
        p4 = cell(length(con_con),1);
        %##
        if any(strcmp(trial_order,SPEED_REF_CHAR))
            condnames = SPEED_OVERRIDE_CHARS;
        else
            condnames = trial_order;
        end
        alltitles = std_figtitle('condnames',condnames);
        %% BOOTSTRAPING ALLERSP
        clust_ersp = cell(size(allersp_sb,1),size(allersp_sb,2));
        clust_maskedersp = cell(size(allersp_sb,1),size(allersp_sb,2));
        for group_i = 1:size(allersp_sb,2)
            for cond_i = 1:size(allersp_sb,1)
                fprintf('Performing Stats for Condition %i & Cluster %i\n',cond_i,cl_i);
                tmp = allersp_sb{cond_i,group_i};
                tmp_mean = mean(tmp,3);
                boot_freq = 1:size(tmp,1);
                boot_subj = 1:size(tmp,3);
                boot_surro = zeros(size(tmp,1),size(tmp,2),BOOT_NITERS);
                surro = zeros(size(tmp,1),size(tmp,2),BOOT_NITERS);
                %- scramble time samples and calculate the average across
                %all times and all frequencies and store that value.
                for n = 1:BOOT_NITERS
                    boot_time = randi(size(tmp,2),[size(tmp,2),1]); % random time samples
                    tmpSurro = mean(tmp(boot_freq,boot_time,boot_subj),3);
                    surro(:,:,n) = tmpSurro; % save 2000 iterations of surrogates 
                end
                %- Pull length(subject) surrogate averages from distribution then calc mean across
                %surrogates 
                for n = 1:BOOT_NITERS
                    bootIdx  = randi(BOOT_NITERS,[size(tmp,3),1]);
                    tmpSurro = mean(surro(:,:,bootIdx),3);
                    boot_surro(:,:,n) = tmpSurro;
                end
                pvalMap = stat_surrogate_pvals(boot_surro,tmp_mean,'both');
                pvalMap(pvalMap>1)=1; 
                [p_masked, ~, ~, ~] = fdr_bh(pvalMap,BOOT_ALPHA,'pdep',1);
                % debri removal
                [labelMap,~] = bwlabeln(p_masked);
                tmpDisp = sort(labelMap(:),'descend');
    %             [occurrence,idx] = hist(tmpDisp,unique(tmpDisp));
                [occurrence,idx,~] = histcounts(tmpDisp,unique(tmpDisp));
                kMask = ismember(labelMap,idx((occurrence<BOOT_CLUST_THRESH)));
                finalMask = p_masked-kMask;
                clust_ersp{cond_i,group_i} = tmp_mean; 
                tmp = clust_ersp{cond_i,group_i}; 
                tmp(~finalMask) = 0;
                clust_maskedersp{cond_i,group_i} = tmp;
            end            
        end
        tmp_plot_struct.subplot_width = 0.13;
        tmp_plot_struct.subplot_height = 0.16;
        tmp_plot_struct.subplot_shift_x = 0.035;
        tmp_plot_struct.subplot_shift_y = 0.05;
        tmp_plot_struct.alltitles = alltitles;
        tmp_plot_struct.clim = GPM_CLIM;
        [fig] = plot_txf_mask_contourf(clust_ersp,alltimes,allfreqs,clust_maskedersp,clust_maskedersp,{},...
            'PLOT_STRUCT',tmp_plot_struct);
        drawnow;
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_bs_ersp_sb.tiff',cl_i,des_i)],'Resolution',300);
        % exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_bs_ersp_sb.pdf',cl_i,des_i)], ...
        % 'ContentType','Vector');
        close(fig);
        %% NO BASELINE
        %-- calculate stats
        [pcond_ersp_crop,pgroup_ersp_crop, ~] = ersp_stats_conds(tmp_study,allersp,allfreqs,alltimes);
        [pcond_gpm_crop,pgroup_gpm_crop, ~] = ersp_stats_conds(tmp_study,allgpm,allfreqs,alltimes);
        %-- calculate per condition means
        p1 = cellfun(@(x) mean(x,3),allersp,'UniformOutput',false);
        p2 = cellfun(@(x) mean(x,3),allgpm,'UniformOutput',false);

        %##  ERSP
        tmp_plot_struct.alltitles = alltitles;
        tmp_plot_struct.group_titles = {groups{groups_ind},'Group Stats'};
        tmp_plot_struct.clim = ERSP_CLIM;
        % tmp_plot_struct.subplot_width = 0.13;
        % tmp_plot_struct.subplot_height = 0.16;
        % tmp_plot_struct.subplot_shift_x = 0.035;
        % tmp_plot_struct.subplot_shift_y = 0.05;
        fig = plot_txf_conds_tftopo(p1,alltimes,allfreqs,pcond_ersp_crop,pgroup_ersp_crop,...
            'PLOT_STRUCT',tmp_plot_struct);
        drawnow;
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_groupersp.tiff',cl_i,des_i)],'Resolution',300);
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_groupersp.pdf',cl_i,des_i)],'ContentType','Vector');
        close(fig);

        %## GPM
%         PLOT_STRUCT_PAR.figure_title = 'GPM corrected';
        tmp_plot_struct.alltitles = alltitles;
        tmp_plot_struct.group_titles = {groups{groups_ind},'Group Stats'};
        tmp_plot_struct.clim = GPM_CLIM;
        fig = plot_txf_conds_tftopo(p2,alltimes,allfreqs,pcond_gpm_crop,pgroup_gpm_crop,...
            'PLOT_STRUCT',tmp_plot_struct);
        drawnow;
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_groupgpm.tiff',cl_i,des_i)],'Resolution',300);
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_groupgpm.pdf',cl_i,des_i)],'ContentType','Vector');
        close(fig);
        %% ACROSS CONDITIONS BASELINED ERSPS
        [pcond_ersp_crop, pgroup_ersp_crop, ~] = ersp_stats_conds(tmp_study,allersp_com,allfreqs,alltimes);
        [pcond_gpm_crop, pgroup_gpm_crop, ~] = ersp_stats_conds(tmp_study,allgpm_com,allfreqs,alltimes);
        %- calculate per condition means
        p1 = cellfun(@(x) mean(x,3),allersp_com,'UniformOutput',false);
        p2 = cellfun(@(x) mean(x,3),allgpm_com,'UniformOutput',false);

        %##  ERSP
        tmp_plot_struct.alltitles = alltitles;
        tmp_plot_struct.group_titles = {groups{groups_ind},'Group Stats'};
        tmp_plot_struct.clim = ERSP_CLIM;
        fig = plot_txf_conds_tftopo(p1,alltimes,allfreqs,pcond_ersp_crop,pgroup_ersp_crop,...
            'PLOT_STRUCT',tmp_plot_struct);
        drawnow;
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_groupersp_com.tiff',cl_i,des_i)],'Resolution',300);
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_groupersp_com.pdf',cl_i,des_i)],'ContentType','Vector');
        close(fig);

        %## GPM
        tmp_plot_struct.alltitles = alltitles;
        tmp_plot_struct.group_titles = {groups{groups_ind},'Group Stats'};
        tmp_plot_struct.clim = GPM_CLIM;
        fig = plot_txf_conds_tftopo(p2,alltimes,allfreqs,pcond_gpm_crop,pgroup_gpm_crop,...
            'PLOT_STRUCT',tmp_plot_struct);
        drawnow;
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_groupgpm_com.tiff',cl_i,des_i)],'Resolution',300);
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_groupgpm_com.pdf',cl_i,des_i)],'ContentType','Vector');
        close(fig);
        %% ============================================================= %%
        data_to_proc = {allersp_f,allgpm_f,...
            allersp_com_f,allgpm_com_f};
        data_chars = {'allersp','allgpm',...
            'allersp_com','allgpm_com'};
        data_clim = {ERSP_CLIM,GPM_CLIM,...
            ERSP_CLIM,GPM_CLIM};
        %--
        args = eeglab_struct2args(ERSP_STAT_PARAMS_COND);
%         args = eeglab_struct2args(ERSP_STAT_PARAMS_GC);
        %--
        tmp_diff_study = pop_statparams(tmp_study,args{:});
        for data_i = 1:length(data_to_proc)
            allersp_in = data_to_proc{data_i};
            allersp_out = cell(size(allersp_in,1)-1,size(allersp_in,2));
            ersp_pcond_out = cell(size(allersp_in,1)-1,size(allersp_in,2));
            ersp_masked_out = cell(size(allersp_in,1)-1,size(allersp_in,2));
            plot_alltitles_out = cell(size(allersp_in,1)-1,size(allersp_in,2));
            for group_i = 1:size(allersp_in,2)
                chk_1 = any(strcmp(TERRAIN_REF_CHAR,trial_order),2);
                chk_2 = any(strcmp(SPEED_REF_CHAR,trial_order),2);
                cond_chars = trial_order;
                if any(chk_1)
                    refErspCond = TERRAIN_REF_CHAR;
                    refErspCond_fext = TERRAIN_REF_CHAR;
                    refErspCond_ind = find(chk_1);
                elseif any(chk_2)
                    refErspCond_ind = find(chk_2);
                    refErspCond = SPEED_OVERRIDE_CHARS{refErspCond_ind};
                    refErspCond_fext = SPEED_REF_CHAR;
                else
                    error('Condition for reference ersp not found in STUDY design: %i',des_i)
                end
                inds_to_comp = setdiff(1:length(alltitles),refErspCond_ind);
                if ~isempty(refErspCond)
                    % mask differenec ersps- check that it's sig. different from zer
                    erspDiff = struct('raw',cell(length(inds_to_comp),1),...
                        'masked',cell(length(inds_to_comp),1),....
                        'pcond',cell(length(inds_to_comp),1),...
                        'pgroup',cell(length(inds_to_comp),1));
                    erspDiff_wind = struct('raw',cell(length(inds_to_comp),1),...
                        'masked',cell(length(inds_to_comp),1),...
                        'pcond',cell(length(inds_to_comp),1),...
                        'pgroup',cell(length(inds_to_comp),1));
                    %- calculate pairwise statistics between conditions of interest
                    for c = inds_to_comp
                        %-
                        fprintf('Computing Pair Stat for %s - %s...\n',refErspCond,cond_chars{c})
                        curr_ersp = allersp_in{c,group_i};
                        ref_ersp = allersp_in{refErspCond_ind,group_i};
                        [tmp,tmpg, ~] = ersp_stats_conds(tmp_diff_study,{curr_ersp;ref_ersp},tmp_hcf,tmp_hct);
                        erspDiff(c).raw = mean(curr_ersp-ref_ersp,3);
                        erspDiff(c).masked = erspDiff(c).raw.*tmp{1,1};
                        erspDiff(c).pcond = tmp{1,1};
%                         erspDiff(c).pgroup = tmpg{1,1};
                        %-
                        curr_ersp_wind = allersp_in{c,group_i}(freq_crop,time_crop,:);
                        ref_ersp_wind = allersp_in{refErspCond_ind,group_i}(freq_crop,time_crop,:);
                        [tmp, tmpg, ~] = ersp_stats_conds(tmp_diff_study,{curr_ersp_wind;ref_ersp_wind},allfreqs,alltimes);
                        erspDiff_wind(c).raw = mean(curr_ersp_wind-ref_ersp_wind,3);
                        erspDiff_wind(c).masked = erspDiff_wind(c).raw.*tmp{1,1};
                        erspDiff_wind(c).pcond = tmp{1,1};
%                         erspDiff_wind(c).pgroup = tmpg{1,1};
                    end
                end
                % if SAVE_STATS
                %     if ~exist([save_dir filesep 'stats_out'])
                %         mkdir([save_dir filesep 'stats_out'])
                %     end
                %     par_save(erspDiff_wind,[save_dir filesep 'stats_out'],sprintf('%i_%s_grouperspdiff_wind.mat',cl_i,data_chars{data_i}));
                %     par_save(erspDiff,[save_dir filesep 'stats_out'],sprintf('%i_%s_grouperspdiff.mat',cl_i,data_chars{data_i}));
                % end
                %-
    %             ersp_raw = {erspDiff.raw};
    %             ersp_pcond = {erspDiff.pcond};
%                 ersp_pgroup = {erspDiff.pgroup};
    %             ersp_masked = {erspDiff.masked};
    %             alltimes = hardcode_times;
    %             allfreqs = hardcode_freqs;
                %--
                ersp_raw = {erspDiff_wind.raw}';
                ersp_pcond = {erspDiff_wind.pcond}';
                ersp_pgroup = {erspDiff_wind.pgroup};
                ersp_masked = {erspDiff_wind.masked}';
                alltimes = tmp_hct(time_crop);
                allfreqs = tmp_hcf(freq_crop);
                %--
                ersp_raw = ersp_raw(~cellfun(@isempty,ersp_raw));
                ersp_pcond = ersp_pcond(~cellfun(@isempty,ersp_pcond));
                ersp_masked = ersp_masked(~cellfun(@isempty,ersp_masked));
                plot_alltitles = cell(size(inds_to_comp));
                for j = 1:length(inds_to_comp)
                    cc = inds_to_comp(j);
                    plot_alltitles{j} = sprintf('%s - %s',alltitles{cc},refErspCond);
                end
                allersp_out(:,group_i) = ersp_raw;
                ersp_pcond_out(:,group_i) = ersp_pcond;
                ersp_masked_out(:,group_i) = ersp_masked;
                plot_alltitles_out(:,group_i) = plot_alltitles;
            end

            %## PLOT
            tmp_plot_struct.alltitles = plot_alltitles_out;
            tmp_plot_struct.clim = data_clim{data_i};
            tmp_plot_struct.group_titles = {groups{groups_ind},'Group Stats'};
            [fig] = plot_txf_mask_contourf(allersp_out,alltimes,allfreqs,ersp_masked_out,ersp_pcond_out,{},...
                'PLOT_STRUCT',tmp_plot_struct);
            drawnow;
            exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spcadiff_%s.tiff',cl_i,des_i,data_chars{data_i})], ...
                'Resolution',300);
            exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spcadiff_%s.pdf',cl_i,des_i,data_chars{data_i})], ...
                'ContentType','Vector');
            close(fig);
        end
    end
end
%%

% hardcode_freqs = [3.00000000000000,3.06742258053202,3.13636042918590,3.20684760039063,3.27891891392104,3.35260997209831,3.42795717737705,3.50499775032772,3.58376974802305,3.66431208283782,3.74666454167101,3.83086780560009,3.91696346997695,4.00499406497544,4.09500307660079,4.18703496817112,4.28113520228174,4.37735026326317,4.47572768014374,4.57631605012836,4.67916506260494,4.78432552369029,4.89184938132775,5.00178975094877,5.11420094171129,5.22913848332777,5.34665915349618,5.46682100594745,5.58968339912332,5.71530702549861,5.84375394156257,5.97508759847400,6.10937287340531,6.24667610159107,6.38706510909672,6.53060924632382,6.67737942226828,6.82744813954852,6.98088953022080,7.13777939239961,7.29819522770088,7.46221627952689,7.62992357221147,7.80139995104498,7.97673012319891,8.15600069957009,8.33930023756540,8.52671928484803,8.71835042406688,8.91428831859121,9.11462975927315,9.31947371226118,9.52892136788816,9.74307619065805,9.96204397035612,10.1859328743077,10.4148535008116,10.6489189337742,10.8882447985713,11.1329493191659,11.3831533765093,11.6389805682547,11.9005572698126,12.1680126967792,12.4414789687669,12.7210911746699,13.0069874393964,13.2993089921002,13.5982002359469,13.9038088194464,14.2162857093901,14.5357852654259,14.8624653163106,15.1964872378750,15.5380160327415,15.8872204118333,16.2442728777155,16.6093498098094,16.9826315515215,17.3643024993308,17.7545511938786,18.1535704131050,18.5615572674787,18.9787132973675,19.4052445725961,19.8413617942425,20.2872803987215,20.7432206642077,21.2094078194496,21.6860721550307,22.1734491371293,22.6717795238361,23.1813094840861,23.7022907192622,24.2349805875331,24.7796422309847,25.3365447056091,25.9059631142147,26.4881787423239,27.0834791971242,27.6921585495426,28.3145174795132,28.9508634245091,29.6015107314125,30.2667808117985,30.9470023007079,31.6425112189893,32.3536511392884,33.0807733557695,33.8242370576498,34.5844095066342,35.3616662183387,36.1563911477894,36.9689768790923,37.7998248193646,38.6493453970245,39.5179582645380,40.4060925057219,41.3141868477056,42.2426898776570,43.1920602643787,44.1627669848850,45.1552895560700,46.1701182715835,47.2077544440297,48.2687106526091,49.3535109963264,50.4626913528889,51.5967996434230,52.7563961031407,53.9420535580883,55.1543577081158,56.3939074162048,57.6613150042995,58.9572065557859,60.2822222247693,61.6370165523021,63.0222587897190,64.4386332292388,65.8868395419959,67.3675931236693,68.8816254478788,70.4296844275240,72.0125347842437,73.6309584261788,75.2857548342250,76.9777414569664,78.7077541144847,80.4766474112439,82.2852951582543,84.1345908047237,86.0254478794103,87.9588004418943,89.9356035439920,91.9568337015388,94.0234893767758,96.1365914715780,98.2971838317667,100.506333762756,102.765132556788,105.074696032019,107.436165083718,109.850706247854,112.319512277352,114.843802731298,117.424824577382,120.063852807891,122.762191069532,125.521172307423,128.342159423546,131.226545950008,134.175756737426,137.191248658784,140.274511329112,143.427067841337,146.650475518672,149.946326683911,153.316249446019,156.761908504399,160.285005971229,163.887282212286,167.570516706663,171.336528925812,175.187179232337,179.124369798993,183.150045548333,187.266195113475,191.474851820462,195.778094692702,200.178049477977,204.676889698534,209.276837724781,213.980165873109,218.789197528387,223.706308291685,228.733927153790,233.874537695100,239.130679312479,244.504948473686,250.000000000000];
% hardcode_times = [58,82,106,130,156,180,204,230,254,278,302,328,352,376,402,426,450,474,500,524,548,574,598,622,646,672,696,720,746,770,794,820,844,868,892,918,942,966,992,1016,1040,1064,1090,1114,1138,1164,1188,1212,1236,1262,1286,1310,1336,1360,1384,1410,1434,1458,1482,1508,1532,1556,1582,1606,1630,1654,1680,1704,1728,1754,1778,1802,1826,1852,1876,1900,1926,1950,1974,2000,2024,2048,2072,2098,2122,2146,2172,2196,2220,2244,2270,2294,2318,2344,2368,2392,2416,2442,2466,2490,2516,2540,2564,2588,2614,2638,2662,2688,2712,2736,2762,2786,2810,2834,2860,2884,2908,2934,2958,2982,3006,3032,3056,3080,3106,3130,3154,3178,3204,3228,3252,3278,3302,3326,3352,3376,3400,3424,3450,3474,3498,3524,3548,3572,3596,3622,3646,3670,3696,3720,3744,3768,3794,3818,3842,3868,3892,3916,3942];
%##
% fig = plot_txf(squeeze(PSC1(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
% exportgraphics(fig,[fpath filesep sprintf('chan%i_allcond_psc1_orig.jpg',CHAN_INT)]);
%##
