%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_speed_kin/.sh

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
ADD_ALL_SUBMODS = false;
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
addpath(SCRIPT_DIR)
addpath(SRC_DIR);
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
% STUDY_DNAME = '10172024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej_speed';
STUDY_DNAME =  '01192025_mim_yaoa_nopowpow_crit_speed';
STUDY_FNAME = 'kin_eeg_epoch_study';
ANALYSIS_DNAME = 'kin_eeg_step_to_step';
%-
cmap_terrain = linspecer(4);
custom_yellow = [254,223,0]/255;
cmap_terrain = [cmap_terrain(3,:);custom_yellow;cmap_terrain(4,:);cmap_terrain(2,:)];
cmap_speed = linspecer(4*3);
cmap_speed = [cmap_speed(1,:);cmap_speed(2,:);cmap_speed(3,:);cmap_speed(4,:)];
%% (PATHS)
studies_fpath = [PATHS.data_dir filesep DATA_SET filesep '_studies'];
%## CLUSTER LOADING
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
% cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb10'];
cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb3'];
% cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb5'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
cluster_k_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];

%## R-STATS LOADING
r_stats_dir = [PATHS.src_dir filesep 'r_scripts' filesep 'sbs_lme_mods'];
%% ===================================================================== %%
% if ~ispc
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME '_UNIX.study'],'filepath',save_dir);
% else
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME '.study'],'filepath',save_dir);
% end
%## LOAD STUDY
if ~ispc
    tmp = load('-mat',[studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep sprintf('%s_UNIX.study',STUDY_FNAME)]);
    STUDY = tmp.STUDY;
else
    tmp = load('-mat',[studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep sprintf('%s.study',STUDY_FNAME)]);
    STUDY = tmp.STUDY;
end
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
%% (LOAD STATISTICS & DATA EXCEL SHEET FROM R) ========================= %%
%## FNAMES
%-- .mat
% fext = 'new_meancondb_nfslidingb6';
% fext = 'new_meandesignb_nfsliding6_stats';
% fext = 'new_sliding6_stats';
%-- r data
fextrd = 'meansd_new_meancondb_nfslidingb6_tbl';
% fextrd = 'meansd_new_meandesignb_nfslidingb6_tbl';
%-- r stats
% fextr = 'new_meancondb_nfslidingb6_stats';
fextr = 'meansd_new_meancondb_nfslidingb6_stats';
% fextr = 'meansd_new_meandesignb_nfslidingb6_stats';

%## IMPORT DATA
% KIN_TABLE = par_load(save_dir,sprintf('sbs_eeg_psd_%s.mat',fext));
%-- r-data
KIN_TABLE = readtable([r_stats_dir filesep sprintf('013022025_lme_eeg_kin_%s.xlsx',fextrd)], ...
    "FileType","spreadsheet","UseExcel",true);
% KIN_TABLE = readtable([r_stats_dir filesep sprintf('013022025_%s.xlsx',fextrd)], ...
%     "FileType","spreadsheet","UseExcel",true);
%-- r-stats
RSTATS_IMPORT = readtable([r_stats_dir filesep sprintf('013022025_lme_eeg_kin_%s.xlsx',fextr)], ...
    "FileType","spreadsheet","UseExcel",true);
X_DIM = 2;

%% MEASURES TO ANALYZE ================================================= %%
%## SLIDING BASELINE MEASURES (RAW)
% EEG_MEASURES = {'mu_avg_theta','std_avg_theta', ...
%     'mu_avg_alpha','std_avg_alpha', ...
%     'mu_avg_beta','std_avg_alpha'};
EEG_MEASURES = {'std_avg_theta_fn1','std_avg_theta_fn2', ...
    'mu_avg_alpha_fn1','mu_avg_alpha_fn2', ...
    'mu_avg_beta_fn1','mu_avg_beta_fn2'};
EEG_MEASURE_LABS = {'Mean 10*log_{10}(PSD_{i})','Std. Dev. 10*log_{10}(PSD_{i})', ...
    'Mean 10*log_{10}(PSD_{i})','Std. Dev. 10*log_{10}(PSD_{i})', ...
    'Mean 10*log_{10}(PSD_{i})','Std. Dev. 10*log_{10}(PSD_{i})'};
EEG_MEASURE_TITLES = {'Sliding Avg. \mu_{\theta}','Sliding Avg. \sigma_{\theta}', ...
    'Sliding Avg. \mu_{\alpha}','Sliding Avg. \sigma_{\alpha}', ...
    'Sliding Avg. \mu_{\beta}','Sliding Avg. \sigma_{\beta}'};

%## MODELS
MODEL_CHARS = {'all_speed_only'};
GROUP_CHARS = {'all'};

%## CLUSTER INFO
%-- 01192025_mim_yaoa_nopowpow_crit_speed (rb3)
cluster_titles = {'','','Right Sensorimotor', ...
    'Precuneus', ...
    'Left Sensorimotor', ...
    'Right Occipital',...
    'Mid Cingulate', ...
    'Left Occipital', ...
    'Left Temporal', ...
    'Left Supplementary Motor',...
    'Right Temporal',...
    'Left Posterior Parietal', ...
    'Right Posterior Parietal'};
cluster_inds_plot = [3,4,5,7,10,12,13];
%% ===================================================================== %%
%## PARAMETERS
%-
% designs = unique(KIN_TABLE.model_n);
% clusters = unique(KIN_TABLE.cluster_n);
% groups = unique(KIN_TABLE.group_n);
% group_chars = unique(KIN_TABLE.group_char);
% cond_chars = unique(KIN_TABLE.cond_char);
speed_ns = unique(KIN_TABLE.speed_n);
clusters = unique(RSTATS_IMPORT.cluster_num);
models = unique(RSTATS_IMPORT.model_char);
%- colors
cmaps_terrain = linspecer(4);
custom_yellow = [254,223,0]/255;
cmaps_terrain = [cmaps_terrain(3,:);custom_yellow;cmaps_terrain(4,:);cmaps_terrain(2,:)];
cmaps_speed = linspecer(4*3);
cmaps_speed = [cmaps_speed(1,:);cmaps_speed(2,:);cmaps_speed(3,:);cmaps_speed(4,:)];
xtick_label_c = {'0.25 m/s','0.50 m/s','0.75 m/s','1.0 m/s'};
%--
SAVE_RES = 300;
TITLE_TXT_SIZE = 14;
IM_RESIZE = 0.9;
AX_W = 0.3;
AX_H = 0.25;
AX_FONT_NAME = 'Arial';
AX_X_SHIFT = 1.5;
AX_Y_SHIFT = 1.4;
AX_INIT_X = 0.1;
AX_INIT_Y = 0.7;
X_DIM = 2;
%--
TITLE_FONT_SIZE = 14;
TITLE_XSHIFT = 0.4;
TITLE_YSHIFT = 0.975;
TITLE_BOX_SZ = [0.4,0.4];
%--
% DO_PLOT_R2 = true;
% REG_TXT_SIZE = 8; % 7
% REG_X_SHIFT = 0.18; % 0.08
% REG_Y_SHIFT = 0.13; % 0.1k
%--
LEG_X_SHIFT = -0.125; %-0.1
LEG_Y_SHIFT =  -0.33; %-0.38
LEG_TXT_SIZE = 9;
LEG_TOKEN_SIZE = 15;
%## 
AXES_DEFAULT_PROPS = {'box','off', ...
    'xtick',[], ...
    'ytick',[],...
    'ztick',[], ...
    'xcolor',[1,1,1], ...
    'ycolor',[1,1,1]};

CL_TITLES = {'','','Right Sensorimotor','Mid Cingulate','Left Temporal','Right Occipital',...
    'Right Premotor','Left Occipital','Left Sensorimotor','Right Posterior Parietal',...
    'Left Posterior Parietal',...
    'Left Supplementary Motor','Right Temporal'};

%## VIOLIN PLOTS
PLOT_STRUCT = struct('color_map',cmaps_speed,...
    'cond_labels',{{'0.25','0.50','0.75','1.0'}},...
    'cond_offsets',[-0.35,-0.1,0.15,0.40],...
    'do_group_labels',false, ...
    'group_labels',{{'All'}},...
    'group_offsets',[0.125,0.475,0.812],...
    'group_lab_yoffset',-0.26,...
    'group_lab_fontweight','normal',...
    'group_lab_fontsize',12,...
    'y_label','',...
    'y_label_fontsize',12,...
    'y_label_fontweight','bold',...
    'ylim',[],...
    'x_label',{''},...
    'x_label_fontsize',12,...
    'x_label_fontweight','bold',...
    'x_label_yoffset',-0.12,...
    'xlim',[],...
    'title',{{''}},...
    'title_fontsize',12,...
    'title_fontweight','normal',...
    'font_size',12,...
    'font_name','Arial',...
    'do_combine_groups',false,...
    'regresslab_txt_size',9,...
    'ax_position',[0,0,1,1],...
    'ax_line_width',1,...
    'xtick_angle',75);
VIOLIN_STRUCT = struct('Width',0.05,...
    'ShowWhiskers',false,...
    'ShowNotches',false,...
    'ShowBox',true,...
    'ShowMedian',true,...
    'Bandwidth',0.1,...
    'QuartileStyle','shadow',...
    'HalfViolin','full',...
    'DataStyle','scatter',...
    'MarkerSize',8,...
    'EdgeColor',[0.5,0.5,0.5],...
    'ViolinAlpha',{{0.2 0.3}},...
    'do_plot_outlier_marks',false,...
    'use_raw_bandwidth',false);
DEF_BRACKET_STRUCT = struct('sig_sign','+',...
    'line_specs',{{'LineStyle','-','LineWidth',2,'Color','k'}},...
    'text_specs',{{'FontSize',10,'FontName','Arial','FontWeight','bold'}},...
    'bracket_conn',[],...
    'conn_offset_y_upper',[],...
    'bracket_offset_y_upper',0,...
    'bracket_offset_y_lower',0,...
    'sig_levels',[0.05,0.01,0.001],...
    'sig_offset_x',0,...
    'sig_offset_y',[]);
DEF_SIGLINE_STRUCT = struct('sig_sign','*',...
    'line_specs',{{'LineStyle','-','LineWidth',2,'Color','k'}},...
    'text_specs',{{'FontSize',10,'FontName','Arial','FontWeight','bold'}},...
    'conn_y',[],...
    'conn_offset_y',[],...
    'sig_levels',[0.05,0.01,0.001],...
    'sig_offset_x',0,...
    'sig_offset_y',0);
%% (ALL SUBJS MODEL) =================================================== %%
%-
m_i = 1;
g_i = 1;
tmp_savedir = [save_dir filesep fext];
mkdir(tmp_savedir);

%##
for c_i = 1:length(clusters)
    %%
    %## INITIATE FIGURE
    fig = figure('color','white');
    set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    annotation('textbox',[TITLE_XSHIFT-(TITLE_BOX_SZ(1)/2/2),TITLE_YSHIFT-(TITLE_BOX_SZ(2)/2),TITLE_BOX_SZ(1),TITLE_BOX_SZ(2)],...
        'String',cluster_titles{double(string(clusters(c_i)))}, ...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle', ...
        'LineStyle','none', ...
        'FontName',AX_FONT_NAME,...
        'FontSize',TITLE_FONT_SIZE, ...
        'FontWeight','Bold', ...
        'Units','normalized');        
    hold on;
    set(gca,AXES_DEFAULT_PROPS{:})
    %-
    x_shift = AX_INIT_X;
    x_cnt = 1;
    y_shift = AX_INIT_Y;
    %##
    vert_shift = 0;
    horiz_shift = 0;
    stats_store = [];
    for e_i = 1:length(EEG_MEASURES)
        %##
        cond_plot_store = [];
        group_plot_store = [];        

        %## SUB-SELECT DATA
        tmp_kint = KIN_TABLE;
        %--
        inds = tmp_kint.model_n == des_i &... 
            tmp_kint.cluster_n == clusters(c_i);        
        %-- get speeds & subjects
        speed_ns = unique(tmp_kint.speed_n);
        subjects = unique(tmp_kint.subj_char);
        %-- ylim calc
        y_lim_calc = [min(tmp_kint.(EEG_MEASURES{e_i}))-1*std(tmp_kint.(EEG_MEASURES{e_i})),...
            max(tmp_kint.(EEG_MEASURES{e_i}))+1*std(tmp_kint.(EEG_MEASURES{e_i}))];
        
        %## ASSIGN STATS
        tmp_stats = RSTATS_IMPORT.cluster_num==double(string(clusters(c_i))) &...
            strcmp(RSTATS_IMPORT.model_char,MODEL_CHARS{m_i}) &...
            strcmp(RSTATS_IMPORT.freq_band_char,EEG_MEASURES{e_i}) &...
            strcmp(RSTATS_IMPORT.kinematic_char,'none') &...
            strcmp(RSTATS_IMPORT.group_char,GROUP_CHARS{g_i});
        tmp_stats = RSTATS_IMPORT(tmp_stats,:);
        tmp_stats = tmp_stats(1,:); %(02/06/2025) JS, bug
        %--
        ran_effs_char = strsplit(tmp_stats.ran_effs_char{1},',');
        ran_effs_n = cellfun(@(x) double(string(x)),strsplit(tmp_stats.ran_effs_n{1},','));
        % 
        %## NORMALIZE DATA USING SUBJECT INTERCEPTS
        subj_chars = unique(tmp_kint.subj_char);
        for subj_i = 1:length(subj_chars)
            ind = strcmp(subj_chars{subj_i},ran_effs_char);
            int = ran_effs_n(ind);
            if ~isempty(int)
                ind = strcmp(subj_chars{subj_i},tmp_kint.subj_char);
                tmp_kint(ind,EEG_MEASURES{e_i}) = tmp_kint(ind,EEG_MEASURES{e_i})-int;
            end
        end
        
        %## EXTRACT STATS
        %- eta/f2
        % F2_sp = tmp_stats.fsq_sp;
        % F2_kin = tmp_stats.fsq_kin;
        fsq_intact = tmp_stats.fsq_intact;

        % F2 = tmp_stats.etasq_sp;
        % F2 = tmp_stats.etasq_kin;
        %- intercepts
        % R2 = tmp_stats.r2_m_int;
        % F2 = tmp_stats.f2_m_int;
        R2 = tmp_stats.r2_c_int;
        % F2 = tmp_stats.f2_c_int;

        %--
        % str = sprintf('eta^2 = %0.2f\ny=(%1.2f)x+(%1.2f)', ...
        %     tmp_stats.etasq_sp,tmp_stats.coeff_sp,tmp_stats.coeff_sp);
        str = sprintf('R^2 = %0.2f\ny=(%1.2f)x+(%1.2f)', ...
            tmp_stats.r2_c_int,tmp_stats.coeff_sp,tmp_stats.coeff_int);
        tmp_stats_struct = struct('anova',{{tmp_stats.pval_sp}},...
            'anova_grp',{{}},...
            'pvals',{{}},...
            'pvals_pairs',{{}},...
            'pvals_grp',{{}},...
            'pvals_grp_pairs',{{}},...
            'regress_pval',{{tmp_stats.pval_sp}},...
            'regress_line',{{[tmp_stats.coeff_int,tmp_stats.coeff_sp]}},...
            'line_type',{'best_fit'},... % ('best_fit' | 'means')
            'regress_xvals',[0,0.25,0.5,0.75,1.0,1.25],...
            'subject_char',[],... % this option when filled prints removal of nan() info
            'group_order',categorical({''}),...
            'display_stats_char',true,... 
            'stats_char',{{str}}, ...
            'stats_char_offsets',[0.5,0]);
        %## PLOT
        ax = axes();
        %--
        PLOT_STRUCT.title = EEG_MEASURE_TITLES(e_i);
        PLOT_STRUCT.y_label = EEG_MEASURE_LABS{e_i}; 
        PLOT_STRUCT.x_label = 'Speed (m/s)';
        % if e_i == 1
        %     PLOT_STRUCT.y_label = EEG_MEASURE_LABS{e_i};
        % else
        %     PLOT_STRUCT.y_label ='';
        % end
        PLOT_STRUCT.ylim = y_lim_calc;
        PLOT_STRUCT.ax_position = [x_shift,y_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE];
        %--
        ax = group_violin(tmp_kint,EEG_MEASURES{e_i},'speed_char','model_char',...
            ax,...
            'VIOLIN_STRUCT',VIOLIN_STRUCT,...
            'PLOT_STRUCT',PLOT_STRUCT,...
            'STATS_STRUCT',tmp_stats_struct,...
            'BRACKET_STRUCT',DEF_BRACKET_STRUCT,...
            'SIGLINE_STRUCT',DEF_SIGLINE_STRUCT);
        if e_i ~= length(EEG_MEASURES)
            xlabel(ax,'');
        end
        %## AX SHIFT
        if x_cnt < X_DIM
            x_shift = x_shift + AX_X_SHIFT*IM_RESIZE*AX_W;
        else
            y_shift = y_shift - AX_Y_SHIFT*IM_RESIZE*AX_H;
            x_shift = AX_INIT_X;
            x_cnt = 0;
        end
        x_cnt = x_cnt + 1;

        %## LEGEND
        % if e_i == length(EEG_MEASURES)
        %     %- lg2                
        %     legend(gca,cond_plot_store);
        %     [lg2,icons,plots,txt]  = legend('boxoff');
        %     tmp = get(lg2,'String');
        %     cnt = 1;
        %     for i = 1:length(cond_plot_store)
        %         tmp{i} = sprintf('%0.2g m/s',double(string(speed_ns(cnt))));
        %         cnt = cnt + 1;
        %     end
        %     set(lg2,'String',tmp,'FontName',AX_FONT_NAME,'FontSize',LEG_TXT_SIZE)
        %     set(lg2,'Orientation','horizontal')
        %     set(lg2,'Units','normalized')
        %     set(lg2,'Position',[AX_INIT_X+LEG_X_SHIFT*IM_RESIZE,...
        %         y_shift+AX_H*IM_RESIZE+LEG_Y_SHIFT*IM_RESIZE,lg2.Position(3),lg2.Position(4)]);
        %     lg2.ItemTokenSize(1) = LEG_TOKEN_SIZE;
        % end
    end
    hold off;
    %##
    exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_eegspeed_plot_allsubj.tiff',string(clusters(c_i)))],...
        'Resolution',300)
    % close(fig)
end
