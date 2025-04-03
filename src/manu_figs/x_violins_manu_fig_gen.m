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
%% (PATHS) ============================================================= %%
%- datset name
DATA_SET = 'MIM_dataset';
%- study name
% STUDY_DNAME = '10172024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej_speed';
% STUDY_DNAME =  '01192025_mim_yaoa_nopowpow_crit_speed';
STUDY_DNAME = '02202025_mim_yaoa_powpow0p3_crit_speed';
STUDY_FNAME = 'kin_eeg_epoch_study';
ANALYSIS_DNAME = 'kin_eeg_step_to_step';
studies_fpath = [PATHS.data_dir filesep DATA_SET filesep '_studies'];

%## CLUSTER LOADING
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
% cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb10'];
% cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb3'];
cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_allcond_rb3'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
cluster_k_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
%## R-STATS LOADING
r_stats_dir = [PATHS.src_dir filesep 'r_scripts' filesep 'sbs_lme_mods'];
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

%% (LOAD STATISTICS & DATA EXCEL SHEET FROM R) ========================= %%
%## FNAMES
%-- r stats
% fextr = 'new_slidingb36';
fextr = 'allcond_perstridefb_mi_nfslidingb36';
% fextr = 'allcond_perstridefb_nfslidingb3';
%## IMPORT DATA
% KIN_TABLE = par_load(save_dir,sprintf('sbs_eeg_psd_%s.mat',fextr));
% %-- r-stats
% RSTATS_IMPORT = readtable([r_stats_dir filesep sprintf('03122025_lme_eeg_kin_%s_stats.xlsx',fextr)], ...
%     "FileType","spreadsheet","UseExcel",true);

%## IMPORT MEANSD DATA
% KIN_TABLE = readtable([r_stats_dir filesep sprintf('03122025_lme_eeg_kin_meansd_%s_tbl.xlsx',fextr)], ...
%     "FileType","spreadsheet","UseExcel",true);
KIN_TABLE = readtable([r_stats_dir filesep sprintf('03312025_lme_eeg_kin_meansd_%s_tbl.xlsx',fextr)], ...
    "FileType","spreadsheet","UseExcel",true);
%-- r-stats
RSTATS_IMPORT = readtable([r_stats_dir filesep sprintf('03312025_lme_eeg_kin_meansd_%s_stats.xlsx',fextr)], ...
    "FileType","spreadsheet","UseExcel",true);
%% MEASURES TO ANALYZE ================================================= %%
%## STATS
try
    STUDY.etc = rmfield(STUDY.etc,'statistics');
end
STUDY = pop_statparams(STUDY,...
    'groupstats','off',...
    'condstats','on',...
    'method','perm',...
    'singletrials','off',...
    'mode','fieldtrip',...
    'fieldtripalpha',NaN,...
    'fieldtripmethod','montecarlo',...
    'fieldtripmcorrect','fdr',...
    'fieldtripnaccu',4000);
stats = STUDY.etc.statistics;
stats.paired{1} = 'on'; % Condition stats
stats.paired{2} = 'off'; % Group stats

%## CLUSTER INFO
%-- 01192025_mim_yaoa_nopowpow_crit_speed (rb3)
cluster_titles = {'Right Cuneus', ...
    'Right Sensorimotor', ...
    'Anterior Cingulate', ...
    'Left Sensorimotor', ...
    'Right Premotor',...
    'Left Posterior Parietal', ...
    'Left Supplementary Motor', ...
    'Right Occipital', ...
    'Mid Cingulate',...
    'Left Temporal',...
    'Left Occipital'};
out = cellfun(@(x) regexp(x,'(.).*\s(...)','tokens'),cluster_titles);
output_titles = cellfun(@(x) strjoin(x,''),out,'UniformOutput',false);
fig_n = 1:length(cluster_titles);
xtick_label_c = {'0.25 m/s','0.50 m/s','0.75 m/s','1.0 m/s'};
cluster_inds_plot = [3,4,5,6,7,8,9,11,12];

%##
speed_xvals = (0:5)*0.25;
c_chars = {'0.25 m/s','0.50 m/s','0.75 m/s','1.0 m/s'};
g_chars_topo = {'Young Adults',{'Older High','Functioning Adults'},{'Older Low','Functioning Adults'}};
g_chars_subp = {'YA','OHFA','OLFA'};
dip_dir = [cluster_k_dir filesep 'topo_dip_inf' filesep 'valid_clusts'];
cmaps_speed = linspecer(4*3);
cmaps_speed = [cmaps_speed(1,:);cmaps_speed(2,:);cmaps_speed(3,:);cmaps_speed(4,:)];
%## EXTRACT PSD DATA
color_dark = cmaps_speed; %color.speed;
color_light = cmaps_speed+0.15; %color.speed_shade;
xtick_label_g = {'0.25','0.50','0.75','1.0'};
x_label = 'speed (m/s)';
cond_offsets = [-0.35,-0.1,0.15,0.40];
%--
des_i = 2;
cl_n = 3;
s_chars = {STUDY.datasetinfo(STUDY.cluster(cl_n).sets).subject};
desdes = cat(1,STUDY.design.variable);
g_chars = {'H1000','H2000','H3000'};
G_ORDER = categorical(g_chars);
%% ===================================================================== %%
%## PARAMETERS
%-
% group_chars = unique(KIN_TABLE.group_char)
designs = unique(KIN_TABLE.model_n);
group_chars = unique(KIN_TABLE.group_char);
cond_chars = unique(KIN_TABLE.cond_char);
clusters = unique(RSTATS_IMPORT.cluster_num);
%--
SAVE_RES = 300;
TITLE_TXT_SIZE = 14;
AX_INIT_X = 0.09;
X_DIM = 2;
%--
TITLE_FONT_SIZE = 14;
TITLE_XSHIFT = 0.4;
TITLE_YSHIFT = 0.975;
TITLE_BOX_SZ = [0.4,0.4];
FIGURE_POSITION =[1,1,6.5,9];
FONT_NAME = 'Arial';
%--
LAB_A_YOFFSET = -0.16;
LAB_A_XOFFSET = -0.125;
LAB_B_YOFFSET = 0.065;
LAB_B_XOFFSET = -0.125;
LAB_C_YOFFSET = 0.06; %0.075
LAB_C_XOFFSET = -0.125;
LAB_D_YOFFSET = 0.09;
LAB_D_XOFFSET = -0.125;
%## 
AXES_DEFAULT_PROPS = {'box','off', ...
    'xtick',[], ...
    'ytick',[],...
    'ztick',[], ...
    'xcolor',[1,1,1], ...
    'ycolor',[1,1,1]};

%## VIOLIN PLOTS
VIO_PLOT_STRUCT = struct('color_map',[],...
    'cond_labels',{{}},...
    'cond_offsets',[-0.35,-0.10,0.15,0.4],...
    'do_group_labels',true, ...
    'group_labels',{{}},...
    'group_offsets',[0.125,0.475,0.812],...
    'group_lab_yoffset',-0.29,...
    'group_lab_fontweight','normal',...
    'group_lab_fontsize',8,...
    'ytick',[], ...
    'ytick_labs',{{''}}, ...
    'y_label',{''},...
    'y_label_fontsize',8,...
    'y_label_fontweight','bold',...
    'ylim',[],...
    'xtick', [], ...
    'xtick_labs', {{''}},...
    'x_label',{''},...
    'x_label_fontsize',8,...
    'x_label_fontweight','bold',...
    'x_label_yoffset',-0.1,...
    'xlim',[],...
    'title',{{''}},...
    'title_fontsize',10,...
    'title_fontweight','normal',...
    'font_size',8,...
    'font_name','Arial',...
    'do_combine_groups',false,...
    'ax_position',[0,0,1,1],...
    'ax_line_width',1,...
    'ax_font_weight','bold', ...
    'xtick_angle',75);
VIO_STRUCT = struct('Width',0.15,...
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
    'do_plot_outlier_marks',true,...
    'use_raw_bandwidth',false);
BRACKET_STRUCT = struct('sig_sign','+',...
    'line_specs',{{'LineStyle','-','LineWidth',2,'Color','k'}},...
    'text_specs',{{'FontSize',8,'FontName','Arial','FontWeight','bold'}},...
    'bracket_conn',[],...
    'conn_offset_y_upper',[],...
    'bracket_offset_y_upper',0,...
    'bracket_offset_y_lower',0,...
    'sig_offset_x',0,...
    'sig_offset_y',[]);
SIGLINE_STRUCT = struct('sig_sign','*',...
    'line_specs',{{'LineStyle','-','LineWidth',2,'Color','k'}},...
    'text_specs',{{'FontSize',8,'FontName','Arial','FontWeight','bold'}},...
    'conn_y',[],...
    'conn_offset_y',[],...
    'sig_offset_x',0,...
    'sig_offset_y',0); 

%## MODELS
COEFF_CHARS_INT = {'(Intercept)','speed_cond_num','group_char1','group_char2', ...
    'speed_cond_num:group_char1','speed_cond_num:group_char2'};
ANV_CHARS_INT = {'(Intercept)','speed_cond_num','group_char','speed_cond_num:group_char'};
ANV_CHARS_GROUP = {'(Intercept)','speed_cond_num','group_char'};
COEFF_CHARS_GROUP = {'(Intercept)','speed_cond_num','group_char1','group_char2'};
%--
%% (ALL SUBJS MODEL) =================================================== %%
meas_ext = 'covistd';
tmp_savedir = [save_dir filesep fextr 'violin_comps'];
mkdir(tmp_savedir);
%##
for cl_i = 1:length(cluster_inds_plot)
    %## INITIATE FIGURE
    %-- initiate params
    cl_ii = find(cluster_inds_plot(cl_i) == double(string(clusters)));
    cl_n = double(string(clusters(cl_ii)));
    atlas_name = cluster_titles{cl_ii};
    
    %## INITIATE FIGURE
    %-- title
    TITLE_FONT_SIZE = 14;
    FIGURE_POSITION =[1,1,6.5,9];
    FONT_NAME = 'Arial';
    TITLE_XSHIFT = 0.4;
    TITLE_YSHIFT = 0.975+fy_shift;
    TITLE_BOX_SZ = [0.4,0.4];
    %-- fig
    fig = figure('color','white', ...
        'Renderer','Painters');    
    annotation('textbox',[TITLE_XSHIFT-(TITLE_BOX_SZ(1)/2/2),TITLE_YSHIFT-(TITLE_BOX_SZ(2)/2), ...
            TITLE_BOX_SZ(1),TITLE_BOX_SZ(2)],...
        'String',atlas_name, ...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle', ...
        'LineStyle','none', ...
        'FontName',FONT_NAME,...
        'FontSize',14, ...
        'FontWeight','Bold', ...
        'Units','normalized');
    set(fig,'Units','inches', ...
        'Position',FIGURE_POSITION, ...
        'PaperUnits','inches', ...
        'PaperSize',[1 1], ...
        'PaperPosition',[0 0 1 1])
    p_sz = get(fig,'Position');
    set(gca,AXES_DEFAULT_PROPS{:});
    hold on;

    %% VIOLIN PLOTS) ================================================== %%
    %--
    YLIM_NTICKS = 5;
    YLIM_SIG_FIGS = 2;
    PRC_YLIM = [3,97];
    YLIM_FAC = 2;
    IM_RESIZE = 0.65;
    AX_W = 0.3;
    AX_H = 0.225;
    AX_FONT_NAME = 'Arial';
    AX_X_SHIFT = 1.4;
    AX_Y_SHIFT = -1.45;
    % F = AX_W*IM_RESIZE;
    % xF = (AX_W*IM_RESIZE*AX_X_SHIFT-F)/2;
    F = AX_W*IM_RESIZE;
    xF = (AX_W*IM_RESIZE*AX_X_SHIFT-F);
    AX_INIT_X = 0.5-(3/2)*F-xF;
    AX_INIT_Y = 0.775;
    X_DIM = 3;
    %-- std measures
    EEG_MEASURES = {'std_avg_theta_fn1', ...
        'std_avg_alpha_fn1', ...
        'std_avg_beta_fn1'};
    EEG_MEASURE_LABS = {'<10*log_{10}(SD(PSD_{N}))>', ...
        '<10*log_{10}(SD(PSD_{N}))>', ...
        '<10*log_{10}(SD(PSD_{N}))>'};
    EEG_MEASURE_TITLES = {'<SD(\theta_{N})>', ...
        '<SD(\alpha_{N})>', ...
        '<SD(\beta_{N})>'};
    %--
    x_shift = AX_INIT_X;   
    y_shift = AX_INIT_Y;
    y_lims = zeros(length(EEG_MEASURES),2);
    ax_s = cell(length(EEG_MEASURES),1);
    x_cnt = 1;    
    for e_i = 1:length(EEG_MEASURES)
        %##
        cond_plot_store = [];
        group_plot_store = [];        
        eeg_measure = EEG_MEASURES{e_i};

        %## SUB-SELECT DATA
        inds = strcmp(KIN_TABLE.model_n,num2str(des_i)) & ...
            KIN_TABLE.cluster_n == cl_n & ...
            ~isnan(KIN_TABLE.(eeg_measure)); %num2str(cl_n);
        tmp_tbl = KIN_TABLE(inds,:); 
        prc_ylim = [round(prctile(tmp_tbl.(eeg_measure),PRC_YLIM(1))-YLIM_FAC*std(tmp_tbl.(eeg_measure)),1),...
                    round(prctile(tmp_tbl.(eeg_measure),PRC_YLIM(2))+YLIM_FAC*std(tmp_tbl.(eeg_measure)),1)];

        %## EXTRACT STATS INFO
        params = [];        
        params.group_chars = {'H1000','H2000','H3000'};
        params.group_order = categorical({'H1000','H2000','H3000'});
        params.model_char_int = 'speed_group_intact_all';
        params.model_char_group = 'speed_group_all';
        params.group_char = 'all';
        %--
        params.anv_chars_int = ANV_CHARS_INT;
        params.anv_chars_group = ANV_CHARS_GROUP;
        params.coeff_chars_int = COEFF_CHARS_INT;
        params.coeff_chars_group = COEFF_CHARS_GROUP;
        %--
        [STATS_STRUCT,CONFINT_STRUCT] = extract_violin_stats(RSTATS_IMPORT,cl_n,eeg_measure,params);

        %## PLOT
        ax = axes();
        %-- set parameters
        % if e_i == length(EEG_MEASURES)
        %     VIO_PLOT_STRUCT.group_labels = g_chars_subp;
        % else
        %     VIO_PLOT_STRUCT.group_labels = {'','',''}; %g_chars_subp;
        % end
        VIO_PLOT_STRUCT.group_labels = {'','',''}; %g_chars_subp;
        VIO_PLOT_STRUCT.color_map = color_dark;
        VIO_PLOT_STRUCT.cond_labels = xtick_label_g;
        VIO_PLOT_STRUCT.title = EEG_MEASURE_TITLES(e_i);
        if e_i == 1
            VIO_PLOT_STRUCT.y_label = EEG_MEASURE_LABS{e_i}; %'10*log_{10}(PSD) - AP. Fit';
        else
            VIO_PLOT_STRUCT.y_label ='';
        end
        VIO_PLOT_STRUCT.ylim = prc_ylim;
        VIO_PLOT_STRUCT.ax_position = [x_shift,y_shift, ...
            AX_W*IM_RESIZE,AX_H*IM_RESIZE];
        VIO_PLOT_STRUCT.x_label = '';
        %-- group violin plot
        ax = group_violin(tmp_tbl,eeg_measure,'speed_n','group_char',...
            ax,...
            'VIOLIN_STRUCT',VIO_STRUCT,...
            'PLOT_STRUCT',VIO_PLOT_STRUCT,...
            'STATS_STRUCT',STATS_STRUCT,...
            'BRACKET_STRUCT',BRACKET_STRUCT,...
            'SIGLINE_STRUCT',SIGLINE_STRUCT, ...
            'CONFINT_STRUCT',CONFINT_STRUCT);
        %-- ax sets
        % if e_i ~= 1
        %     ylabel('');
        % end
        ylabel('');
        y_lims(e_i,:) = get(ax,'YLim');
        ax_s{e_i} = ax;

        %## AX SHIFT
        if x_cnt < X_DIM
            x_shift = x_shift + AX_X_SHIFT*IM_RESIZE*AX_W;
        else
            y_shift = y_shift + AX_Y_SHIFT*IM_RESIZE*AX_H;
            x_shift = AX_INIT_X;
            x_cnt = 0;
        end
        x_cnt = x_cnt + 1; 
    end
    %## Y-LIMIT SETTING
    % YLIM_NTICKS = 5;
    % YLIM_SIG_FIGS = 2;
    %-- ylim
    u = round(max(y_lims(:,2),[],1),1,'significant');
    l = round(min(y_lims(:,1),[],1),1,'significant'); 
    bb = (u-l)/(YLIM_NTICKS-1);
    bbr = round(bb,1,'significant');
    dbb = abs(bbr-bb)*(YLIM_NTICKS-1);
    l = dbb+l;
    %-- yticks    
    tmp = unique(round(linspace(l,u,YLIM_NTICKS), ...
        YLIM_SIG_FIGS,'significant'));
    yticks = tmp;        
    ytick_labs = cellstr(string(tmp));
    for e_i = 1:length(EEG_MEASURES)
        set(ax_s{e_i},'YLim',[l,u], ...
            'YTick',yticks, ...
            'YTickLabel',ytick_labs);        
    end
    hold on;
    %% VIOLIN PLOTS) ================================================== %%
    %-- cov measures
    EEG_MEASURES = {'cov_i_avg_theta_fn1', ...
        'cov_i_avg_alpha_fn1', ...
        'cov_i_avg_beta_fn1'};
    % EEG_MEASURE_LABS = {'log10(<SD(PSD_{\theta,N})>/<<PSD_{\theta,N}>>)', ...
    %     'log10(<SD(PSD_{\alpha,N})>/<<PSD_{\alpha,N}>>)', ...
    %     'log10(<SD(PSD_{\beta,N})>/<<PSD_{\beta,N}>>)'};
    % EEG_MEASURE_LABS = {'log10(<SD(logPSD_{N})>/<<logPSD_{N}>>)', ...
    %     'log10(<SD(logPSD_{N})>/<<logPSD_{N}>>)', ...
    %     'log10(<SD(logPSD_{N})>/<<logPSD_{N}>>)'};
    EEG_MEASURE_LABS = {'log10(SD(logPSD_{N})/<logPSD_{N}>)', ...
        'log10(SD(logPSD_{N})/<logPSD_{N}>)', ...
        'log10(SD(logPSD_{N})/<logPSD_{N}>)'};
    EEG_MEASURE_TITLES = {'<COV(\theta_{N})>', ...
        '<COV(\alpha_{N})>', ...
        '<COV(\beta_{N})>'};
    %--
    % AX_INIT_X = AX_INIT_X + AX_W*IM_RESIZE*AX_X_SHIFT;
    AX_INIT_Y = AX_INIT_Y - AX_H*IM_RESIZE*AX_X_SHIFT;
    x_shift = AX_INIT_X;   
    y_shift = AX_INIT_Y;
    y_lims = zeros(length(EEG_MEASURES),2);
    ax_s = cell(length(EEG_MEASURES),1);
    x_cnt = 1;
    for e_i = 1:length(EEG_MEASURES)
        %##
        cond_plot_store = [];
        group_plot_store = [];        
        eeg_measure = EEG_MEASURES{e_i};

        %## SUB-SELECT DATA        
        inds = strcmp(KIN_TABLE.model_n,num2str(des_i)) & ...
            KIN_TABLE.cluster_n == cl_n & ...
            ~isnan(KIN_TABLE.(eeg_measure)); %num2str(cl_n);
        tmp_tbl = KIN_TABLE(inds,:); 
        prc_ylim = [round(prctile(tmp_tbl.(eeg_measure),PRC_YLIM(1))-YLIM_FAC*std(tmp_tbl.(eeg_measure)),1),...
                    round(prctile(tmp_tbl.(eeg_measure),PRC_YLIM(2))+YLIM_FAC*std(tmp_tbl.(eeg_measure)),1)];

        %## EXTRACT STATS INFO
        params = [];        
        params.group_chars = {'H1000','H2000','H3000'};
        params.group_order = categorical({'H1000','H2000','H3000'});
        params.model_char_int = 'speed_group_intact_all';
        params.model_char_group = 'speed_group_all';
        params.group_char = 'all';
        %--
        params.anv_chars_int = ANV_CHARS_INT;
        params.anv_chars_group = ANV_CHARS_GROUP;
        params.coeff_chars_int = COEFF_CHARS_INT;
        params.coeff_chars_group = COEFF_CHARS_GROUP;
        %--
        [STATS_STRUCT,CONFINT_STRUCT] = extract_violin_stats(RSTATS_IMPORT,cl_n,eeg_measure,params);

        %## PLOT
        ax = axes();
        %-- set parameters
        % if e_i == length(EEG_MEASURES)
        %     VIO_PLOT_STRUCT.group_labels = g_chars_subp;
        % else
        %     VIO_PLOT_STRUCT.group_labels = {'','',''}; %g_chars_subp;
        % end
        VIO_PLOT_STRUCT.group_labels = g_chars_subp;
        VIO_PLOT_STRUCT.color_map = color_dark;
        VIO_PLOT_STRUCT.cond_labels = xtick_label_g;
        VIO_PLOT_STRUCT.title = EEG_MEASURE_TITLES(e_i);
        if e_i == 1
            VIO_PLOT_STRUCT.y_label = EEG_MEASURE_LABS{e_i}; %'10*log_{10}(PSD) - AP. Fit';
        else
            VIO_PLOT_STRUCT.y_label ='';
        end
        VIO_PLOT_STRUCT.ylim = prc_ylim;
        VIO_PLOT_STRUCT.ax_position = [x_shift,y_shift, ...
            AX_W*IM_RESIZE,AX_H*IM_RESIZE];
        VIO_PLOT_STRUCT.x_label = 'Speed (m/s)';
        VIO_PLOT_STRUCT.x_label_yoffset = -0.2;
        %-- group violin plot
        ax = group_violin(tmp_tbl,eeg_measure,'speed_n','group_char',...
            ax,...
            'VIOLIN_STRUCT',VIO_STRUCT,...
            'PLOT_STRUCT',VIO_PLOT_STRUCT,...
            'STATS_STRUCT',STATS_STRUCT,...
            'BRACKET_STRUCT',BRACKET_STRUCT,...
            'SIGLINE_STRUCT',SIGLINE_STRUCT, ...
            'CONFINT_STRUCT',CONFINT_STRUCT);
        %-- ax sets
        y_lims(e_i,:) = get(ax,'YLim');
        ax_s{e_i} = ax;
        % if e_i ~= 1
        %     ylabel('');
        % end          
        %## AX SHIFT
        if x_cnt < X_DIM
            x_shift = x_shift + AX_X_SHIFT*IM_RESIZE*AX_W;
        else
            y_shift = y_shift + AX_Y_SHIFT*IM_RESIZE*AX_H;
            x_shift = AX_INIT_X;
            x_cnt = 0;
        end
        x_cnt = x_cnt + 1; 
    end
    %## Y-LIMIT SETTING
    YLIM_NTICKS = 5;
    YLIM_SIG_FIGS = 2;
    %-- ylim
    u = round(max(y_lims(:,2),[],1),1,'significant');
    l = round(min(y_lims(:,1),[],1),1,'significant'); 
    bb = (u-l)/(YLIM_NTICKS-1);
    bbr = round(bb,1,'significant');
    dbb = abs(bbr-bb)*(YLIM_NTICKS-1);
    l = dbb+l;
    %-- yticks
    tmp = unique(round(linspace(l,u,YLIM_NTICKS), ...
        YLIM_SIG_FIGS,'significant'));
    tmp(tmp < 1e-8 & tmp > -1e-8) = 0;
    yticks = tmp;
    ytick_labs = cellstr(string(tmp));
    for e_i = 1:length(EEG_MEASURES)
        set(ax_s{e_i},'YLim',[l,u], ...
            'YTick',yticks, ...
            'YTickLabel',ytick_labs);        
    end
    hold off;
    %%
    fname = sprintf('cl%s_%s_manscript_plot',string(cl_n),meas_ext);
    exportgraphics(fig,[tmp_savedir filesep fname '.pdf'],...
        'ContentType','Vector')
    exportgraphics(fig,[tmp_savedir filesep fname '.tif'],...
        'Resolution',SAVE_RES)
    % close(fig)
end