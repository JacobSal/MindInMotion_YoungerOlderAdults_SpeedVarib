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

%%
fpaths = {STUDY.datasetinfo.filepath};
fextr = 'slidingb3';
dat = par_load(fpaths{1},sprintf('psd_output_%s.mat',fextr));
%--
fooof_freqs = dat.freqs;
% basel_chars = {'slidingb3','slidingb6','slidingb12','slidingb18','slidingb24','slidingb30','slidingb36'};
basel_chars = {'slidingb36'};

% basel_chars = {'slidingb1','slidingb3','slidingb6','slidingb12','slidingb18', ...
%     'slidingb24','slidingb30','slidingb36'};
dat_out_structs = cell(1,length(basel_chars));
%## LOAD
for b_i = 1:length(basel_chars)
    dat_out_structs{b_i} = par_load([cluster_k_dir filesep 'kin_eeg_step_to_step' filesep sprintf('raw_psd_dat_%s.mat',basel_chars{b_i})]);
end
%% (LOAD STATISTICS & DATA EXCEL SHEET FROM R) ========================= %%
%## FNAMES
%-- r stats
fextr = 'new_slidingb36';

%## IMPORT DATA
% KIN_TABLE = par_load(save_dir,sprintf('sbs_eeg_psd_%s.mat',fextr));
% %-- r-stats
% RSTATS_IMPORT = readtable([r_stats_dir filesep sprintf('03122025_lme_eeg_kin_%s_stats.xlsx',fextr)], ...
%     "FileType","spreadsheet","UseExcel",true);

%## IMPORT MEANSD DATA
% KIN_TABLE = readtable([r_stats_dir filesep sprintf('03122025_lme_eeg_kin_meansd_%s_tbl.xlsx',fextr)], ...
%     "FileType","spreadsheet","UseExcel",true);
KIN_TABLE = readtable([r_stats_dir filesep sprintf('02202025_lme_eeg_kin_meansd_%s_tbl.xlsx',fextr)], ...
    "FileType","spreadsheet","UseExcel",true);
%-- r-stats
RSTATS_IMPORT = readtable([r_stats_dir filesep sprintf('03122025_lme_eeg_kin_meansd_%s_stats.xlsx',fextr)], ...
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
    'fieldtripnaccu',2000);
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
violin_ylimits = {{[],[],[],[],[],[],[],[],[],[],[],[]};...
            {[],[],[],[],[],[],[],[],[],[],[],[]};...
            {[],[],[],[],[],[],[],[],[],[],[],[]}};
%--
speed_xvals = (0:5)*0.25;
c_chars = {'0.25 m/s','0.50 m/s','0.75 m/s','1.0 m/s'};
g_chars_topo = {'Young Adults','Older High Functioning Adults','Older Low Functioning Adults'};
g_chars_subp = {'YA','OHFA','OLFA'};
dip_dir = [cluster_k_dir filesep 'topo_dip_inf' filesep 'all'];
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
% c_chars = desdes(strcmp({desdes.label},'cond'));
% c_chars = {c_chars.value};
% g_chars = desdes(strcmp({desdes.label},'group'));
% g_chars = {g_chars.value};
% g_chars = g_chars{1};
% g_chars = {'H3000','H1000','H2000'};
% G_ORDER = categorical(g_chars);
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
% IM_RESIZE = 0.8;
% AX_W = 0.3;
% AX_H = 0.25;
% AX_FONT_NAME = 'Arial';
% AX_X_SHIFT = 1.7;
% AX_Y_SHIFT = -1.4;
% AX_INIT_X = 0.09;
AX_INIT_X = 0.09;
% AX_INIT_Y = 0.7;
% AX_INIT_Y = 0.09; %0.08
X_DIM = 2;
%--
TITLE_FONT_SIZE = 14;
TITLE_XSHIFT = 0.4;
TITLE_YSHIFT = 0.975;
TITLE_BOX_SZ = [0.4,0.4];
FIGURE_POSITION =[1,1,6.5,9];
FONT_NAME = 'Arial';
%--
DIP_IM_DPI = 1000;
AX_INIT_HORIZ_TOPO = 0.085;
AX_INIT_VERT_TOPO = 0.765;
AX_INIT_VERT_DIP = 0.83; %0.79; %0.8 --- %7.2; % inches for some reason, maybe a bug with normal units
AX_INIT_HORIZ_DIP = 0.3846; %2.5; % inches for some reason, maybe a bug with normal units
%--
LAB_A_YOFFSET = -0.16;
LAB_A_XOFFSET = -0.125;
LAB_B_YOFFSET = 0.065;
LAB_B_XOFFSET = -0.125;
LAB_C_YOFFSET = 0.06; %0.075
LAB_C_XOFFSET = -0.125;
LAB_D_YOFFSET = 0.09;
LAB_D_XOFFSET = -0.125;
%--
% DO_PLOT_R2 = true;
% REG_TXT_SIZE = 8; % 7
% REG_X_SHIFT = 0.18; % 0.08
% REG_Y_SHIFT = 0.13; % 0.1k
%--
% LEG_X_SHIFT = -0.125; %-0.1
% LEG_Y_SHIFT =  -0.33; %-0.38
% LEG_TXT_SIZE = 9;
% LEG_TOKEN_SIZE = 15;
%## 
AXES_DEFAULT_PROPS = {'box','off', ...
    'xtick',[], ...
    'ytick',[],...
    'ztick',[], ...
    'xcolor',[1,1,1], ...
    'ycolor',[1,1,1]};

%## VIOLIN PLOTS
PLOT_STRUCT = struct('color_map',[],...
    'cond_labels',{{}},...
    'cond_offsets',[-0.35,-0.10,0.15,0.4],...
    'do_group_labels',true, ...
    'group_labels',{{}},...
    'group_offsets',[0.125,0.475,0.812],...
    'group_lab_yoffset',-0.285,...
    'group_lab_fontweight','normal',...
    'group_lab_fontsize',10,...
    'y_label',{''},...
    'y_label_fontsize',10,...
    'y_label_fontweight','bold',...
    'ylim',[],...
    'x_label',{''},...
    'x_label_fontsize',9,...
    'x_label_fontweight','bold',...
    'x_label_yoffset',-0.1,...
    'xlim',[],...
    'title',{{''}},...
    'title_fontsize',12,...
    'title_fontweight','normal',...
    'font_size',10,...
    'font_name','Arial',...
    'do_combine_groups',false,...
    'ax_position',[0,0,1,1],...
    'ax_line_width',1,...
    'xtick_angle',75);
VIOLIN_STRUCT = struct('Width',0.15,...
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
    'text_specs',{{'FontSize',10,'FontName','Arial','FontWeight','bold'}},...
    'bracket_conn',[],...
    'conn_offset_y_upper',[],...
    'bracket_offset_y_upper',0,...
    'bracket_offset_y_lower',0,...
    'sig_offset_x',0,...
    'sig_offset_y',[]);
SIGLINE_STRUCT = struct('sig_sign','*',...
    'line_specs',{{'LineStyle','-','LineWidth',2,'Color','k'}},...
    'text_specs',{{'FontSize',10,'FontName','Arial','FontWeight','bold'}},...
    'conn_y',[],...
    'conn_offset_y',[],...
    'sig_offset_x',0,...
    'sig_offset_y',0); 

%## SLIDING BASELINE MEASURES (MEANSD)
EEG_MEASURES = {'mu_avg_theta_fn1','std_avg_theta_fn1', ...
    'mu_avg_alpha_fn1','std_avg_alpha_fn1', ...
    'mu_avg_beta_fn1','std_avg_beta_fn1'};
EEG_MEASURE_LABS = {'10*log_{10}(PSD_{N})','10*log_{10}(PSD_{N})', ...
    '10*log_{10}(PSD_{N})','10*log_{10}(PSD_{N})', ...
    '10*log_{10}(PSD_{N})','10*log_{10}(PSD_{N})'};
EEG_MEASURE_TITLES = {'Mean of Mean \theta','Mean of Sliding Std. Dev.  \theta', ...
    'Mean of Mean. \alpha','Mean of Sliding Std. Dev.  \alpha', ...
    'Mean of Mean \beta','Mean of Sliding Std. Dev.  \beta'};
meas_ext = 'std';
violin_ylimits = {{[],[],[],[],[],[],[],[],[],[],[],[]};...
            {[],[],[],[],[],[],[],[],[],[],[],[]};...
            {[],[],[],[],[],[],[],[],[],[],[],[]}};
%--
%## MODELS
% MODEL_CHARS = {'speed_only_all'};
% GROUP_CHARS = {'all'};
% COEFF_CHARS = {'(Intercept)','speed_cond_num'};
% ANV_CHARS_GROUP = {'(Intercept)','speed_cond_num'};
% COEFF_CHARS_GROUP = {'(Intercept)','speed_cond_num'};
%--
% ANV_CHARS_GROUP = {'(Intercept)','speed_cond_num','group_char'};
% COEFF_CHARS_GROUP = {'(Intercept)','speed_cond_num','group_charH2000','group_charH3000'};
% COEFF_CHARS_INT = {'(Intercept)','speed_cond_num','group_charH2000','group_charH3000', ...
%     'speed_cond_num:group_charH2000','speed_cond_num:group_charH3000'};
% ANV_CHARS_INT = {'(Intercept)','speed_cond_num','group_char','speed_cond_num:group_char'};
% ANV_CHARS_GROUP = {'(Intercept)','speed_cond_num','group_char'};
% COEFF_CHARS_GROUP = {'(Intercept)','speed_cond_num','group_charH2000','group_charH3000'};
%--
COEFF_CHARS_INT = {'(Intercept)','speed_cond_num','group_char1','group_char2', ...
    'speed_cond_num:group_char1','speed_cond_num:group_char2'};
ANV_CHARS_INT = {'(Intercept)','speed_cond_num','group_char','speed_cond_num:group_char'};
ANV_CHARS_GROUP = {'(Intercept)','speed_cond_num','group_char'};
COEFF_CHARS_GROUP = {'(Intercept)','speed_cond_num','group_char1','group_char2'};
%--
%% (ALL SUBJS MODEL) =================================================== %%
tmp_savedir = [save_dir filesep fextr];
mkdir(tmp_savedir);
%#%--
% x_shift = AX_INIT_X;
% x_cnt = 1;
% y_shift = AX_INIT_Y;
% vert_shift = 0;
% horiz_shift = 0;
% stats_store = [];
for cl_i = 1:length(cluster_inds_plot)
    %%
    %## INITIATE FIGURE
    %-- initiate params
    cl_ii = find(cluster_inds_plot(cl_i) == double(string(clusters)));
    cl_n = double(string(clusters(cl_ii)));
    atlas_name = cluster_titles{cl_ii};
    
    %## INITIATE FIGURE
    fig = figure('color','white', ...
        'Renderer','Painters');
    TITLE_XSHIFT = 0.4;
    TITLE_YSHIFT = 0.975;
    TITLE_BOX_SZ = [0.4,0.4];
    annotation('textbox',[TITLE_XSHIFT-(TITLE_BOX_SZ(1)/2/2),TITLE_YSHIFT-(TITLE_BOX_SZ(2)/2),TITLE_BOX_SZ(1),TITLE_BOX_SZ(2)],...
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
    
    %## TOPOGRAPHY PLOT
    IM_RESIZE = 0.225;
    ax_position = [AX_INIT_HORIZ_TOPO,AX_INIT_VERT_TOPO,0,0];
    local_plot_topography(fig,STUDY,cl_n, ...
        group_chars,g_chars,g_chars_topo, ...
        ax_position,IM_RESIZE);
    
    %## DIPOLE PLOT
    IM_RESIZE = 1.1;
    dip_fig_path = [dip_dir filesep sprintf('%i_dipplot_alldipspc_top.fig',cl_n)];
    ax_position = [AX_INIT_HORIZ_DIP,AX_INIT_VERT_DIP,0,0];
    label_position = [AX_INIT_X+LAB_A_XOFFSET+(0.1/2),1+LAB_A_YOFFSET+(0.1/2),.1,.1];
    local_plot_dipole_slices(fig,dip_fig_path,IM_RESIZE, ...
        p_sz,ax_position,label_position);

    
    %## EXTRACT PSD DATA =============================================== %%
    %%
    IM_RESIZE = 0.8;
    AX_W = 0.3;
    AX_H = 0.25;
    AX_FONT_NAME = 'Arial';
    AX_X_SHIFT = 1.7;
    AX_Y_SHIFT = -1.4;
    AX_INIT_X = 0.09;
    AX_INIT_Y = 0.7;    
    DESIGNS = {{'flat','low','med','high'},{'0p25','0p5','0p75','1p0'}};
    d_i = 1;
    des_i = 2;
    %--
    ax_store = [];
    psd_avg_char = [];
    y_lim_store = zeros(length(cluster_inds_plot),2);
    ycnt = 1;
    %--
    x_shift = AX_INIT_X;
    x_cnt = 1;
    y_shift = AX_INIT_Y;
    leg_store = [];
    %--
    for g_i = 1:length(g_chars)
        % psd_dat_out1 = dat_out_structs{d_i}.psd_dat;      
        % psd_dat_out2 = dat_out_structs{d_i}.psd_std_dat;
        % psd_dat_out = psd_dat_out2./psd_dat_out1;
        %--
        % psd_dat_out = dat_out_structs{d_i}.psd_dat; 
        % %--
        psd_dat_out = dat_out_structs{d_i}.psd_std_dat;
        if cl_i == 1
            % psd_avg_char = [psd_avg_char,'mean'];
            psd_avg_char = [psd_avg_char,'std'];
            % psd_avg_char = [psd_avg_char,'cov'];
        end
        cond_dat_out = dat_out_structs{d_i}.cond_dat;
        subj_dat_out = dat_out_structs{d_i}.subj_dat;
        group_dat_out = dat_out_structs{d_i}.group_dat;
        %--
        tmp_dat = squeeze(psd_dat_out(:,:,:,cluster_inds_plot(cl_i))); %[subject, frequency, epoch/splice, channel/component]
        tmp_dat = reshape(permute(tmp_dat,[3,1,2]),size(tmp_dat,1)*size(tmp_dat,3),size(tmp_dat,2)); %[subject x epoch/splice, frequency];
        chk = all(~isnan(tmp_dat),2);
        tmp_dat = tmp_dat(chk,:);
        chk = ~all(tmp_dat==0,2);
        tmp_dat = tmp_dat(chk,:);
        % if all(chk)
        %     chk = ~all(tmp_dat==0,2);
        % end
        % tmp_dat = tmp_dat(chk,:);
        % sum(chk)
        %--
        tmp_cond = squeeze(cond_dat_out(:,:,cluster_inds_plot(cl_i))); %[subject, frequency, epoch/splice, channel/component]    
        tmp_cond = reshape(permute(tmp_cond,[2,1]),[size(cond_dat_out,1)*size(cond_dat_out,2),1]);
        chk = cellfun(@isempty,tmp_cond);
        tmp_cond = tmp_cond(~chk);
        conds = unique(tmp_cond);
        %--
        tmp_subj = squeeze(subj_dat_out(:,:,cluster_inds_plot(cl_i))); %[subject, frequency, epoch/splice, channel/component]    
        tmp_subj = reshape(permute(tmp_subj,[2,1]),[size(subj_dat_out,1)*size(subj_dat_out,2),1]);
        % chk = ~all(isnan(tmp_subj),2);
        chk = all(~isnan(tmp_subj),2);
        tmp_subj = tmp_subj(chk,:);
        chk = ~all(tmp_subj==0,2);
        tmp_subj = tmp_subj(chk,:);
        % if all(chk)
        %     chk = ~all(tmp_subj==0,2);
        % end
        % tmp_subj = tmp_subj(chk,:);
        subjs = unique(tmp_subj);
        %--
        tmp_group = squeeze(group_dat_out(:,:,cluster_inds_plot(cl_i))); %[subject, frequency, epoch/splice, channel/component]    
        tmp_group = reshape(permute(tmp_group,[2,1]),[size(group_dat_out,1)*size(group_dat_out,2),1]);
        chk = cellfun(@isempty,tmp_group);
        tmp_group = tmp_group(~chk);
        %--  get good indices
        conds = unique(tmp_cond);        
        group_chars = unique(tmp_group);
        indsc = cellfun(@(x) any(strcmp(x,DESIGNS{des_i})),conds);
        conds = conds(indsc);
        indsc = cellfun(@(x) any(strcmp(x,conds)),tmp_cond);        
        indsg = strcmp(tmp_group,group_chars{g_i});
        indsg = indsg & indsc;
        %--
        tmp_cond = tmp_cond(indsg,:);
        tmp_subj = tmp_subj(indsg,:);
        tmp_dat = tmp_dat(indsg,:);
        % tmp_group = tmp_cond(indsg,:);
        %## STATISTICS
        %-- Ho : all samples come from the same distribution
        %-- Ha : all samples come from different distributions
        tmp_psd_in = cell(length(conds),1);
        for c_ii = 1:length(conds)
            indc = strcmp(tmp_cond,conds{c_ii});
            tmp = nan(size(tmp_dat,2),length(subjs));
            for s_i = 1:length(subjs)
                inds = tmp_subj == subjs(s_i);
                chk = indc & inds;
                tmp(:,s_i) = mean(tmp_dat(chk,:),1);
                % tmp(:,s_i) = std(tmp_dat(chk,:),[],1);        
                % tmp(:,s_i) = prctile(tmp_dat(chk,:),75,1) - prctile(tmp_dat(chk,:),25,1);
                if c_ii == 1 && s_i == 1 && cl_i == 1
                    psd_avg_char = [psd_avg_char,'mean'];
                    % psd_avg_char = [psd_avg_char,'std'];
                    % psd_avg_char = [psd_avg_char,'prct'];
                end
            end
            % tmp = tmp(:,all(tmp ~= 0,1));
            tmp = tmp(:,all(~isnan(tmp),1));
            tmp_psd_in{c_ii,1} = tmp; %tmp_dat(inds,:);
        end
        [pcond, pgroup, pinter, ~, ~, ~] = ...
            std_stat(tmp_psd_in, stats);    
        pcond=pcond{1} < 0.05;
        %## PLOT
        ax = axes();
        for c_i = 1:length(conds)
            %--
            % ind = find(cluster_inds_plot(cl_i) == cluster_inds_plot);
            PLOT_STRUCT.title = {sprintf('%s',cluster_titles{cluster_inds_plot(cl_i)})};
            PLOT_STRUCT.ax_position = [x_shift,y_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE];
            PLOT_STRUCT.xlim = [3,40];
            % PLOT_STRUCT.ylim = [-2.5,5]; %sort([prctile([tmp_psd_in{:}],99,'all'),prctile([tmp_psd_in{:}],1,'all')]);
            mu = mean(cat(2,tmp_psd_in{:}),[2,1]);
            sd = std(cat(2,tmp_psd_in{:}),[],[2,1]);
            y_lim_store(ycnt,:) = [mu-1.75*sd,mu+1.75*sd];
            PLOT_STRUCT.ylim = y_lim_store(ycnt,:);
            ycnt = ycnt+1;
           
            disp(PLOT_STRUCT.ylim);
            %--
            LINE_STRUCT.line_avg_fcn = @(x) mean(x,2);
            LINE_STRUCT.line_color = cmaps_speed(c_i,:);
            LINE_STRUCT.line_alpha = 0.7;
            LINE_STRUCT.line_label = xtick_label_c{c_i};
            %--
            LINE_STRUCT.do_err_shading = true;
            LINE_STRUCT.err_color = cmaps_speed(c_i,:)+0.15;
            LINE_STRUCT.err_alpha = 0.3;
            LINE_STRUCT.err_upr_bnd_fcn = @(x) mean(x,2) + std(x,[],2);
            LINE_STRUCT.err_lwr_bnd_fcn = @(x) mean(x,2) - std(x,[],2);
            %--
            [ax,Pa,Li] = plot_psd(ax,tmp_psd_in{c_i,1},fooof_freqs, ...
                'LINE_STRUCT',LINE_STRUCT, ...
                'PLOT_STRUCT',PLOT_STRUCT);
            if c_i == 1
                ax_store = [ax_store, ax];
            end
            if cl_i == 1
                leg_store = [leg_store, Li];
            end
            hold on;
        end        
        [axsignif,Pa] = plot_psd_stats(ax,fooof_freqs,pcond, ...
            'background','Frequency (Hz)');
        if cl_i < length(cluster_inds_plot)-1
            xlabel('');
        end
        %--
        y_lim_store = [min(y_lim_store,[],'all'),max(y_lim_store,[],'all')];
        for aa = 1:length(ax_store)
            set(ax_store(aa),'YLim',y_lim_store)
        end
        
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

    %## LEGEND
    % %- lg2                
    % legend(gca,leg_store);
    % [lg2,icons,plots,txt]  = legend('boxoff');
    % tmp = get(lg2,'String');
    % cnt = 1;
    % for i = 1:length(leg_store)
    %     tmp{i} = sprintf('%s',leg_chars{cnt});
    %     cnt = cnt + 1;
    % end
    % set(lg2,'String',tmp,'FontName',AX_FONT_NAME,'FontSize',LEG_TXT_SIZE)
    % set(lg2,'Orientation','horizontal')
    % set(lg2,'Units','normalized')
    % set(lg2,'Position',[AX_INIT_X+LEG_X_SHIFT*IM_RESIZE*AX_W,...
    %     y_shift+AX_H*IM_RESIZE+LEG_Y_SHIFT*IM_RESIZE*AX_H,lg2.Position(3),lg2.Position(4)]);
    % lg2.ItemTokenSize(1) = LEG_TOKEN_SIZE;
    % hold off;

    

    %## VIOLIN PLOTS) ================================================== %%
    %-
    x_shift = AX_INIT_X;
    x_cnt = 1;
    y_shift = AX_INIT_Y;
    for e_i = 1:length(EEG_MEASURES)
        %##
        cond_plot_store = [];
        group_plot_store = [];        

        %## SUB-SELECT DATA
        tmp_kint = KIN_TABLE;
        %-- select cluster
        % inds = tmp_kint.model_n == des_i &... 
        %     tmp_kint.cluster_n == clusters(c_i); 
        inds = tmp_kint.cluster_n == clusters(cl_ii); 
        %-- get speeds & subjects
        % speed_ns = unique(tmp_kint.speed_n);
        subjects = unique(tmp_kint.subj_char);   
        %-- subset kintable
        tmp_kint = tmp_kint(inds,:);
        %--
        % cnt = 1;
        % y_shift = 0;    
        % x_shift= 0;
        % prc_ylim = zeros(length(EEG_MEASURES),2);
        % inds = strcmp(KIN_TABLE.model_n,num2str(des_i)) & KIN_TABLE.cluster_n == num2str(cl_n);
        % tmp_tbl = KIN_TABLE(inds,:); 
        % for m_i = 1:length(EEG_MEASURES)
        %     if isempty(violin_ylimits{m_i}{cl_ii})
        %         prc_ylim(m_i,:) = [floor(prctile(tmp_tbl.(EEG_MEASURES{m_i}),1))-floor(std(tmp_tbl.(EEG_MEASURES{m_i}))),...
        %             ceil(prctile(tmp_tbl.(EEG_MEASURES{m_i}),99))+ceil(std(tmp_tbl.(EEG_MEASURES{m_i})))*1.5];
        %         disp(prc_ylim)
        %     else
        %         prc_ylim(m_i,:) = violin_ylimits{m_i}{cl_ii};
        %     end
        % end
        %--
        inds = strcmp(KIN_TABLE.model_n,num2str(des_i)) & KIN_TABLE.cluster_n == cl_n; %num2str(cl_n);
        tmp_tbl = KIN_TABLE(inds,:); 
        prc_ylim = [floor(prctile(tmp_tbl.(EEG_MEASURES{e_i}),1))-floor(std(tmp_tbl.(EEG_MEASURES{e_i}))),...
                    ceil(prctile(tmp_tbl.(EEG_MEASURES{e_i}),99))+ceil(std(tmp_tbl.(EEG_MEASURES{e_i})))*1.5];
               
        %## EXTRACT STATS INFO        
        tmp_stats = RSTATS_IMPORT.cluster_num==double(string(cl_n)) &...
            strcmp(RSTATS_IMPORT.model_char,'speed_group_intact_all') &...
            strcmp(RSTATS_IMPORT.freq_band_char,EEG_MEASURES{e_i}) &...
            strcmp(RSTATS_IMPORT.group_char,'all');
        tmp_stats = RSTATS_IMPORT(tmp_stats,:);
        %--
        tmp_ac = strsplit(tmp_stats.anv_chars{1},',');
        tmp_anv = cellfun(@(x) double(string(x)),strsplit(tmp_stats.anv_pvals{1},','));
        %-- anova p-values
        anvs = zeros(length(ANV_CHARS_INT),1);
        for cc = 1:length(ANV_CHARS_INT)
            ind = strcmp(ANV_CHARS_INT{cc},tmp_ac);
            if ~isempty(ind)
                anvs(cc) = tmp_anv(ind);              
            else
                fprintf("Coefficient %s not found.\n",ANV_CHARS_INT{cc})
            end            
        end
        
        %## CHECK FOR SIGNIFICANT INTERACTION
        if anvs(4) > 0.05
            %## USE NON-INTERACTION MODEL
            tmp_stats = RSTATS_IMPORT.cluster_num==double(string(cl_n)) &...
                strcmp(RSTATS_IMPORT.model_char,'speed_group_all') &...
                strcmp(RSTATS_IMPORT.freq_band_char,EEG_MEASURES{e_i}) &...
                strcmp(RSTATS_IMPORT.group_char,'all');
            tmp_stats = RSTATS_IMPORT(tmp_stats,:);
            %--
            tmp_ac = strsplit(tmp_stats.anv_chars{1},',');
            tmp_anv = cellfun(@(x) double(string(x)),strsplit(tmp_stats.anv_pvals{1},','));
            %-- anova p-values
            anvs = zeros(length(ANV_CHARS_GROUP),1);
            for cc = 1:length(ANV_CHARS_GROUP)
                ind = strcmp(ANV_CHARS_GROUP{cc},tmp_ac);
                if ~isempty(ind)
                    anvs(cc) = tmp_anv(ind);              
                else
                    fprintf("Coefficient %s not found.\n",ANV_CHARS_GROUP{cc})
                end            
            end
            tmp_cc = strsplit(tmp_stats.coeff_chars{1},',');
            tmp_fsq_chars = strsplit(tmp_stats.fsq_chars{1},',');
            tmp_ci_chars = strsplit(tmp_stats.confint_chars{1},',');
            tmp_coeffs = cellfun(@(x) double(string(x)),strsplit(tmp_stats.coeffs{1},','));
            tmp_fsq = cellfun(@(x) double(string(x)),strsplit(tmp_stats.fsq_vals{1},','));
            % tmp_em = cellfun(@(x) double(string(x)),strsplit(tmp_stats.emmeans{1},','));
            tmp_ci_lwr = cellfun(@(x) double(string(x)),strsplit(tmp_stats.confint_lwr{1},','));
            tmp_ci_upr = cellfun(@(x) double(string(x)),strsplit(tmp_stats.confint_upr{1},','));
            %-- model coefficients       
            coeffs = zeros(length(COEFF_CHARS_GROUP),1);        
            for cc = 1:length(COEFF_CHARS_GROUP)
                ind = strcmp(COEFF_CHARS_GROUP{cc},tmp_cc);
                if ~isempty(ind)
                    coeffs(cc) = tmp_coeffs(ind);              
                else
                    fprintf("Coefficient %s not found.\n",COEFF_CHARS_GROUP{cc})
                end            
            end        
            %-- cohens f^2 values
            fsq_chars = strcmp(ANV_CHARS_GROUP,'(Intercept)');
            fsq_chars = ANV_CHARS_GROUP(~fsq_chars);
            fsqs = zeros(length(fsq_chars),1);
            for cc = 1:length(fsq_chars)
                ind = strcmp(fsq_chars{cc},tmp_fsq_chars);
                if ~isempty(ind)
                    fsqs(cc) = tmp_fsq(ind);
                else
                    fprintf("Coefficient %s not found.\n",ANV_CHARS_GROUP{cc})
                end
            end
            %-- confidence intervals
            ci_chars = strcmp(g_chars,'(Intercept)');
            ci_chars = g_chars(~ci_chars);
            cis = zeros(length(ci_chars),1,2);
            for cc = 1:length(ci_chars)
                ind = strcmp(ci_chars{cc},tmp_ci_chars);
                if ~isempty(ind)
                    cis(cc,1,:) = [tmp_ci_lwr(ind),tmp_ci_upr(ind)];
                else
                    fprintf("Coefficient %s not found.\n",g_chars{cc})
                end
            end
            %--
            % ran_effs_char = strsplit(tmp_stats.ran_effs_char{1},',');
            % ran_effs_n = cellfun(@(x) double(string(x)),strsplit(tmp_stats.ran_effs_n{1},','));
            %--
            if anvs(3) < 0.05 && anvs(3) > 0.01
                strg = '^{+}';
            elseif anvs(3) <= 0.01 && anvs(3) > 0.001
                strg = '^{++}';
            elseif anvs(3) <= 0.001
                strg = '^{+++}';
            else
                strg = '^{ns}';
            end
            %--
            if anvs(2) < 0.05 && anvs(2) > 0.01
                strs = '^{*}';
            elseif anvs(2) <= 0.01 && anvs(2) > 0.001
                strs = '^{**}';
            elseif anvs(2) <= 0.001
                strs = '^{***}';
            else
                strs = '^{ns}';
            end

            %## ASSIGN STATS
            str = {[sprintf('%sf_{s}^{2}=%1.2f    %sf_{g}^{2}=%1.2f\nR^2=%1.2f', ...
                strs,fsqs(1), ...
                strg,fsqs(2), ...
                tmp_stats.r2_c_int)],'',''};
            %--
            if anvs(3) > 0.05 && anvs(2) > 0.05
                chkd = false;
            else
                chkd = true;
            end
            CONFINT_STRUCT = struct('do_display',chkd, ...
                    'y_bnds',cis, ...
                    'x_vals',repmat(0,[3,1]), ...
                    'errbar_struct',struct('line_specs',{{'LineStyle','-','LineWidth',2,'Color','k'}}, ...
                        'err_bar_width',0.25));
            % regl =[[coeffs(1),coeffs(2)]; ...
            %         [coeffs(1)+coeffs(3),coeffs(2)]; ...
            %         [coeffs(1)+coeffs(4),coeffs(2)]];
            regl = [[coeffs(1)+coeffs(3),coeffs(2)]; ...
                    [coeffs(1)+coeffs(4),coeffs(2)]; ...
                    [coeffs(1),coeffs(2)]];
            tssc = struct('var_type','continuous', ...
                'anova',anvs(2), ...
                'multc_pvals',[],...
                'multc_pvals_pairs',[],...
                'regress_line',regl, ... 
                'line_type',{'best_fit'},... % ('best_fit' | 'means') % continuous predictors
                'regress_xvals',(0:5)*0.25,... % continuous predictors
                'order',{{}});
            tssg = struct('var_type','categorical', ...
                'anova',anvs(3), ...
                'multc_pvals',[],...
                'multc_pvals_pairs',[],...
                'regress_line',[],... 
                'line_type',{'best_fit'},... % ('best_fit' | 'means') % continuous predictors
                'regress_xvals',0,... % continuous predictors
                'order',{G_ORDER});
            STATS_STRUCT = struct('sig_levels',[0.05,0.01,0.001],...
                'cond_stats',tssc, ...
                'group_stats',tssg, ...
                'stats_char',struct('str',{str}, ...
                    'offsets',[-0.1,-0.05], ...
                    'font_size',9, ...
                    'do_display',true));
        else
            %## USE INTERACTION MODEL
            tmp_cc = strsplit(tmp_stats.coeff_chars{1},',');
            tmp_fsq_chars = strsplit(tmp_stats.fsq_chars{1},',');
            tmp_ci_chars = strsplit(tmp_stats.confint_chars{1},',');
            tmp_coeffs = cellfun(@(x) double(string(x)),strsplit(tmp_stats.coeffs{1},','));
            tmp_fsq = cellfun(@(x) double(string(x)),strsplit(tmp_stats.fsq_vals{1},','));
            % tmp_em = cellfun(@(x) double(string(x)),strsplit(tmp_stats.emmeans{1},','));
            tmp_ci_lwr = cellfun(@(x) double(string(x)),strsplit(tmp_stats.confint_lwr{1},','));
            tmp_ci_upr = cellfun(@(x) double(string(x)),strsplit(tmp_stats.confint_upr{1},','));
            %-- model coefficients       
            coeffs = zeros(length(COEFF_CHARS_INT),1);        
            for cc = 1:length(COEFF_CHARS_INT)
                ind = strcmp(COEFF_CHARS_INT{cc},tmp_cc);
                if ~isempty(ind)
                    coeffs(cc) = tmp_coeffs(ind);              
                else
                    fprintf("Coefficient %s not found.\n",COEFF_CHARS_INT{cc})
                end            
            end        
            %-- cohens f^2 values
            fsq_chars = strcmp(ANV_CHARS_INT,'(Intercept)');
            fsq_chars = ANV_CHARS_INT(~fsq_chars);
            fsqs = zeros(length(fsq_chars),1);
            for cc = 1:length(fsq_chars)
                ind = strcmp(fsq_chars{cc},tmp_fsq_chars);
                if ~isempty(ind)
                    fsqs(cc) = tmp_fsq(ind);
                else
                    fprintf("Coefficient %s not found.\n",ANV_CHARS_INT{cc})
                end
            end
            %-- confidence intervals
            ci_chars = strcmp(g_chars,'(Intercept)');
            ci_chars = g_chars(~ci_chars);
            cis = zeros(length(ci_chars),1,2);
            for cc = 1:length(ci_chars)
                ind = strcmp(ci_chars{cc},tmp_ci_chars);
                if ~isempty(ind)
                    cis(cc,1,:) = [tmp_ci_lwr(ind),tmp_ci_upr(ind)];
                else
                    fprintf("Coefficient %s not found.\n",g_chars{cc})
                end
            end
            %--
            if anvs(4) < 0.05 && anvs(4) > 0.01
                stri = '*';
            elseif anvs(4) <= 0.01 && anvs(4) > 0.001
                stri = '**';
            elseif anvs(4) <= 0.001
                stri = '***';
            else
                stri = '^{ns}';
            end
            %## ASSIGN STATS
            str = {sprintf('%sf_{s:g}^{2}=%1.2f\nR^2=%1.2f', ...
                stri,fsqs(3),tmp_stats.r2_c_int),'',''};
            txt_sz = 9;
            offs = [-0.1,-0.05];
            %--
            % str = {sprintf('%sm_{ya}=%1.2f  m_{ohf}=%1.2f  m_{olf}=%1.2f\nR^2=%1.2f', ...
            %     stri,coeffs(2), ...
            %     coeffs(2)+coeffs(5), ...
            %     coeffs(2)+coeffs(6), ...
            %     tmp_stats.r2_c_int),'',''};
            % txt_sz = 9;
            % offs = [-0.11,-0.05];
            %--
            % str = {sprintf('%sy=(%1.1f)x+(%1.1f)\nR^2=%1.2f',stri,coeffs(2),coeffs(1),tmp_stats.r2_c_int), ...
            %     sprintf('y=(%1.1f)x+(%1.1f)',coeffs(2)+coeffs(5),coeffs(1)+coeffs(3)), ...
            %     sprintf('y=(%1.1f)x+(%1.1f)',coeffs(2)+coeffs(6),coeffs(1)+coeffs(4))};
            % txt_sz = 7;
            % offs = [-0.11,-0.05];
            %--
            if anvs(4) > 0.05 && anvs(3) > 0.05
                chkd = false;
            else
                chkd = true;
            end
            CONFINT_STRUCT = struct('do_display',chkd, ...
                    'y_bnds',cis, ...
                    'x_vals',repmat(0,[3,1]), ...
                    'errbar_struct',struct('line_specs',{{'LineStyle','-','LineWidth',2,'Color','k'}}, ...
                        'err_bar_width',0.25));
            % regl = [[coeffs(1)+coeffs(3),coeffs(2)]; ...
            %         [coeffs(1)+coeffs(4),coeffs(2)]; ...
            %         [coeffs(1),coeffs(2)]];
            regl = [[coeffs(1)+coeffs(3),coeffs(2)+coeffs(5)]; ...
                    [coeffs(1)+coeffs(4),coeffs(2)+coeffs(6)]; ...
                    [coeffs(1),coeffs(2)]];
            tssc = struct('var_type','continuous', ...
                'anova',anvs(4), ...
                'multc_pvals',[],...
                'multc_pvals_pairs',[],...
                'regress_line',regl, ...
                'line_type',{'best_fit'},... % ('best_fit' | 'means') % continuous predictors
                'regress_xvals',(0:5)*0.25,... % continuous predictors
                'order',{{}});
            tssg = struct('var_type','categorical', ...
                'anova',anvs(4), ...
                'multc_pvals',[],...
                'multc_pvals_pairs',[],...
                'regress_line',[],... 
                'line_type',{'best_fit'},... % ('best_fit' | 'means') % continuous predictors
                'regress_xvals',0,... % continuous predictors
                'order',{G_ORDER});
            STATS_STRUCT = struct('sig_levels',[0.05,0.01,0.001],...
                'cond_stats',tssc, ...
                'group_stats',tssg, ...
                'stats_char',struct('str',{str}, ...
                    'offsets',[-0.1,-0.05], ...
                    'font_size',txt_sz, ...
                    'do_display',true));
        end

        %## PLOT
        ax = axes();
        %-- set parameters
        PLOT_STRUCT.group_labels = g_chars_subp;
        PLOT_STRUCT.color_map = color_dark;
        PLOT_STRUCT.cond_labels = xtick_label_g;
        PLOT_STRUCT.title = EEG_MEASURE_TITLES(e_i);
        if e_i == 1
            PLOT_STRUCT.y_label ='10*log_{10}(PSD) - AP. Fit';
        else
            PLOT_STRUCT.y_label ='';
        end
        PLOT_STRUCT.ylim = prc_ylim;
        PLOT_STRUCT.ax_position = [x_shift,y_shift, ...
            AX_W*IM_RESIZE,AX_H*IM_RESIZE];
        %-- group violin plot
        ax = group_violin(tmp_kint,EEG_MEASURES{e_i},'speed_n','group_char',...
            ax,...
            'VIOLIN_STRUCT',VIOLIN_STRUCT,...
            'PLOT_STRUCT',PLOT_STRUCT,...
            'STATS_STRUCT',STATS_STRUCT,...
            'BRACKET_STRUCT',BRACKET_STRUCT,...
            'SIGLINE_STRUCT',SIGLINE_STRUCT, ...
            'CONFINT_STRUCT',CONFINT_STRUCT);
        % ax = group_violin(tmp_fooof_t,EEG_MEASURES{e_i},'cond_id','group_id',...
        %     ax,...
        %     'VIOLIN_PARAMS',VIOLIN_STRUCT,...
        %     'PLOT_STRUCT',VIO_PLOT_STRUCT,...
        %     'STATS_STRUCT',tmp_stats,...
        %     'BRACKET_STRUCT',BRACKET_STRUCT,...
        %     'SIGLINE_STRUCT',SIGLINE_STRUCT);
        if e_i ~= 1
            ylabel('');
        end
        % x_shift = x_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
        % cnt = cnt + 1;    

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
    hold off;
    %##
    % exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_std_allsubj.tiff',string(clusters(c_i)))],...
    %     'Resolution',300)
    exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_%s_muscript_plot.tiff',string(clusters(c_i)),meas_ext)],...
        'Resolution',SAVE_RES)
    % close(fig)
end