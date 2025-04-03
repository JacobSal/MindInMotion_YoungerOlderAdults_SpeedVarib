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
fextr = 'perstridefb_mi_nfslidingb36';
dat = par_load(fpaths{1},sprintf('psd_output_%s.mat',fextr));
%--
fooof_freqs = dat.freqs;
basel_chars = {'perstridefb_mi_nfslidingb36'};
% basel_chars = {'perstridefb_nfslidingb3'};

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
KIN_TABLE = readtable([r_stats_dir filesep sprintf('03122025_lme_eeg_kin_meansd_%s_tbl.xlsx',fextr)], ...
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
violin_ylimits = {{[],[],[],[],[],[],[],[],[],[],[],[]};...
            {[],[],[],[],[],[],[],[],[],[],[],[]};...
            {[],[],[],[],[],[],[],[],[],[],[],[]}};
%--
speed_xvals = (0:5)*0.25;
c_chars = {'0.25 m/s','0.50 m/s','0.75 m/s','1.0 m/s'};
g_chars_topo = {'Young Adults','Older High Func. Adults','Older Low Func. Adults'};
% g_chars_topo = {'Young Adults',{'Older High','Functioning Adults'},{'Older Low','Functioning Adults'}};
g_chars_subp = {'YA','OHFA','OLFA'};
% dip_dir = [cluster_k_dir filesep 'topo_dip_inf' filesep 'all'];
dip_dir = [cluster_k_dir filesep 'topo_dip_inf' filesep 'valid_clusts'];
cmaps_speed = linspecer(4*3);
cmaps_speed = [cmaps_speed(1,:);cmaps_speed(2,:);cmaps_speed(3,:);cmaps_speed(4,:)];
%## EXTRACT PSD DATA
color_dark = cmaps_speed; %color.speed;
color_light = cmaps_speed+0.15; %color.speed_shade;
xtick_label_g = {'0.25','0.50','0.75','1.0'};
x_label = 'speed (m/s)';
cond_offsets = [-0.35,-0.1,0.15,0.40];
%% ===================================================================== %%
%## PARAMETERS
%-
designs = unique(KIN_TABLE.model_n);
group_chars = unique(KIN_TABLE.group_char);
cond_chars = unique(KIN_TABLE.cond_char);
clusters = unique(RSTATS_IMPORT.cluster_num);

%--
SAVE_RES = 300;

%--

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

%## PSD PLOTS
PLOT_STRUCT = struct('y_label',{'10*log_{10}(PSD)'},...
    'y_label_fontsize',10,...
    'y_label_fontweight','bold',...
    'ylim',[],...
    'x_label',{'Frequency (Hz)'},...
    'x_label_fontsize',10,...
    'x_label_fontweight','bold',...
    'x_label_yoffset',0,...
    'xtick_labs',{{}}, ...
    'xticks',[], ...
    'xlim',[],...
    'title',{{''}},...
    'title_fontsize',12,...
    'title_fontweight','normal',...
    'font_size',10,...
    'font_name','Arial',...
    'ax_position',[0,0,1,1],...
    'ax_line_width',1,...
    'xtick_angle',45);
LINE_STRUCT = struct('do_line_avg',false, ...
    'line_width',2, ...
    'line_style','-', ...
    'line_alpha',1, ...
    'line_color',[1,1,1], ...
    'line_label',{'label'}, ...
    'line_avg_fcn',@mean, ...
    'do_err_shading',true, ...
    'err_alpha',0.6, ...
    'err_color',[0.5,0.5,0.5], ...
    'err_edge_color',[], ...
    'err_bnd_vec',[], ...
    'err_line_style',':', ...
    'err_line_width',3);
%% (ALL SUBJS MODEL) =================================================== %%
tmp_savedir = [save_dir filesep fextr '_2by_manufigs'];
mkdir(tmp_savedir);
%--
FIG_INIT_Y = 0;
FIG_INIT_X = 0;
FIG_Y_SHIFT = -0.4;
FIG_X_SHIFT = 0;
FIG_Y_DIM = 2;
%--
fy_cnt = FIG_Y_DIM;
fy_shift = FIG_INIT_Y;
fx_shift = FIG_INIT_X;
%##
for cl_i = 1:length(cluster_inds_plot)
    %% INITIATE TEMPS
    %-- initiate params
    cl_ii = find(cluster_inds_plot(cl_i) == double(string(clusters)));
    cl_n = double(string(clusters(cl_ii)));
    atlas_name = cluster_titles{cl_ii};

    %## INITIATE FIGURE
    if fy_cnt < FIG_Y_DIM
        % fx_shift = fx_shift + FIG_X_SHIFT;
        fy_shift = fy_shift + FIG_Y_SHIFT;
    else
        %-- figure shift
        fy_shift = FIG_Y_INIT;
        fy_cnt = 0;

        %## FIGURE INIT
        fig = figure('color','white', ...
            'Renderer','Painters');
        %-- title
        TITLE_FONT_SIZE = 14;
        FIGURE_POSITION =[1,1,6.5,9];
        FONT_NAME = 'Arial';
        TITLE_XSHIFT = 0.4;
        TITLE_YSHIFT = 0.975+fy_shift;
        TITLE_BOX_SZ = [0.4,0.4];
        annotation('textbox',[TITLE_XSHIFT-(TITLE_BOX_SZ(1)/2/2),TITLE_YSHIFT-(TITLE_BOX_SZ(2)/2), ...
                TITLE_BOX_SZ(1),TITLE_BOX_SZ(2)],...
            'String',atlas_name, ...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle', ...
            'LineStyle','none', ...
            'FontName',FONT_NAME,...
            'FontSize',TITLE_FONT_SIZE, ...
            'FontWeight','Bold', ...
            'Units','normalized');
        set(fig,'Units','inches', ...
            'Position',FIGURE_POSITION, ...
            'PaperUnits','inches', ...
            'PaperSize',[1 1], ...
            'PaperPosition',[0 0 1 1])
        p_sz = get(fig,'Position');
        set(gca,AXES_DEFAULT_PROPS{:});
    end
    fy_cnt = fy_cnt + 1;
    hold on;
    
    %## TOPOGRAPHY PLOT
    IM_RESIZE = 0.2;
    AX_INIT_Y = 0.765+fy_shift;
    AX_INIT_X = 0.085+fx_shift;
    ax_position = [AX_INIT_X,AX_INIT_Y,0,0];
    local_plot_topography(fig,STUDY,cl_n, ...
        group_chars,g_chars,g_chars_topo, ...
        ax_position,IM_RESIZE);
    
    %## DIPOLE PLOT
    LAB_A_YOFFSET = -0.16;
    LAB_A_XOFFSET = 0;
    DIP_IM_DPI = 1000;
    AX_INIT_Y = 0.6+fy_shift; %0.79; %0.8 --- %7.2; % inches for some reason, maybe a bug with normal units
    AX_INIT_X = 0.05+fx_shift; %2.5; % inches for some reason, maybe a bug with normal units
    IM_RESIZE = 1.8;
    dip_fig_path = [dip_dir filesep sprintf('%i_dipplot_alldipspc_angle.fig',cl_n)];
    ax_position = [AX_INIT_X,AX_INIT_Y,0,0];
    
    local_plot_dipole_angle(fig,dip_fig_path,IM_RESIZE, ...
        p_sz,ax_position);

    %## TITLE
    LAB_POS = [AX_INIT_X+LAB_A_XOFFSET+(0.1/2),1+LAB_A_YOFFSET+(0.1/2),.1,.1];
    annotation(fig,'textbox',LAB_POS,...
        'String','(A)', ...
        'HorizontalAlignment','left',...
        'VerticalAlignment','top', ...
        'LineStyle','none', ...
        'FontName','Arial',...
        'FontSize',14, ...
        'FontWeight','Bold', ...
        'Units','normalized');
    %% EXTRACT PSD DATA =============================================== %%    
    %## (STD-MEAN) GET PSD DATA
    dat_out_struct = dat_out_structs{1};
    dat_calcs = {'std','mean'};
    conds_out = {'0p25','0p5','0p75','1p0'};
    groups_out = g_chars;
    [psd_dat_out] = extract_psd_sbs(dat_out_struct,dat_calcs,cl_n,conds_out,groups_out);
    
    %## ALL SUBJECTS & CONDITIONS PER GROUP COMPARISON ====================
    IM_RESIZE = 0.75;
    AX_W = 0.35;
    AX_H = 0.225;
    AX_X_SHIFT = 1.3;
    AX_Y_SHIFT = -1.2;
    AX_INIT_X = 0.375+fx_shift; 
    AX_INIT_Y = 0.75+fy_shift; 
    LEG_X_SHIFT = 0.2; 
    LEG_Y_SHIFT =  0; 
    Y_LIM_SCALE = 1.5;
    X_DIM = 1;
    %--
    % psd_dat_in = cell(size(psd_dat_out,2),1);
    % for g_i = 1:size(psd_dat_out,2)
    %     % psd_dat_in{g_i,1} = cat(2,psd_dat_out{:,g_i});
    %     %--
    %     tmp = cat(3,psd_dat_out{:,g_i});
    %     psd_dat_in{g_i,1} = mean(tmp,3);
    % end
    % params.stats = stats;
    % params.stats.condstats = 'on';
    % params.stats.groupstats = 'off';
    % params.stats.paired = {'off','off'};
    %--
    psd_dat_in = psd_dat_out; 
    params.stats = stats;
    params.stats.condstats = 'on';
    params.stats.groupstats = 'on';
    params.stats.paired = {'on','off'};
    
    %## PLOT
    x_shift = AX_INIT_X;
    y_shift = AX_INIT_Y;
    tmp_plot_struct = PLOT_STRUCT;
    tmp_line_struct = LINE_STRUCT;
    %--
    tmp_plot_struct.ax_position = [x_shift,y_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE];
    % tmp_plot_struct.title = {'Group Comparison'};
    if fy_cnt == 1
        tmp_plot_struct.title = {'Group Comparison'};        
    else
        tmp_plot_struct.title = {''};        
    end
    %{'Combined','Comparisons'};
    tmp_plot_struct.xlim = [3,40];    
    mu = mean(cat(2,psd_dat_in{:}),[2,1]);
    sd = std(cat(2,psd_dat_in{:}),[],[2,1]);
    tmp_plot_struct.ylim = [mu-Y_LIM_SCALE*sd,mu+Y_LIM_SCALE*sd];
    %--
    tmp_line_struct.line_avg_fcn = @(x) mean(x,2);
    tmp_line_struct.line_alpha = 0.7;
    tmp_line_struct.do_err_shading = true;
    tmp_line_struct.err_alpha = 0.2;
    %--
    params.cmaps = linspecer(length(g_chars_subp));
    params.cmaps_scond = cmaps_speed;
    params.cmaps_sgroup = linspecer(length(g_chars_subp));
    params.xtick_label = g_chars_subp;        
    params.freqs = fooof_freqs;
    params.leg_txt_size = 9;
    params.leg_position = [AX_INIT_X+LEG_X_SHIFT*IM_RESIZE*AX_W,...
        y_shift+AX_H*IM_RESIZE+LEG_Y_SHIFT*IM_RESIZE*AX_H];
    params.leg_token_size = 15;
    %--
    ax = axes();
    [paramso] = local_psd_plot_grp(ax,psd_dat_in,params,tmp_plot_struct,tmp_line_struct);
    % [paramso] = local_psd_plot_inter(ax,psd_dat_in,params,tmp_plot_struct,tmp_line_struct);
    xlabel(ax,'');
    ylabel(ax,'<SD(10*log_{10}(PSD_{N}))>');
    
    %% EXTRACT PSD DATA =============================================== %%    
    %## (STD-MEAN) GET PSD DATA
    dat_out_struct = dat_out_structs{1};
    dat_calcs = {'std','mean'};
    conds_out = {'0p25','0p5','0p75','1p0'};
    groups_out = g_chars;
    [psd_dat_out] = extract_psd_sbs(dat_out_struct,dat_calcs,cl_n,conds_out,groups_out);
    
    %## ALL SUBJECTS & CONDITIONS PER GROUP COMPARISON =================
    %--
    psd_dat_in = cell(size(psd_dat_out,1),1);
    for c_i = 1:size(psd_dat_out,1)
        psd_dat_in{c_i,1} = cat(2,psd_dat_out{c_i,:});
    end
    params.stats = stats;
    params.stats.condstats = 'on';
    params.stats.groupstats = 'off';
    params.stats.paired = {'on','off'};
    %--
    % psd_dat_in = psd_dat_out;
    % params.stats = stats;
    % params.stats.condstats = 'on';
    % params.stats.groupstats = 'on';
    % params.stats.paired = {'on','off'};
    
    %## PLOT
    x_shift = AX_INIT_X + AX_W*IM_RESIZE*AX_X_SHIFT;
    y_shift = AX_INIT_Y;
    tmp_plot_struct = PLOT_STRUCT;
    tmp_line_struct = LINE_STRUCT;
    %--
    tmp_plot_struct.ax_position = [x_shift,y_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE];
    if fy_cnt == 1
        tmp_plot_struct.title = {'Condition Comparison'};        
    else
        tmp_plot_struct.title = {''};
    end
    tmp_plot_struct.xlim = [3,40];
    mu = mean(cat(2,psd_dat_in{:}),[2,1]);
    sd = std(cat(2,psd_dat_in{:}),[],[2,1]);
    tmp_plot_struct.ylim = [mu-Y_LIM_SCALE*sd,mu+Y_LIM_SCALE*sd];
    %--
    tmp_line_struct.line_avg_fcn = @(x) mean(x,2);
    tmp_line_struct.line_alpha = 0.7;
    tmp_line_struct.do_err_shading = true;
    tmp_line_struct.err_alpha = 0.2;
    %--
    params.cmaps = cmaps_speed; %linspecer(length(g_chars_subp));
    params.cmaps_scond = cmaps_speed;
    params.cmaps_sgroup = linspecer(length(g_chars_subp));
    params.xtick_label = c_chars; %g_chars_subp;        
    params.freqs = fooof_freqs;
    params.leg_txt_size = 9;
    params.leg_position = [AX_INIT_X+LEG_X_SHIFT*IM_RESIZE*AX_W,...
        y_shift+AX_H*IM_RESIZE+LEG_Y_SHIFT*IM_RESIZE*AX_H];
    params.leg_token_size = 15;
    %--
    ax = axes();
    [paramso] = local_psd_plot_grp(ax,psd_dat_in,params,tmp_plot_struct,tmp_line_struct);
    % [paramso] = local_psd_plot_inter(ax,psd_dat_in,params,tmp_plot_struct,tmp_line_struct);
    xlabel(ax,'');
    ylabel(ax,'');
    %% EXTRACT PSD DATA =============================================== %%    
    %## (MEAN-MEAN) GET PSD DATA
    dat_out_struct = dat_out_structs{1};
    dat_calcs = {'mean','mean'};
    conds_out = {'0p25','0p5','0p75','1p0'};
    groups_out = g_chars;
    [psd_dat_out] = extract_psd_sbs(dat_out_struct,dat_calcs,cl_n,conds_out,groups_out);
    
    %## ALL SUBJECTS & CONDITIONS PER GROUP COMPARISON =================
    %--
    % psd_dat_in = cell(size(psd_dat_out,2),1);
    % for g_i = 1:size(psd_dat_out,2)
    %     % psd_dat_in{g_i,1} = cat(2,psd_dat_out{:,g_i});
    %     %--
    %     tmp = cat(3,psd_dat_out{:,g_i});
    %     psd_dat_in{g_i,1} = mean(tmp,3);
    % end
    % params.stats = stats;
    % params.stats.condstats = 'on';
    % params.stats.groupstats = 'off';
    % params.stats.paired = {'off','off'};
    %--
    psd_dat_in = psd_dat_out; 
    params.stats = stats;
    params.stats.condstats = 'on';
    params.stats.groupstats = 'on';
    params.stats.paired = {'on','off'};
    
    %## PLOT
    x_shift = AX_INIT_X; % + AX_W*IM_RESIZE*AX_X_SHIFT;
    y_shift = AX_INIT_Y + AX_H*IM_RESIZE*AX_Y_SHIFT;
    tmp_plot_struct = PLOT_STRUCT;
    tmp_line_struct = LINE_STRUCT;
    %--
    tmp_plot_struct.ax_position = [x_shift,y_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE];
    tmp_plot_struct.title = {''}; 
    tmp_plot_struct.xlim = [3,40];    
    mu = mean(cat(2,psd_dat_in{:}),[2,1]);
    sd = std(cat(2,psd_dat_in{:}),[],[2,1]);
    tmp_plot_struct.ylim = [mu-Y_LIM_SCALE*sd,mu+Y_LIM_SCALE*sd];
    %--
    tmp_line_struct.line_avg_fcn = @(x) mean(x,2);
    tmp_line_struct.line_alpha = 0.7;
    tmp_line_struct.do_err_shading = true;
    tmp_line_struct.err_alpha = 0.2;
    %--
    params.cmaps = linspecer(length(g_chars_subp));
    params.cmaps_scond = cmaps_speed;
    params.cmaps_sgroup = linspecer(length(g_chars_subp));
    params.xtick_label = g_chars_subp;        
    params.freqs = fooof_freqs;
    params.leg_txt_size = 9;
    params.leg_position = [AX_INIT_X+LEG_X_SHIFT*IM_RESIZE*AX_W,...
        y_shift+AX_H*IM_RESIZE+LEG_Y_SHIFT*IM_RESIZE*AX_H];
    params.leg_token_size = 15;
    %--
    ax = axes();
    % [paramso] = local_psd_plot_grp(ax,psd_dat_in,params,tmp_plot_struct,tmp_line_struct);
    [paramso] = local_psd_plot_inter(ax,psd_dat_in,params,tmp_plot_struct,tmp_line_struct);
    % xlabel(ax,'');
    xlabel(ax,'Frequency (Hz)');
    ylabel(ax,'<<10*log_{10}(PSD_{N})>>');
    
    %% EXTRACT PSD DATA =============================================== %%    
    %## (MEAN-MEAN) GET PSD DATA
    dat_out_struct = dat_out_structs{1};
    dat_calcs = {'mean','mean'};
    conds_out = {'0p25','0p5','0p75','1p0'};
    groups_out = g_chars;
    [psd_dat_out] = extract_psd_sbs(dat_out_struct,dat_calcs,cl_n,conds_out,groups_out);
    
    %## ALL SUBJECTS & CONDITIONS PER GROUP COMPARISON =================
    %--
    psd_dat_in = cell(size(psd_dat_out,1),1);
    for c_i = 1:size(psd_dat_out,1)
        psd_dat_in{c_i,1} = cat(2,psd_dat_out{c_i,:});
    end
    params.stats = stats;
    params.stats.condstats = 'on';
    params.stats.groupstats = 'off';
    params.stats.paired = {'on','off'};
    %--
    % psd_dat_in = psd_dat_out; 
    % params.stats = stats;
    % params.stats.condstats = 'on';
    % params.stats.groupstats = 'on';
    % params.stats.paired = {'on','off'};
    
    %## PLOT
    x_shift = AX_INIT_X + AX_W*IM_RESIZE*AX_X_SHIFT;
    y_shift = AX_INIT_Y + AX_H*IM_RESIZE*AX_Y_SHIFT;
    tmp_plot_struct = PLOT_STRUCT;
    tmp_line_struct = LINE_STRUCT;
    %--
    tmp_plot_struct.ax_position = [x_shift,y_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE];
    tmp_plot_struct.title = {''};
    tmp_plot_struct.xlim = [3,40];
    mu = mean(cat(2,psd_dat_in{:}),[2,1]);
    sd = std(cat(2,psd_dat_in{:}),[],[2,1]);
    tmp_plot_struct.ylim = [mu-Y_LIM_SCALE*sd,mu+Y_LIM_SCALE*sd];
    %--
    tmp_line_struct.line_avg_fcn = @(x) mean(x,2);
    tmp_line_struct.line_alpha = 0.7;
    tmp_line_struct.do_err_shading = true;
    tmp_line_struct.err_alpha = 0.2;
    %--
    params.cmaps = cmaps_speed; 
    params.cmaps_scond = cmaps_speed;
    params.cmaps_sgroup = linspecer(length(g_chars_subp));
    params.xtick_label = c_chars;       
    params.freqs = fooof_freqs;
    params.leg_txt_size = 9;
    params.leg_position = [AX_INIT_X+LEG_X_SHIFT*IM_RESIZE*AX_W,...
        y_shift+AX_H*IM_RESIZE+LEG_Y_SHIFT*IM_RESIZE*AX_H];
    params.leg_token_size = 15;
    %--
    ax = axes();
    [paramso] = local_psd_plot_grp(ax,psd_dat_in,params,tmp_plot_struct,tmp_line_struct);
    % [paramso] = local_psd_plot_inter(ax,psd_dat_in,params,tmp_plot_struct,tmp_line_struct);
    % xlabel(ax,'');
    ylabel(ax,'');
    hold off;
    %%
    fname = sprintf('cl%s_%s_manscript_plot',string(cl_n),meas_ext);
    exportgraphics(fig,[tmp_savedir filesep fname '.pdf'],...
        'ContentType','Vector')
    exportgraphics(fig,[tmp_savedir filesep fname '.tif'],...
        'Resolution',SAVE_RES)
    % close(fig)
end