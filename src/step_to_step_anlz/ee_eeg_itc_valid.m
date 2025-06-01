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
%-
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
            'variable2','cond','values2',{'flat','low','med','high'},'pairing2','on','vartype2','categorical',...
            'variable1','group','values1',{'H1000','H2000','H3000'},'pairing1','off','vartype1','categorical'},...
            {'subjselect',{},...
            'variable2','cond','values2',{'0p25','0p5','0p75','1p0'},'pairing2','on','vartype2','continuous',...
            'variable1','group','values1',{'H1000','H2000','H3000'},'pairing1','off','vartype1','categorical'}};

%## ersp plot per cluster per condition
args = eeglab_struct2args(ERSP_STAT_PARAMS_COND);
CL_STUDY = pop_statparams(CL_STUDY,args{:});
args = eeglab_struct2args(ERSP_PARAMS);
CL_STUDY = pop_erspparams(CL_STUDY,args{:});
CL_STUDY.cache = [];
for des_i = 1:length(STUDY_DESI_PARAMS)
    [CL_STUDY] = std_makedesign(CL_STUDY,[],des_i,STUDY_DESI_PARAMS{des_i}{:});
end

%% PARAMS =========================================================== %%
%-- load dat table
% fext = 'phasec_notw_mw';
fext = 'phasec_notw_mw_based';
% itc_so = par_load([save_dir filesep sprintf('itc_table_%s_new.mat',fext)]);
itc_so = par_load([save_dir filesep sprintf('itc_table_%s.mat',fext)]);
%--
twp = struct('timewarpms',itc_so.twg_times(1,:));
% twp = itc_so.twg_times(1,:);
%--
TW_CHARS = {'RHS','RTO','LHS','LTO','RHS'};
% FREQ_BOUND = [3,60];
% TIME_BOUND = [twp.timewarpms(1),twp.timewarpms(end)];
FREQ_BOUND = [3,250];
TIME_BOUND = [-500,3000];
alltimes = itc_so.itc_times{1};
allfreqs = itc_so.itc_freqs{1};
tinds = alltimes > TIME_BOUND(1) & alltimes < TIME_BOUND(2);
finds = allfreqs > FREQ_BOUND(1) & allfreqs < FREQ_BOUND(2);
tcrop = alltimes(tinds);
fcrop = allfreqs(finds);
%%
%##
fext = 'phasec_notw_mw_based';
clusters = unique(itc_so.cluster_n);
subjects = unique(itc_so.subj_char);
groups = unique(itc_so.group_char);
%## PARAMETERS
cluster_inds_plot = [3,4,6,8];
COND_CHARS = {'0p25','0p5','0p75','1p0'};
alltitles = {'0.25 m/s','0.50 m/s','0.75 m/s','1.0 m/s'};
cluster_titles = {'Right Posterior Parietal', ...
    'Right Sensorimotor', ...
    'Left Precuneus', ... %'Anterior Cingulate', ...
    'Left Sensorimotor', ...
    'Right Premotor',...
    'Left Posterior Parietal', ...
    'Left Supplementary Motor', ...
    'Right Occipital', ...
    'Mid Cingulate',...
    'Left Temporal',...
    'Left Occipital'};

%##
pgsz = [1,1,6.5,9];
psr = 9/6.5;
label_struct = struct('FontName','Arial', ...
    'FontSize',8, ...
    'FontWeight','bold');
ax_struct = struct('FontName','Arial', ...
    'FontSize',8, ...
    'FontWeight','normal', ...
    'OuterPosition',[0,0,1,1], ...
    'LineWidth',1, ...
    'Position',[0.06,0.7,0.13,0.16]);
PLOT_STRUCT = struct( ...
    'title',{'ERSP'},...
    'title_props',label_struct, ...
    'xlabel','Gait Cycle Time (ms)',...
    'xlabel_props',label_struct, ...
    'ylabel','Frequency (Hz)',...
    'ylabel_props',label_struct, ...
    'xticklabel_times',[twp.timewarpms],...
    'xticklabel_chars',{TW_CHARS},...
    'xticklabel_angle',45,...
    'ax_props',ax_struct, ...
    'clim',[],...
    'freq_lims',[],...
    'time_lims',[],...    
    'contourf_grain',ceil((500/pi())),...
    'alpha_multiple',0.6,...
    'do_display_sgtitle',true, ...
    'sgtitle_char',{''}, ...
    'sgtitle_shift',[0,0.65], ...
    'sgtitle_boxsz',[0.1,0.1], ...
    'sgtitle_props',struct(...
        'LineStyle','none',...
        'FontName','Arial', ...
        'FontSize',12,...
        'FontWeight','bold',...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top',...
        'Units','normalized'), ...
    'do_display_colorbar',true, ...
    'cbar_shift',[0.8*psr,psr*(1/psr)],...
    'cbar_label_shift',[1.65*psr,1.27*(1/psr)],...
    'cbar_ticks',[],...
    'cbar_ticklabs',{{''}},...
    'cbar_label','\Delta Power (dB)',...
    'cbar_label_props',struct( ...
        'LineStyle','none',...
        'FontName','Arial', ...
        'FontSize',8,...
        'FontWeight','bold',...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle',...
        'Units','normalized',...
        'Rotation',270), ...
    'do_display_bandmarks',true);
AXES_DEFAULT_PROPS = {'box','off', ...
    'xtick',[], ...
    'ytick',[],...
    'ztick',[], ...
    'xcolor',[1,1,1], ...
    'ycolor',[1,1,1]};
BOOT_STRUCT = struct(...
    'niters',1000, ...
    'alpha',0.05, ...
    'cluster_thresh',300);
%--
ieee_sz = [8.5-(0.65*2),11-(0.7*2)];
FIGURE_POSITION =[1,1,ieee_sz];
AX_DIM = [0.13,0.16];
AX_SHIFT = [1.2,-1.25];
AX_INIT_X = 0.2;
AX_INIT_Y = 0.775;
X_DIM = 4;
fy_shift = 0;

%##
tmp_save_dir = [save_dir filesep 'itc_raw_plots'];
mkdir(tmp_save_dir);

%% (MORE EFFICIENT PLOTTING) ============================================
%## LOOP
for cl_i = 1:length(cluster_inds_plot)
    %##
    blext = 'none';
    stat_ext = 'none';    
    %--
    cl_ii = find(cluster_inds_plot(cl_i) == double(string(clusters)));
    cl_n = clusters(cl_ii);
    atlas_name = cluster_titles{cl_ii};

    %## EXTRACT DATA
    allitc = cell(length(COND_CHARS),length(groups));
    cropitc = cell(length(COND_CHARS),length(groups));
    %--
    for g_i = 1:length(groups)
        for c_i = 1:length(COND_CHARS)
            inds = itc_so.cluster_n == cl_n & ...
                strcmp(itc_so.group_char,groups{g_i}) & ...
                strcmp(itc_so.cond_char,COND_CHARS{c_i});
            tt = itc_so(inds,:);
            %--
            % itc_dat = cellfun(@abs,tt.itc_dat,'UniformOutput',false);
            % itc_dat = cat(3,itc_dat{:});
            %--
            % itc_dat = cellfun(@abs,tt.itc_dat_slide,'UniformOutput',false);
            % itc_dat = cat(3,itc_dat{:});

            %--
            itc_dat = cellfun(@abs,tt.ersp_dat_slide,'UniformOutput',false);
            itc_dat = cat(3,itc_dat{:});
            %--
            % itc_dat = cellfun(@abs,tt.erspb_dat_slide,'UniformOutput',false);
            % itc_dat = cat(3,itc_dat{:});

            %--
            allitc{c_i,g_i} = itc_dat;
            cropitc{c_i,g_i} = itc_dat(finds,tinds,:);
        end
    end
    %--
    itco_c = cropitc;
    %--
    % [itco_f,itco_c] = eeglab_baseln(allitc,alltimes,allfreqs,TIME_BOUND,FREQ_BOUND, ...
        % 'DO_COMMON_BASE',true, ...
        % 'DO_SUBJ_BASE',true);
    % blext = 'com';
    % blext = 'sub';
    % blext = 'cs';

    %-- clims
    % clim = cellfun(@(x) [prctile(x,3,'all'),prctile(x,97,'all')],itco_c, ...
    %     'UniformOutput',false);
    % clim = mean(cat(1,clim{:}),1);
    clim = [0,0.25];

   %## INITIATE FIGURE
    x_shift = AX_INIT_X;
    y_shift = AX_INIT_Y;
    x_cnt = 1;
    y_cnt = 1;
    %-- fig
    TITLE_XSHIFT = 0.4;
    TITLE_YSHIFT = 0.975+fy_shift;
    TITLE_BOX_SZ = [0.4,0.4];
    fig = figure('color','white');
    set(fig,'Units','inches', ...
        'Position',FIGURE_POSITION, ...
        'PaperUnits','inches', ...
        'PaperSize',[1 1], ...
        'PaperPosition',[0 0 1 1])
    %-- title
    annotation('textbox',[TITLE_XSHIFT-(TITLE_BOX_SZ(1)/2/2),TITLE_YSHIFT-(TITLE_BOX_SZ(2)/2), ...
            TITLE_BOX_SZ(1),TITLE_BOX_SZ(2)],...
        'String',atlas_name, ...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle', ...
        'LineStyle','none', ...
        'FontName','Arial',...
        'FontSize',14, ...
        'FontWeight','Bold', ...
        'Units','normalized');
    p_sz = get(fig,'Position');
    set(gca,AXES_DEFAULT_PROPS{:});
    %-- axes
    hold on;
    %--
    % rsubj = randi(size(itco_c{1,1},3),1);
    for g_i = 1:length(groups)
        for c_i = 1:length(COND_CHARS)
            %--
            % allersp = bs_ersp{c_i,g_i};
            % allersp_mask = bs_masked{c_i,g_i};
            % allersp_pcond = bs_pval{c_i,g_i};
            %--
            allersp = mean(itco_c{c_i,g_i},3);
            allersp_mask = mean(itco_c{c_i,g_i},3);
            allersp_pcond = zeros(size(itco_c{c_i,g_i},1),size(itco_c{c_i,g_i},2));

            %## PLOT
            tmp_plot_struct = PLOT_STRUCT;
            tmp_plot_struct.colormap = linspecer(); %cmp(ceil(length(cmp)/2):length(cmp),:);
            tmp_plot_struct.alpha_multiple=0.9;
            tmp_plot_struct.title = alltitles{c_i};
            tmp_plot_struct.ax_props.FontSize = 8;
            tmp_plot_struct.contourf_grain = 30;
            tmp_plot_struct.ax_props.Position = [x_shift,y_shift,AX_DIM(1),AX_DIM(2)];
            tmp_plot_struct.clim = clim;
            % tmp_plot_struct.clim = [0,0.12]; %clim
            % tmp_plot_struct.cbar_ticklabs = {'0','0.04','0.08','0.12'};
            % tmp_plot_struct.cbar_ticks = [0,0.04,0.08,0.12];
            if c_i > 1
                tmp_plot_struct.do_display_bandmarks = false;
                tmp_plot_struct.ylabel = '';
            end
            if c_i < length(COND_CHARS)
                tmp_plot_struct.do_display_colorbar = false;
            end
            if g_i < length(groups)
                tmp_plot_struct.xticklabel_chars = {''};
                tmp_plot_struct.xlabel = '';
            elseif c_i > 1 && g_i == length(groups)
                tmp_plot_struct.xlabel = '';
            end
            ax = axes();
            [ax] = plot_contourf_cell(ax,allersp,tcrop,fcrop,...
                allersp_mask,allersp_pcond, ...
                'PLOT_STRUCT',tmp_plot_struct);
            %## AX SHIFT
            if x_cnt < X_DIM
                x_shift = x_shift + AX_SHIFT(1)*AX_DIM(1);
            else
                y_shift = y_shift + AX_SHIFT(2)*AX_DIM(2);
                x_shift = AX_INIT_X;
                x_cnt = 0;
                y_cnt = y_cnt + 1;
            end
            x_cnt = x_cnt + 1;
        end
    end
    drawnow;
    exportgraphics(fig,[tmp_save_dir filesep sprintf('cl%i_%s_%s_%s_clim.png',cl_n,fext,stat_ext,blext)],'Resolution',300);
    % close(fig);
end

%% MOD INDEX =========================================================== %%
%--
ALPHA = 0.05;
mi_pfreq_vec = itc_dato.itc_dat_struct.mi_phase;
mi_afreq_vec = itc_dato.itc_dat_struct.mi_freqs;

%## PARAMETERS
cluster_inds_plot = [3,4,6,8];
clusters = unique(itc_so.cluster_n);
subjects = unique(itc_so.subj_char);
groups = unique(itc_so.group_char);
COND_CHARS = {'0p25','0p5','0p75','1p0'};
alltitles = {'0.25 m/s','0.50 m/s','0.75 m/s','1.0 m/s'};
cluster_titles = {'Right Posterior Parietal', ...
    'Right Sensorimotor', ...
    'Left Precuneus', ... %'Anterior Cingulate', ...
    'Left Sensorimotor', ...
    'Right Premotor',...
    'Left Posterior Parietal', ...
    'Left Supplementary Motor', ...
    'Right Occipital', ...
    'Mid Cingulate',...
    'Left Temporal',...
    'Left Occipital'};

%##
pgsz = [1,1,6.5,9];
psr = 9/6.5;
label_struct = struct('FontName','Arial', ...
    'FontSize',8, ...
    'FontWeight','bold');
ax_struct = struct('FontName','Arial', ...
    'FontSize',8, ...
    'FontWeight','normal', ...
    'OuterPosition',[0,0,1,1], ...
    'LineWidth',1, ...
    'Position',[0.06,0.7,0.13,0.16]);
PLOT_STRUCT = struct( ...
    'title',{'ERSP'},...
    'title_props',label_struct, ...
    'xlabel','Gait Cycle Time (ms)',...
    'xlabel_props',label_struct, ...
    'ylabel','Frequency (Hz)',...
    'ylabel_props',label_struct, ...
    'xticklabel_times',(min(mi_pfreq_vec):10:max(mi_pfreq_vec)),...
    'xticklabel_chars',{cellstr(string((min(mi_pfreq_vec):5:max(mi_pfreq_vec))))},...
    'xticklabel_angle',45,...
    'ax_props',ax_struct, ...
    'clim',[],...
    'freq_lims',[min(mi_afreq_vec),max(mi_afreq_vec)],...
    'time_lims',[min(mi_pfreq_vec),max(mi_pfreq_vec)],...    
    'contourf_grain',ceil((500/(2*pi()))),...
    'alpha_multiple',0.6,...
    'do_display_sgtitle',true, ...
    'sgtitle_char',{''}, ...
    'sgtitle_shift',[0,0.65], ...
    'sgtitle_boxsz',[0.1,0.1], ...
    'sgtitle_props',struct(...
        'LineStyle','none',...
        'FontName','Arial', ...
        'FontSize',12,...
        'FontWeight','bold',...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top',...
        'Units','normalized'), ...
    'do_display_colorbar',true, ...
    'cbar_shift',[0.8*psr,psr*(1/psr)],...
    'cbar_label_shift',[1.65*psr,1.27*(1/psr)],...
    'cbar_ticks',[],...
    'cbar_ticklabs',{{''}},...
    'cbar_label','\Delta Power (dB)',...
    'cbar_label_props',struct( ...
        'LineStyle','none',...
        'FontName','Arial', ...
        'FontSize',8,...
        'FontWeight','bold',...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle',...
        'Units','normalized',...
        'Rotation',270), ...
    'do_display_bandmarks',true);
AXES_DEFAULT_PROPS = {'box','off', ...
    'xtick',[], ...
    'ytick',[],...
    'ztick',[], ...
    'xcolor',[1,1,1], ...
    'ycolor',[1,1,1]};
% BOOT_STRUCT = struct(...
%     'niters',1000, ...
%     'alpha',0.05, ...
%     'cluster_thresh',300);
%--
GROUPT_CHARS = {'Younger Adults','Older Higher Functioning Adults','Older Lower Functioning Adults'};
GROUPT_SHIFT = [0,0.5];
GROUPT_PROPS = {...
    'LineStyle','none',...
    'FontName','Arial', ...
    'FontSize',10,...
    'FontWeight','bold',...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','top',...
    'Units','normalized'...
    };
%--
ieee_sz = [8.5-(0.65*2),11-(0.7*2)];
FIGURE_POSITION = ieee_sz;
AX_DIM = [0.13,0.16];
AX_SHIFT = [1.2,-1.25];
AX_INIT_X = 0.2;
AX_INIT_Y = 0.775;
X_DIM = 4;
fy_shift = 0;
%##
tmp_save_dir = [save_dir filesep 'mi_plots'];
mkdir(tmp_save_dir);

%% (MORE EFFICIENT PLOTTING) ============================================
%## LOOP
for cl_i = 1:length(cluster_inds_plot)
    %##
    blext = 'none';
    stat_ext = 'none';    
    %--
    cl_ii = find(cluster_inds_plot(cl_i) == double(string(clusters)));
    cl_n = clusters(cl_i);
    atlas_name = cluster_titles{cl_n};

    %## EXTRACT DATA
    allitc = cell(length(COND_CHARS),length(groups));
    cropitc = cell(length(COND_CHARS),length(groups));
    %--
    for g_i = 1:length(groups)
        for c_i = 1:length(COND_CHARS)
            inds = itc_so.cluster_n == cl_n & ...
                strcmp(itc_so.group_char,groups{g_i}) & ...
                strcmp(itc_so.cond_char,COND_CHARS{c_i});
            tt = itc_so(inds,:);
            % itc_dat = cellfun(@real,tt.itc_dat,'UniformOutput',false);
            % itc_dat = cellfun(@abs,tt.mi_dat,'UniformOutput',false);
            itc_dat = tt.mi_dat;
            itc_dat = cat(3,itc_dat{:});
            %--
            allitc{c_i,g_i} = itc_dat;
            % cropitc{c_i,g_i} = itc_dat(finds,tinds,:);
        end
    end
    %--
    itco_c = allitc;
    %--
    % [itco_f,itco_c] = eeglab_baseln(allitc,alltimes,allfreqs,TIME_BOUND,FREQ_BOUND, ...
        % 'DO_COMMON_BASE',true, ...
        % 'DO_SUBJ_BASE',true);
    % blext = 'com';
    % blext = 'sub';
    % blext = 'cs';
    %--
    % itco_c = cellfun(@(x) std(x,[],3)/size(x,3),itco_c,'UniformOutput',false);

    %##
    stats = CL_STUDY.etc.statistics;
    stats.condstats = 'on';
    stats.groupstats = 'off';
    stats.paired = {'on','off'};
    %--
    % stats = CL_STUDY.etc.statistics;
    % stats.condstats = 'off';
    % stats.groupstats = 'on';
    % stats.paired = {'on','off'};
    %--
    [pcond,pgroup,pinter,scond,sgroup,sinter] = std_stat(itco_c, stats);
    stat_ext = [stats.condstats,stats.groupstats];
    %##
    %-- clims
    % clim = cellfun(@(x) [prctile(x,3,'all'),prctile(x,97,'all')],itco_c, ...
    %     'UniformOutput',false);
    % clim = mean(cat(1,clim{:}),1);
    % clim = double([-max(abs(clim)),max(abs(clim))]);
    clim = [0,5e-5];

    %## INITIATE FIGURE
    x_shift = AX_INIT_X;
    y_shift = AX_INIT_Y;
    x_cnt = 1;
    y_cnt = 1;
    %-- fig
    ieee_sz = [8.5-(0.65*2),11-(0.7*2)];
    TITLE_XSHIFT = 0.4;
    TITLE_YSHIFT = 0.975+fy_shift;
    TITLE_BOX_SZ = [0.4,0.4];
    fig = figure('color','white');
    set(fig,'Units','inches', ...
        'Position',ieee_sz, ...
        'PaperUnits','inches', ...
        'PaperSize',[1 1], ...
        'PaperPosition',[0 0 1 1])
    %-- title
    annotation('textbox',[TITLE_XSHIFT-(TITLE_BOX_SZ(1)/2/2),TITLE_YSHIFT-(TITLE_BOX_SZ(2)/2), ...
            TITLE_BOX_SZ(1),TITLE_BOX_SZ(2)],...
        'String',atlas_name, ...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle', ...
        'LineStyle','none', ...
        'FontName','Arial',...
        'FontSize',14, ...
        'FontWeight','Bold', ...
        'Units','normalized');
    p_sz = get(fig,'Position');
    set(gca,AXES_DEFAULT_PROPS{:});
    %-- axes
    hold on;
    %--
    % rsubj = randi(size(itco_c{1,1},3),1);
    for g_i = 1:length(groups)
        for c_i = 1:length(COND_CHARS)            
            % tmpdat = itco_c{c_i,g_i}(:,:,rsubj);
            %--
            allersp = mean(itco_c{c_i,g_i},3);
            allersp_mask = mean(itco_c{c_i,g_i},3);
            allersp_pcond = zeros(size(allersp));
            %--
            % if ~isempty(pcond)
            %     tmpp = double(pcond{g_i});
            %     tmps = double(scond{g_i});
            % elseif ~isempty(pgroup)
            %     tmpp = double(pgroup{c_i});
            %     tmps = double(sgroup{c_i});
            % end
            % allersp = mean(itco_c{c_i,g_i},3);
            % allersp_mask = allersp.*tmpp;
            % allersp_pcond = tmps;
            
            %## PLOT
            tmp_plot_struct = PLOT_STRUCT;
            tmp_plot_struct.colormap = linspecer(); %cmp(ceil(length(cmp)/2):length(cmp),:);
            tmp_plot_struct.alpha_multiple=0.9;
            tmp_plot_struct.title = alltitles{c_i};
            tmp_plot_struct.ax_props.FontSize = 8;
            tmp_plot_struct.contourf_grain = 30;
            tmp_plot_struct.ax_props.Position = [x_shift,y_shift,AX_DIM(1),AX_DIM(2)];
            tmp_plot_struct.clim = clim;
            % tmp_plot_struct.clim = [0,0.12]; %clim
            % tmp_plot_struct.cbar_ticklabs = {'0','0.04','0.08','0.12'};
            % tmp_plot_struct.cbar_ticks = [0,0.04,0.08,0.12];
            if c_i > 1
                tmp_plot_struct.do_display_bandmarks = false;
                tmp_plot_struct.ylabel = '';
            end
            if c_i < length(COND_CHARS)
                tmp_plot_struct.do_display_colorbar = false;
            end
            if g_i < length(groups)
                tmp_plot_struct.xticklabel_chars = {''};
                tmp_plot_struct.xlabel = '';
            elseif c_i > 1 && g_i == length(groups)
                tmp_plot_struct.xlabel = '';
            end
            %--
            if c_i == 1
                GROUPTITLE_BOXSIZE = [0.5,0.1];
                xx = 0.5+(-GROUPTITLE_BOXSIZE/2)+GROUPT_SHIFT(1);
                yy = y_shift+tmp_plot_struct.ax_props.Position(2)*GROUPT_SHIFT(2);
                a = annotation(gcf,'textbox',[xx,yy,GROUPTITLE_BOXSIZE],...
                    'String',GROUPT_CHARS{g_i}, ...
                    GROUPT_PROPS{:});
            end
            %--
            ax = axes();
            [ax] = plot_contourf_cell(ax,allersp',mi_pfreq_vec,mi_afreq_vec,...
                allersp_mask,allersp_pcond, ...
                'PLOT_STRUCT',tmp_plot_struct);
            %## AX SHIFT
            if x_cnt < X_DIM
                x_shift = x_shift + AX_SHIFT(1)*AX_DIM(1);
            else
                y_shift = y_shift + AX_SHIFT(2)*AX_DIM(2);
                x_shift = AdX_INIT_X;
                x_cnt = 0;
                y_cnt = y_cnt + 1;
            end
            x_cnt = x_cnt + 1;
        end
    end
    drawnow;
    exportgraphics(fig,[tmp_save_dir filesep sprintf('cl%i_mi_%s_%s_%s.tiff',cl_n,fext,stat_ext,blext)],'Resolution',300);
    % close(fig);
end

%% IMPORT R FLASSO ===================================================== %%
clusters_ext = [3,4,6,8];
% conds_ext = [5,6,7,8];
conds_ext = [1,2,3,4];
% fext = 'itc_rdata_table_phasec_notw_mw_based';
% fext = 'itc_rdata_table_phasec_notw_mw_based_flasso_results_bsz5';
% fext = 'itc_rdata_table_phasec_notw_mw_based_fl_res_bsz5_nob';
% fext = 'itc_rdata_flasso_out_bsz5_nob_sliding';
fext = 'itc_rdata_flasso_125f_out_bsz5_nob_sliding';
itc_dat_masks = par_load(save_dir,sprintf('itc_rdata_table_%s_cl%s_c%s.mat',fext,strjoin(string(conds_ext),''),strjoin(string(clusters_ext),'')));
%--
COND_CHARS = {'0p25','0p5','0p75','1p0'};
FREQ_BOUND = [3,60];
TIME_BOUND = [twp.timewarpms(1),twp.timewarpms(end)];
% FREQ_BOUND = [3,250];
% TIME_BOUND = [-500,3000];
itc_times = double(itc_dat_masks.itc_times{1});
itc_freqs = double(itc_dat_masks.itc_freqs{1});
tinds = itc_times > TIME_BOUND(1) & itc_times < TIME_BOUND(2);
finds = itc_freqs > FREQ_BOUND(1) & itc_freqs < FREQ_BOUND(2);
tcrop = itc_times(tinds);
fcrop = itc_freqs(finds);
%--
ALPHA = 0.05;

cluster_inds_plot = clusters_ext;
clusters = (3:11); %unique(itc_dat_masks.cluster_n);
subjects = unique(itc_dat_masks.subj_char);
groups = unique(itc_dat_masks.group_char);

%## PARAMETERS
alltitles = {'0.25 m/s','0.50 m/s','0.75 m/s','1.0 m/s'};
c_alltitles = {'','0.50 ms^{-1}-0.25 ms^{-1}','0.75 ms^{-1}-0.25 ms^{-1}','1.0 ms^{-1}-0.25 ms^{-1}'};
cluster_titles = {'Right Posterior Parietal', ...
    'Right Sensorimotor', ...
    'Left Precuneus', ... %'Anterior Cingulate', ...
    'Left Sensorimotor', ...
    'Right Premotor',...
    'Left Posterior Parietal', ...
    'Left Supplementary Motor', ...
    'Right Occipital', ...
    'Mid Cingulate',...
    'Left Temporal',...
    'Left Occipital'};
%##
pgsz = [1,1,6.5,9];
psr = 9/6.5;
label_struct = struct('FontName','Arial', ...
    'FontSize',8, ...
    'FontWeight','bold');
ax_struct = struct('FontName','Arial', ...
    'FontSize',8, ...
    'FontWeight','normal', ...
    'OuterPosition',[0,0,1,1], ...
    'LineWidth',1, ...
    'Position',[0.06,0.7,0.13,0.16]);
PLOT_STRUCT = struct( ...
    'title',{'ERSP'},...
    'title_props',label_struct, ...
    'xlabel','Gait Cycle Time (ms)',...
    'xlabel_props',label_struct, ...
    'ylabel','Frequency (Hz)',...
    'ylabel_props',label_struct, ...
    'xticklabel_times',[twp.timewarpms],...
    'xticklabel_chars',{TW_CHARS},...
    'xticklabel_angle',45,...
    'ax_props',ax_struct, ...
    'clim',[],...
    'freq_lims',[],...
    'time_lims',[],...    
    'contourf_grain',ceil(500/(2*pi())),...
    'alpha_multiple',0.6,...
    'do_display_sgtitle',true, ...
    'sgtitle_char',{''}, ...
    'sgtitle_shift',[0,0.65], ...
    'sgtitle_boxsz',[0.1,0.1], ...
    'sgtitle_props',struct(...
        'LineStyle','none',...
        'FontName','Arial', ...
        'FontSize',12,...
        'FontWeight','bold',...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top',...
        'Units','normalized'), ...
    'do_display_colorbar',true, ...
    'cbar_shift',[0.8*psr,psr*(1/psr)],...
    'cbar_label_shift',[1.65*psr,1.27*(1/psr)],...
    'cbar_ticks',[],...
    'cbar_ticklabs',{{''}},...
    'cbar_label','\Delta Power (dB)',...
    'cbar_label_props',struct( ...
        'LineStyle','none',...
        'FontName','Arial', ...
        'FontSize',8,...
        'FontWeight','bold',...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle',...
        'Units','normalized',...
        'Rotation',270), ...
    'do_display_bandmarks',true);
AXES_DEFAULT_PROPS = {'box','off', ...
    'xtick',[], ...
    'ytick',[],...
    'ztick',[], ...
    'xcolor',[1,1,1], ...
    'ycolor',[1,1,1]};
BOOT_STRUCT = struct(...
    'niters',1000, ...
    'alpha',0.05, ...
    'cluster_thresh',300);

%##
tmp_save_dir = [save_dir filesep 'itc_plots_flasso'];
mkdir(tmp_save_dir);

%% (MORE EFFICIENT PLOTTING) =========================================== %%
%## LOOP
% COND_CHARS = {'flat','low','med','high'};
COND_CHARS = {'0p25','0p5','0p75','1p0'};
MOD_CUTOFF = 0.08;
for cl_i = 1:length(cluster_inds_plot)
    %%
    
    fext = 'mod_group';
    blext = 'none';
    stat_ext = 'none';    
    %--
    cl_ii = find(cluster_inds_plot(cl_i) == double(string(clusters)));
    cl_n = clusters(cl_ii);
    atlas_name = cluster_titles{cl_ii};

    %## LOAD STATS
    fexts = 'itc_rdata_flasso_125f_out_bsz5_nob_sliding';
    %--
    % stat_ext = 'int_rafdr';
    % stat_ext = 'int_rindvfdr';
    % stat_ext = 'int_rrawp';
    % flasso_fpath = [r_stats_dir filesep fexts filesep sprintf('statmatint_cl%i.mat',cl_n)];
    %--
    stat_ext = 'grp_rafdr';
    % stat_ext = 'grp_rindvfdr';
    flasso_fpath = [r_stats_dir filesep fexts filesep sprintf('statmatgrp_cl%i.mat',cl_n)];
    %--
    % tmp = load(flasso_fpath,'-mat');
    tmp = load(flasso_fpath);
    tmp = tmp.stat_mat;
    %--
    estm = reshape(tmp.estimate,[length(tmp.freqs),length(tmp.times),length(tmp.coeff_c)]);
    pval = reshape(tmp.fdrp,[length(tmp.freqs),length(tmp.times),length(tmp.pinds)]);
    % pval = reshape(tmp.rawp,[length(tmp.freqs),length(tmp.times),length(tmp.coeff_c)]);
    %--
    astat = reshape(tmp.astat,[length(tmp.freqs),length(tmp.times),length(tmp.acoeff_c)]);
    apval = reshape(tmp.afdrp(:,1),[length(tmp.freqs),length(tmp.times),length(tmp.apinds)]);
    % apval = reshape(tmp.arawp,[length(tmp.freqs),length(tmp.times),length(tmp.acoeff_c)]);
    %--
    TIME_BOUND = [min(tmp.times),max(tmp.times)];
    FREQ_BOUND = [min(tmp.freqs),max(tmp.freqs)];
    tinds = itc_times > TIME_BOUND(1) & itc_times < TIME_BOUND(2);
    finds = itc_freqs > FREQ_BOUND(1) & itc_freqs < FREQ_BOUND(2);
    tcrop = itc_times(tinds);
    fcrop = itc_freqs(finds);
    %--
    pval(isnan(pval)) = 1;
    estm(isnan(estm)) = 0;
    %--
    % stat_ext = [stat_ext,'_anv'];
    % stattitles = {'Speed','Group','Speed:Group'};
    % allpvs = {apval(:,:,1),apval(:,:,2),apval(:,:,3)};
    % allest = {astat(:,:,1),astat(:,:,2),astat(:,:,3)};
    % clim_asi = [0,0,0];
    % clim_tit = {'F Stat','F Stat','F Stat'};
    %--
    % stat_ext = [stat_ext,'_est'];
    % stattitles = {'\DeltaSpeed_{YA}','OHFA-YA','OLFA-YA','\DeltaSpeed_{OHFA-YA}','\DeltaSpeed_{OLFA-YA}'};
    % allpvs = {pval(:,:,2),pval(:,:,3),pval(:,:,4),pval(:,:,5),pval(:,:,6)};
    % allest = {estm(:,:,2),estm(:,:,3),estm(:,:,4),estm(:,:,5),estm(:,:,6)};
    % clim_asi = [1,2,2,3,3];
    % clim_tit = {'slope (m_{YA})','intercept (b_{OHFA}-b_{YA})','intercept (b_{OLFA}-b_{YA})','slope (m_{OHFA}-m_{YA})','slope (m_{OLFA}-m_{YA})'};
    %--
    stat_ext = [stat_ext,'_anv'];
    stattitles = {'Speed','Group','OHFA-YA','OLFA-YA'};
    allpvs = {apval(:,:,1),apval(:,:,2),pval(:,:,3),pval(:,:,4)};
    allest = {astat(:,:,1),astat(:,:,2),estm(:,:,3),estm(:,:,4)};
    clim_asi = [4,5,2,2];
    clim_tit = {'F Stat','F Stat','intercept (ITC)','intercept (ITC)'};
    %--
    % stat_ext = [stat_ext,'_est'];
    % stattitles = {'\DeltaSpeed_{YA}','OHFA-YA','OLFA-YA'};
    % allpvs = {pval(:,:,2),pval(:,:,3),pval(:,:,4)};
    % allest = {estm(:,:,2),estm(:,:,3),estm(:,:,4)};
    % clim_asi = [1,2,2,3,3];
    % clim_tit = {'slope (m_{YA})','intercept (b_{OHFA}-b_{YA})','intercept (b_{OLFA}-b_{YA})'};
    
    %## EXTRACT DATA
    allitc = cell(length(COND_CHARS),length(groups));
    allmod = cell(length(COND_CHARS),length(groups));
    allmodb = cell(length(COND_CHARS),length(groups));
    %--
    for g_i = 1:length(groups)
        for c_i = 1:length(COND_CHARS)
            inds = itc_dat_masks.cluster_n == cl_n & ...
                strcmp(itc_dat_masks.group_char,groups{g_i}) & ...
                strcmp(itc_dat_masks.cond_char,COND_CHARS{c_i});
            tt = itc_dat_masks(inds,:);
            %--
            itc_mod = tt.itc_mod;
            itc_mod = cellfun(@(x) x(finds,tinds),itc_mod,'UniformOutput',false);
            itc_mod = cat(3,itc_mod{:});
            %--
            % itc_mod = tt.itc_dat;
            % itc_mod = cellfun(@(x) x(finds,tinds),itc_mod,'UniformOutput',false);
            % itc_mod = cat(3,itc_mod{:});
            % fext = 'itc_group';
            %--
            allmod{c_i,g_i} = itc_mod;
        end
    end
    %--
    itco_c = allmod;
    
    %-- common
    % blc = mean(cat(3,itco_c{:}),3);
    % itco_c = cellfun(@(x) x-blc,itco_c,'UniformOutput',false);

    %## CONDITION & GROUP SPECIFIC BASELINES
    % tmp = cell(size(itco_c));
    % for g_i = 1:size(itco_c,2)
    %     datmu = mean(cat(3,itco_c{:,g_i}),3);
    %     tmp(:,g_i) = cellfun(@(x) x-datmu,itco_c(:,g_i),'UniformOutput',false);
    % end
    %--
    % tmp = cell(size(itco_c));
    % for c_i = 1:size(itco_c,1)
    %     datmu = mean(cat(3,itco_c{c_i,:}),3);
    %     tmp(c_i,:) = cellfun(@(x) x-datmu,itco_c(c_i,:),'UniformOutput',false);
    % end
    %--
    % tmp = cell(size(itco_c));
    % blg = cellfun(@(x) mean(x,3),itco_c(:,3),'UniformOutput',false);
    % for c_i = 1:size(itco_c,1)
    %     datmu = blg{c_i};
    %     tmp(c_i,:) = cellfun(@(x) x-datmu,itco_c(c_i,:),'UniformOutput',false);
    % end
    %-- baseline to 0.25 m/s for each group
    % tmp = cell(size(itco_c));
    % blc = cellfun(@(x) mean(x,3),itco_c(1,:),'UniformOutput',false);
    % for g_i = 1:size(itco_c,2)
    %     datmu = blc{g_i};
    %     tmp(:,g_i) = cellfun(@(x) x-datmu,itco_c(:,g_i),'UniformOutput',false);
    % end
    % blext = "condd0p25";
    % ttits = c_alltitles;
    %-- assign
    % itco_c = tmp;
    
    %## STANDARD ERROR
    % itco_c = cellfun(@(x) std(x,[],3)/size(x,3),itco_c,'UniformOutput',false);

    %##
    %-- clims
    % clim = cellfun(@(x) [prctile(x,5,'all'),prctile(x,95,'all')],itco_c, ...
    %     'UniformOutput',false);
    % clim = mean(cat(1,clim{:}),1);
    %-- baseline to 0.25 m/s clims
    % clim = [0,0.0225];
    % clim_ticks = [0,0.005,0.01,0.015,0.020];
    % clim_tlabs = {'0','0.0050','0.010','0.015','0.020'};
    %-- baseline to common
    % clim = [-0.012,0.012];
    % clim_ticks = [-0.012,-0.006,0,0.006,0.012];
    % clim_tlabs = {'-0.012','-0.0060','0','0.0060','0.012'};
    %-- no baseline
    % clim = [0,0.20];
    % clim_ticks = [0,0.05,0.1,0.15,0.20];
    % clim_tlabs = {'0','0.050','0.10','0.15','0.20'};
    %--
    clim = [0.145,0.185];
    clim_ticks = [0.15,0.16,0.17,0.18];
    clim_tlabs = {'0.15','0.16','0.17','0.18'};
    %## INITIATE FIGURE
    ieee_sz = [8.5-(0.65*2),11-(0.7*2)];
    FIGURE_POSITION =ieee_sz;
    AX_DIM = [0.13,0.16];
    % AX_SHIFT = [1.2,-1.5];
    AX_SHIFT = [1.3,-1.1];
    AX_INIT_X = 0.2;
    AX_INIT_Y = 0.75;
    X_DIM = 4;
    fy_shift = 0;
    %--
    x_shift = AX_INIT_X;
    y_shift = AX_INIT_Y;
    x_cnt = 1;
    y_cnt = 1;
    %-- fig
    TITLE_XSHIFT = 0.4;
    TITLE_YSHIFT = 0.975+fy_shift;
    TITLE_BOX_SZ = [0.4,0.4];
    fig = figure('color','white');
    set(fig,'Units','inches', ...
        'Position',[1,1,6.5,9], ...
        'PaperUnits','inches', ...
        'PaperSize',[1 1], ...
        'PaperPosition',[0 0 1 1])
    %-- title
    annotation('textbox',[TITLE_XSHIFT-(TITLE_BOX_SZ(1)/2/2),TITLE_YSHIFT-(TITLE_BOX_SZ(2)/2), ...
            TITLE_BOX_SZ(1),TITLE_BOX_SZ(2)],...
        'String',atlas_name, ...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle', ...
        'LineStyle','none', ...
        'FontName','Arial',...
        'FontSize',14, ...
        'FontWeight','Bold', ...
        'Units','normalized');
    p_sz = get(fig,'Position');
    set(gca,AXES_DEFAULT_PROPS{:});
    %-- axes
    hold on;
    for g_i = 1:length(groups)
        for c_i = 1:length(COND_CHARS)
            %##
            % stat_ext = 'mod';
            tmpst = zeros(size(itco_c{c_i,g_i},1),size(itco_c{c_i,g_i},2));
            % tmpst = pval < 0.05;
            % tmpst = tmpst(:,:,1);

            %##
            tmp = itco_c{c_i,g_i};
            allersp = mean(tmp,3);
            % allersp_mask = mean(tmp,3).*tmpst;
            % allersp_mask = mean(tmp,3);
            % allersp_pcond = double(tmpst);
            allersp_pcond = tmpst;
            %--
            cmp = linspecer();

            %## PLOT
            tmp_plot_struct = PLOT_STRUCT;
            tmp_plot_struct.colormap = linspecer(); %cmp(ceil(length(cmp)/2):length(cmp),:);
            tmp_plot_struct.alpha_multiple = 0.9;
            if g_i == 1
                tmp_plot_struct.title = alltitles{c_i};
                tmp_plot_struct.title_props.FontSize = 10;
            else
                tmp_plot_struct.title = '';
            end
            % tmp_plot_struct.title = c_alltitles{c_i};
            tmp_plot_struct.ax_props.FontSize = 8;
            tmp_plot_struct.ax_props.Box = 'on';
            tmp_plot_struct.ax_props.LineWidth = 2;
            tmp_plot_struct.contourf_grain = 30;
            tmp_plot_struct.ax_props.Position = [x_shift,y_shift,AX_DIM(1),AX_DIM(2)];
            tmp_plot_struct.clim = clim;
            tmp_plot_struct.cbar_ticklabs = clim_tlabs; %{'0','0.04','0.08','0.12'};
            tmp_plot_struct.cbar_ticks = clim_ticks; %[0,0.04,0.08,0.12];
            tmp_plot_struct.cbar_label_shift = [-0.45,0.7];
            tmp_plot_struct.cbar_label = 'ITC';
            
            if c_i > 1
                tmp_plot_struct.do_display_bandmarks = false;
                tmp_plot_struct.ylabel = '';
            end
            if c_i == length(COND_CHARS) && g_i == 1
                tmp_plot_struct.do_display_colorbar = true;
            else
                tmp_plot_struct.do_display_colorbar = false;
            end
            if g_i < length(groups)
                tmp_plot_struct.xticklabel_chars = {''};
                tmp_plot_struct.xlabel = '';
            elseif c_i > 1 && g_i == length(groups)
                tmp_plot_struct.xlabel = '';
            end
            %--
            if c_i == 2
                % GROUPT_CHARS = {'Younger Adults','Older Higher Functioning Adults','Older Lower Functioning Adults'};
                % GROUPT_SHIFT = [0,0.675];
                % GROUPT_PROPS = {...
                %     'LineStyle','none',...
                %     'FontName','Arial', ...
                %     'FontSize',12,...
                %     'FontWeight','bold',...
                %     'HorizontalAlignment','center', ...
                %     'VerticalAlignment','top',...
                %     'Units','normalized'...
                %     };
                % GROUPTITLE_BOXSIZE = [0.5,0.1];
                % xx = 0.5+(-GROUPTITLE_BOXSIZE(1)/2)+GROUPT_SHIFT(1);
                % yy = tmp_plot_struct.ax_props.Position(2)+tmp_plot_struct.ax_props.Position(4)*GROUPT_SHIFT(2);
                % a = annotation(gcf,'textbox',[xx,yy,GROUPTITLE_BOXSIZE],...
                %     'String',GROUPT_CHARS{g_i}, ...
                %     GROUPT_PROPS{:});
                %--
                GROUPT_CHARS = {'YA','OHFA','OLFA'};
                GROUPT_SHIFT = [-0.15,0.55];
                GROUPT_PROPS = {...
                    'LineStyle','none',...
                    'FontName','Arial', ...
                    'FontSize',12,...
                    'FontWeight','bold',...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','top',...
                    'Units','normalized'...
                    };
                GROUPTITLE_BOXSIZE = [0.1,0.1];
                xx = AX_INIT_X+GROUPT_SHIFT(1);
                yy = tmp_plot_struct.ax_props.Position(2)+tmp_plot_struct.ax_props.Position(4)*GROUPT_SHIFT(2);
                a = annotation(gcf,'textbox',[xx,yy,GROUPTITLE_BOXSIZE],...
                    'String',GROUPT_CHARS{g_i}, ...
                    GROUPT_PROPS{:});
            end
            ax = axes();
            [ax] = plot_contourf_cell(ax,allersp,itc_times(tinds),itc_freqs(finds),...
                allersp_pcond, ...
                'PLOT_STRUCT',tmp_plot_struct);
            %## AX SHIFT
            if x_cnt < X_DIM
                x_shift = x_shift + AX_SHIFT(1)*AX_DIM(1);
            else
                y_shift = y_shift + AX_SHIFT(2)*AX_DIM(2);
                x_shift = AX_INIT_X;
                x_cnt = 0;
                y_cnt = y_cnt + 1;
            end
            x_cnt = x_cnt + 1;
        end
    end

    %## ADD STATS
    X_DIM = 5;
    % AX_INIT_X = 0.05;
    AX_INIT_X = 0.135;
    x_shift = AX_INIT_X+AX_DIM(1)/2;
    y_shift = y_shift + (-0.075);
    x_cnt = 1;
    y_cnt = 1;
    % stattitles = {'\Delta Faster Speeds','OHFA-YA','OLFA-YA'};
    % allpvs = {pval(:,:,1),pval(:,:,2),pval(:,:,3)};
    % allest = {estm(:,:,1),estm(:,:,2),estm(:,:,3)};
    
    for s_i = 1:length(allpvs)
        %##
        p_masked = allpvs{s_i}<0.05;
        %## MATLAB BASED FDR
        % FDR()
        % [p_masked,~,~,adjp] = fdr_bh(allpvs{s_i},0.05,'pdep');
        % [adjp,p_masked] = fdr(allpvs{s_i},0.05,'nonparameteric');        
        % [inds,pcrit] = local_fdr(reshape(allpvs{s_i},1,[]),0.05,false);
        % p_masked = allpvs{s_i} < pcrit;
        %##
        allersp = allest{s_i};
        allersp_pcond = allpvs{s_i};
        % allersp_mask = allersp.*p_masked; %(allpvs{s_i}<0.05);
        % allersp_pcond = double(p_masked); %allpvs{s_i}<0.05);
        if all(isnan(allersp),[1,2])
            continue;
        end
        % allersp_mask = allersp;
        % allersp_pcond = zeros(size(allersp));
        %-- clim
        % clim = cellfun(@(x) [prctile(x,5,'all'),prctile(x,95,'all')],{allersp}, ...
        %     'UniformOutput',false);
        % clim = round(mean(cat(1,clim{:}),1),2,'significant');
        % clim_ticks = linspace(clim(1),clim(2),5);
        % clim_tlabs = cellstr(string(clim_ticks));
        %--
        switch clim_asi(s_i)
            case 1
                %-- speed
                clim = [0.010,0.040];
                clim_ticks = [0.010,0.020,0.03,0.04];
                clim_tlabs = cellstr(string(clim_ticks)); %{'0','0.01','0.02','0.03'};
            case 2
                %-- group
                clim = [0,0.009];
                clim_ticks = [0,0.003,0.006,0.009];
                clim_tlabs = cellstr(string(clim_ticks)); %{'0','0.01','0.02','0.03'};
            case 3
                %-- comps
                clim = [-0.012,0];
                clim_ticks = [-0.012,-0.006,0];
                clim_tlabs = cellstr(string(clim_ticks));
            case 4
                %-- a speed
                clim = [0,120];
                clim_ticks = [0,40,80,120];
                clim_tlabs = cellstr(string(clim_ticks));
            case 5
                %-- a group
                clim = [0,4];
                clim_ticks = [0,1,2,3,4];
                clim_tlabs = cellstr(string(clim_ticks));
            otherwise
                %-- clim
                clim = cellfun(@(x) [prctile(x,5,'all'),prctile(x,95,'all')],{allersp}, ...
                    'UniformOutput',false);
                clim = round(mean(cat(1,clim{:}),1),2,'significant');
                clim_ticks = linspace(clim(1),clim(2),5);
                clim_tlabs = cellstr(string(clim_ticks));
        end
        %-- speed
        % clim = [0,0.03];
        % clim_ticks = [0,0.01,0.02,0.03];
        % clim_tlabs = {'0','0.01','0.02','0.03'};
        %-- group
        % clim = [-0.001,0.03];
        % clim_ticks = [0,0.01,0.02,0.03];
        % clim_tlabs = {'0','0.01','0.02','0.03'};

        %##
        tmp_plot_struct = PLOT_STRUCT;
        tmp_plot_struct.colormap = linspecer(); %cmp(ceil(length(cmp)/2):length(cmp),:);
        tmp_plot_struct.alpha_multiple=0.8;
        tmp_plot_struct.title = stattitles{s_i};
        tmp_plot_struct.title_props.FontSize = 10;
        % tmp_plot_struct.title = c_alltitles{c_i};
        
        tmp_plot_struct.ax_props.FontSize = 8;
        tmp_plot_struct.ax_props.Box = 'on';
        tmp_plot_struct.ax_props.LineWidth = 2;
        tmp_plot_struct.contourf_grain = 30;
        tmp_plot_struct.ax_props.Position = [x_shift,y_shift,AX_DIM(1),AX_DIM(2)];
        tmp_plot_struct.clim = clim;
        %-- cbar
        tmp_plot_struct.cbar_label_props.VerticalAlignment = 'bottom';
        tmp_plot_struct.cbar_label_props.Rotation = 0;
        tmp_plot_struct.cbar_props.Location = 'SouthOutside';
        tmp_plot_struct.cbar_tickangle = 45;
        tmp_plot_struct.cbar_shift = [1,0.7];
        tmp_plot_struct.cbar_label_shift = [-0.15,-0.36];
        tmp_plot_struct.cbar_ticklabs = clim_tlabs; %{'0','0.04','0.08','0.12'};
        tmp_plot_struct.cbar_ticks = clim_ticks; %[0,0.04,0.08,0.12];
        tmp_plot_struct.cbar_label = clim_tit{s_i}; %'ITC Estimate';
        tmp_plot_struct.xlabel = '';
        tmp_plot_struct.xticklabel_chars = {'','','','',''};
        if s_i > 1
            tmp_plot_struct.do_display_bandmarks = false;
            tmp_plot_struct.ylabel = '';
        end
        % if s_i < length(allpvs)
        %     tmp_plot_struct.do_display_colorbar = false;
        % end
        % if s_i > 1
        %     tmp_plot_struct.xlabel = '';
        % end
        % tmp_plot_struct.xlabel = '';
        % tmp_plot_struct.xticklabel_chars = {''};
        %--
        if s_i == 1
            GROUPT_SHIFT = [0,0.69];
            GROUPT_PROPS = {...
                'LineStyle','none',...
                'FontName','Arial', ...
                'FontSize',12,...
                'FontWeight','bold',...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','top',...
                'Units','normalized'...
                };
            GROUPTITLE_BOXSIZE = [0.5,0.1];
            xx = 0.5+(-GROUPTITLE_BOXSIZE(1)/2)+GROUPT_SHIFT(1);
            yy = tmp_plot_struct.ax_props.Position(2)+tmp_plot_struct.ax_props.Position(4)*GROUPT_SHIFT(2);
            a = annotation(gcf,'textbox',[xx,yy,GROUPTITLE_BOXSIZE],...
                'String','Statistics', ...
                GROUPT_PROPS{:});
        end
        %##
        ax = axes();
        [ax] = plot_contourf_cell(ax,allersp,itc_times(tinds),itc_freqs(finds),...
            allersp_pcond, ...
            'PLOT_STRUCT',tmp_plot_struct);
        %## AX SHIFT
        if x_cnt < X_DIM
            x_shift = x_shift + AX_SHIFT(1)*AX_DIM(1);
        else
            y_shift = y_shift + AX_SHIFT(2)*AX_DIM(2);
            x_shift = AX_INIT_X;
            x_cnt = 0;
            y_cnt = y_cnt + 1;
        end
        x_cnt = x_cnt + 1;
    end

    drawnow;
    exportgraphics(fig,[tmp_save_dir filesep sprintf('cl%i_%s_%s_%s.png',cl_n,fext,stat_ext,blext)], ...
        'Resolution',300);
    % exportgraphics(fig,[tmp_save_dir filesep sprintf('cl%i_%s_%s_%s.pdf',cl_n,fext,stat_ext,blext)], ...
    %     'ContentType','vector');
    % close(fig);
end
