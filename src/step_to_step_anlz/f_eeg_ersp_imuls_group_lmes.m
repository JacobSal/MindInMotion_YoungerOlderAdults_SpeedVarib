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
clear java
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
        STUDY_DIR = fileparts(SCRIPT_DIR); % change this if in sub folder
        SRC_DIR = fileparts(fileparts(STUDY_DIR));
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
    SRC_DIR = fileparts(fileparts(STUDY_DIR));
end
%## Add Study, Src, & Script Paths
addpath(SCRIPT_DIR)
addpath(SRC_DIR);
addpath(STUDY_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
% [SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca_speed');
%% (PARAMETERS) ======================================================== %%
%- statistics & conditions
speed_trials = {'0p25','0p5','0p75','1p0'};
terrain_trials = {'flat','low','med','high'};
COND_DESIGNS = {speed_trials,terrain_trials};
%-
cmap_terrain = linspecer(4);
custom_yellow = [254,223,0]/255;
cmap_terrain = [cmap_terrain(3,:);custom_yellow;cmap_terrain(4,:);cmap_terrain(2,:)];
cmap_speed = linspecer(4*3);
cmap_speed = [cmap_speed(1,:);cmap_speed(2,:);cmap_speed(3,:);cmap_speed(4,:)];
%* ERSP PARAMS
ERSP_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','cluster',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
SPEC_PARAMS = struct('freqrange',[1,200],...
    'subject','',...
    'specmode','psd',...
    'freqfac',4,...
    'logtrials','on',...
    'comps','all',...
    'plot_freqrange',[4,60],...
    'plot_ylim',[-35,-8],...
    'subtractsubjectmean','on',...
    'plotmode','normal');
ERSP_PARAMS = struct('subbaseline','off',...
    'timerange',[],...
    'ersplim',[-2,2],...
    'freqfac',4,...
    'cycles',[3,0.8],...
    'freqrange',[1,200]);
%% (PATHS)
%- datset name
DATA_SET = 'MIM_dataset';
%- study name
STUDY_DNAME = '10172024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej_speed';
STUDY_FNAME = 'kin_eeg_epoch_study';
ANALYSIS_DNAME = 'kin_eeg_ersp_step_to_step';
%-
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
%- load cluster
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb10'];
% cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb3'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
cluster_k_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
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
%% MEASURES TO ANALYZE ================================================= %%
EEG_MEASURES = {'avg_theta_post_swing2','avg_theta_post_stance2','avg_theta_post_swing1', ...
    'avg_beta_post_swing1','avg_beta_post_stance2', ...
    'avg_alpha_post_swing1','avg_alpha_post_stance2'};
EEG_MEASURES_LABS = {'Mean \theta swing2','Mean \theta stance2','Mean \theta swing1', ...
    'Mean \beta swing1','Mean \beta stance2', ...
    'Mean \alpha swing1','Mean \alpha stance2'};
MEASURE_SYMBOLS = {'\theta','\theta','\theta', ...
    '\beta','\beta', ...
    '\alpha','\alpha'};
%-
% KINNAMES = {'step_dur','var_step_dur_1','var_step_dur_2'};
% KINNAMES_LABS = {'Step Duration (s)','SD_{step}-SD_{mu} (s)','SD_{std}/(sqrt(SD_{step}-SD_{mu})^2)'};
KINNAMES = {'step_width_mm','ap_exc_mm'};
KINNAMES_LABS = {'Step Width (m)','AP Exc. (m)'};
%% (LOAD STATISTICS & DATA EXCEL SHEET FROM R) ========================= %%
%-
% KIN_TABLE = par_load(save_dir,'step_by_step_eeg_ersp_table.mat');
KIN_TABLE = par_load(save_dir,'step_by_step_eeg_ersp_trialbase_table.mat');
for i = 1:length(KINNAMES)
    KIN_TABLE = KIN_TABLE(KIN_TABLE.(KINNAMES{i}) >= 0.01,:);
end
%- 
r_stats_dir = [PATHS.src_dir filesep '2_STUDY' filesep 'mim_yaoa_speed_kin' filesep 'r_scripts' filesep 'sbs_lme_mods'];
RSTATS_IMPORT = readtable([r_stats_dir filesep 'moderation_eeg_ersp_kin_speed_intact_nonorm_0p005cut.xlsx'], ...
    "FileType","spreadsheet","UseExcel",true);
%% ===================================================================== %%
%## SPEED MANUSCRIPT GROUP PLOT
designs = unique(KIN_TABLE.model_n);
% clusters = unique(KIN_TABLE.cluster_n);
groups = unique(KIN_TABLE.group_n);
group_chars = unique(KIN_TABLE.group_char);
cond_chars = unique(KIN_TABLE.cond_char);
speed_ns = unique(KIN_TABLE.speed_n);
clusters = unique(RSTATS_IMPORT.cluster_num);
models = unique(RSTATS_IMPORT.model_char);
% models = [models,'all'];
%##
PLOT_STRUCT = struct('color_map',[],...
    'cond_labels',{{'0.25','0.50','0.75','1.0'}},...
    'cond_offsets',[-0.35,-0.1,0.15,0.40],...
    'group_labels',{{'YA','OHFA','OLFA'}},...
    'group_offsets',[0.125,0.475,0.812],...
    'group_lab_yoffset',-0.23,...
    'group_lab_fontweight','normal',...
    'group_lab_fontsize',10,...
    'y_label','',...
    'y_label_fontsize',10,...
    'y_label_fontweight','bold',...
    'ylim',[],...
    'x_label','Speed (m/s)',...
    'x_label_fontsize',10,...
    'x_label_fontweight','bold',...
    'x_label_yoffset',-0.155,...
    'xlim',[],...
    'title',{{''}},...
    'title_fontsize',10,...
    'title_fontweight','bold',...
    'font_size',9,...
    'font_name','Arial',...
    'do_combine_groups',false,...
    'regresslab_txt_size',7,...
    'ax_position',[0,0,1,1],...
    'ax_line_width',1,...
    'xtick_angle',75);
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],...
        'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
% %-
% cluster_titles = {'','','Right Posterior Parietal','Mid Cingulate','Left Temporal','Right Occipital',...
%     'Right Supplementary Motor','Left Occipital','Left Sensorimotor','Cuneus',...
%     'Left Posterior Parietal',...
%     'Left Supplementary Motor','Right Temporal'};
cluster_titles = {'','','Right Sensorimotor','Mid Cingulate','Left Temporal','Right Occipital',...
    'Right Premotor','Left Occipital','Left Sensorimotor','Right Posterior Parietal',...
    'Left Posterior Parietal',...
    'Left Supplementary Motor','Right Temporal'};
% output_titles = {'','','lsm','rppa','midc','rcun','rsm','lsupm','rocp','locp','ltemp','lppa','rtemp'};
% fig_n = [0,0,7,10,11,10,8,12,0,0,0,9,0];
%-
des_i = 2;
s_chars = {STUDY.datasetinfo.subject};
g_chars = {'Younger Adults','Older High Function','Older Low Function'}; %unique(KIN_TABLE.group_name); %KIN_TABLE.group_char;
% for i = 1:length(STUDY.design(des_i).variable)
%     if strcmp(STUDY.design(des_i).variable(i).label,'cond')
%         c_chars = STUDY.design(des_i).variable(i).value;
%     elseif strcmp(STUDY.design(des_i).variable(i).label,'group')
%         g_chars = STUDY.design(des_i).variable(i).value;
%     end
% end
%## PARAMS
%-
COLORS = cmap_speed;

TITLE_TXT_SIZE = 14;
ALPHA = 0.05;
DO_PLOT_RANDOM_VARS = false;
%-
AX_W = 0.35;
AX_H = 0.25;
AX_TXT_SIZE = 10;
AX_FONT_NAME = 'Arial';
AX_X_SHIFT = 1.2;
AX_Y_SHIFT = 1.55;
AX_INIT_X = 0.1;
AX_INIT_Y = 0.725;
%-
LINE_WIDTH_REFF = 1;
LINE_WIDTH_MEFF = 3;
%-
DO_PLOT_R2 = true;
REG_TXT_SIZE = 7;
REG_X_SHIFT = 0.08;
REG_Y_SHIFT = 0.08;
%-
LEG_X_SHIFT = -0.05;
LEG_Y_SHIFT =  -0.38;
LEG_TXT_SIZE = 9;
LEG_TOKEN_SIZE = 20;
%-
% IM_RESIZE = 0.5;
% X_DIM = 4; % number of measures
% EEG_MEASURES = {'avg_beta_post_swing1','avg_beta_post_stance2', ...
%     'avg_alpha_post_swing1','avg_alpha_post_stance2'};
% EEG_MEASURES_LABS = {'Mean \beta swing1','Mean \beta stance2', ...
%     'Mean \alpha swing1','Mean \alpha stance2'};
% MEASURE_SYMBOLS = {'\beta','\beta', ...
%     '\alpha','\alpha'};
% tmp_savedir = [save_dir filesep 'ersp_kin_beta_alpha'];
%-
IM_RESIZE = 0.7;
X_DIM = 3; % number of measures
EEG_MEASURES = {'avg_theta_post_swing2','avg_theta_post_stance2','avg_theta_post_swing1', ...
    };
EEG_MEASURES_LABS = {'Mean \theta swing2','Mean \theta stance2','Mean \theta swing1', ...
    };
MEASURE_SYMBOLS = {'\theta','\theta','\theta', ...
    };
tmp_savedir = [save_dir filesep 'ersp_kin_theta'];
%-
% IM_RESIZE = 0.3;
% HZ_DIM = 7;
%-
% tmp_savedir = [save_dir filesep '2dplot_intact_group'];
mkdir(tmp_savedir);
%% (PER GROUP MODEL) =================================================== %%
for cl_i = 1:length(clusters)
    %##
    for kin_i = 1:length(KINNAMES)
        %## FIGURE SETUP
        atlas_name = cluster_titles{double(string(clusters(cl_i)))};
        fig = figure('color','white');
        sgtitle(atlas_name,'FontName',AX_FONT_NAME,'FontSize',TITLE_TXT_SIZE,'FontWeight','bold','Interpreter','none');
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        hold on;
        set(gca,AXES_DEFAULT_PROPS{:})
        %-
        x_shift = AX_INIT_X;
        y_cnt = 1;
        x_cnt = 1;
        y_shift = AX_INIT_Y;
        %##
        for grp_i = 1:length(groups)
            vert_shift = 0;
            horiz_shift = 0;
            stats_store = [];
            for meas_i = 1:length(EEG_MEASURES)
                %##
                cond_plot_store = [];
                group_plot_store = [];
                %##
                inds = KIN_TABLE.model_n == des_i &... %designs(des_i) &...
                        KIN_TABLE.cluster_n == clusters(cl_i);
                muc = mean(KIN_TABLE.(KINNAMES{kin_i}));
                stdc = std(KIN_TABLE.(KINNAMES{kin_i}));
                inds_keep = KIN_TABLE.(KINNAMES{kin_i}) > muc-3*stdc & KIN_TABLE.(KINNAMES{kin_i}) < muc+3*stdc;                
                tmp_tbl = KIN_TABLE(inds&inds_keep,:);
                %- determine axes limits
                y_lim_calc = [min(tmp_tbl.(EEG_MEASURES{meas_i}))-std(tmp_tbl.(EEG_MEASURES{meas_i})),...
                    max(tmp_tbl.(EEG_MEASURES{meas_i}))+std(tmp_tbl.(EEG_MEASURES{meas_i}))];
                x_lim_calc = [min(tmp_tbl.(KINNAMES{kin_i}))-std(tmp_tbl.(KINNAMES{kin_i})),...
                    max(tmp_tbl.(KINNAMES{kin_i}))+std(tmp_tbl.(KINNAMES{kin_i}))];
                %##
                if ~isempty(grp_i)
                    inds = KIN_TABLE.model_n == des_i &... %designs(des_i) &...
                        KIN_TABLE.cluster_n == clusters(cl_i) &...
                        KIN_TABLE.group_n == groups(grp_i);
                else
                    inds = KIN_TABLE.model_n == des_i &... %designs(des_i) &...
                        KIN_TABLE.cluster_n == clusters(cl_i);
                end
                %## CUT-OUT KINEMATIC VALUES >/< 3 STD'S FROM MEAN
                muc = mean(KIN_TABLE.(KINNAMES{kin_i}));
                stdc = std(KIN_TABLE.(KINNAMES{kin_i}));
                inds_keep = KIN_TABLE.(KINNAMES{kin_i}) > muc-3*stdc & KIN_TABLE.(KINNAMES{kin_i}) < muc+3*stdc;                
                tmp_tbl = KIN_TABLE(inds&inds_keep,:);                
                %- get speeds & subjects
                speed_ns = unique(tmp_tbl.speed_n);
                subjects = unique(tmp_tbl.subj_char);
                
                %## ASSIGN STATS
                tmp_stats = RSTATS_IMPORT.cluster_num==double(string(clusters(cl_i))) &...
                    strcmp(RSTATS_IMPORT.model_char,'interaction') &...
                    strcmp(RSTATS_IMPORT.freq_band_char,EEG_MEASURES{meas_i}) &...
                    strcmp(RSTATS_IMPORT.kinematic_char,KINNAMES{kin_i}) &...
                    strcmp(RSTATS_IMPORT.group_char,string(group_chars(grp_i)));
                tmp_stats = RSTATS_IMPORT(tmp_stats,:);
                %- intercepts
                % R2 = tmp_stats.r2_m_int;
                F2 = tmp_stats.f2_m_int;
                R2 = tmp_stats.r2_c_int;
                % F2 = tmp_stats.f2_c_int;
                %- speed
                % R2 = tmp_stats.r2_m_nsp;
                % F2 = tmp_stats.f2_m_nsp;
                % R2 = tmp_stats.r2_c_nsp;
                % F2 = tmp_stats.f2_c_nsp;
                %- kin
                % R2 = tmp_stats.r2_m_nkin;
                % F2 = tmp_stats.f2_m_nkin;
                % R2 = tmp_stats.r2_c_nkin;
                % F2 = tmp_stats.f2_c_nkin;
                %- intact
                % R2 = tmp_stats.r2_m_nintact;
                % F2 = tmp_stats.f2_m_nintact;
                % R2 = tmp_stats.r2_c_nintact;
                % F2 = tmp_stats.f2_c_nintact;
                anova_p_kin = tmp_stats.pval_kin;
                anova_p_sp = tmp_stats.pval_sp;
                anova_p_intact = tmp_stats.pval_intact;
                slope_sp = tmp_stats.coeff_sp;
                slope_kin = tmp_stats.coeff_kin;
                slope_intact = tmp_stats.coeff_intact;
                inter_mn = tmp_stats.coeff_int;
                %## SCATTER
                %- plot
                ax = axes();
                hold on;
                for cond_i = 1:length(speed_ns)
                    inds = tmp_tbl.speed_n==speed_ns(cond_i); %double(string(cond_chars(cond_i)));
                    data = tmp_tbl(inds,:);
                    [vals,inds] = sort(data.(KINNAMES{kin_i}));
                    data = data(inds,:);
                    ss = scatter(data,KINNAMES{kin_i},EEG_MEASURES{meas_i},'DisplayName',sprintf('%0.2f',speed_ns(cond_i)));
                    ss.CData = COLORS(cond_i,:);
                    ss.SizeData = 15;
                    % ss.MarkerFaceAlpha = 'flat';
                    % ss.AlphaData = repmat(0.3,[size(data,1),1]);
                    ss.MarkerEdgeAlpha = 0.15; %0.6;
                    if meas_i == 1
                        cond_plot_store = [cond_plot_store, ss];
                    end
                end
                %## LINEAR MODEL FIT
                hold on;
                for cond_i = 1:length(speed_ns)
                    inds = tmp_tbl.speed_n==double(string(speed_ns(cond_i)));
                    data = tmp_tbl(inds,:);
                    [vals,inds] = sort(data.(KINNAMES{kin_i}));
                    data = data(inds,:);
                    if anova_p_intact < ALPHA
                        %## PLOT RANDOM EFFECTS
                        % rnd_colors = linspecer(size(bretable,1));
                        % for subj_i = 1:size(bretable,1)
                        %     xx = unique(data.(KINNAMES{var_i}));
                        %     yy = double(string(speed_ns(cond_i)));
                        %     % y = [];
                        %     % for i = 1:length(x)
                        %     %     y(i) = x(i)*slope_kin + xs*slope_sp + xs*x(i)*slope_intact + bretable.Estimate(subj_i) + inter_mn ;
                        %     % end
                        %     zz = xx.*slope_kin+yy*slope_sp+xx.*(yy*slope_intact)+inter_mn;
                        %     pps = plot(xx,zz,...
                        %         'DisplayName',sprintf('subj. %s',bretable.Level{subj_i}),...
                        %         'LineWidth',LINE_WIDTH_REFF);
                        %     % pp.Color = rnd_colors(subj_i)*0.8;
                        %     pps.Color = [COLORS(cond_i,:)*0.5,0.35];
                        %     pps.LineStyle = ':';
                        % end
                        %## PLOT MAIN EFFECTS
                        xx = unique(data.(KINNAMES{kin_i}));
                        yy = double(string(speed_ns(cond_i)));
                        % zz = [];
                        % for i = 1:length(xx)
                        %     zz(i) = xx(i)*slope_kin + yy*slope_sp + yy*xx(i)*slope_intact + inter_mn ;
                        % end
                        zz = xx.*slope_kin+yy*slope_sp+xx.*(yy*slope_intact)+inter_mn;
                        pp = plot(xx,zz,...
                            'DisplayName',sprintf('pval_{kin,%s}=%0.2f\npval_{speed,%s}=%0.2f\npval_{intact,%s}=%0.2f',...
                            MEASURE_SYMBOLS{meas_i},anova_p_kin,MEASURE_SYMBOLS{meas_i},anova_p_sp,MEASURE_SYMBOLS{meas_i},anova_p_intact),...
                            'LineWidth',LINE_WIDTH_MEFF);
                        pp.Color = COLORS(cond_i,:)*0.8;
                        
                        %## GIVE CORRELATION COEFFICIENT
                        if DO_PLOT_R2
                            chkl = find((anova_p_intact <= [0.05,0.01,0.001]),1,'last');
                            if ~isempty(chkl)
                                ast = repmat('*',[1,chkl]);
                                % str = sprintf('%s R^2=%0.2g\nm_{speed}=%0.2g\nm_{kin}=%0.2g\nm_{sp:kin}=%0.2g',ast,R2,slope_sp,slope_kin,slope_intact);
                                str = sprintf('%s R^2=%0.2g  f^2=%0.2g\ny=%0.1g*sp+%0.1g*kin+%0.1g*sp:kin+%0.1g',ast,R2,F2,slope_sp,slope_kin,slope_intact,inter_mn);
                                
                                r2_pos = [x_shift+REG_X_SHIFT*IM_RESIZE,y_shift+REG_Y_SHIFT*IM_RESIZE,0.2,0.2];
                                annotation('textbox',r2_pos,...
                                        'String',str,'HorizontalAlignment','center',...
                                        'VerticalAlignment','middle','LineStyle','none','FontName',AX_FONT_NAME,...
                                        'FontSize',REG_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                        'BackgroundColor','none');
                            end
                        end
                    end
                end
                
                %## LEGEND
                if meas_i == 1 && grp_i == 1
                    %- lg2
                    legend(gca,cond_plot_store);
                    [lg2,icons,plots,txt]  = legend('boxoff');
                    tmp = get(lg2,'String');
                    cnt = 1;
                    for i = 1:length(cond_plot_store)
                        tmp{i} = sprintf('%0.2g',double(string(speed_ns(cnt))));
                        cnt = cnt + 1;
                    end
                    set(lg2,'String',tmp,'FontName',AX_FONT_NAME,'FontSize',LEG_TXT_SIZE)
                    set(lg2,'Orientation','horizontal')
                    set(lg2,'Units','normalized')
                    set(lg2,'Position',[AX_INIT_X+LEG_X_SHIFT*IM_RESIZE,...
                        AX_INIT_Y+AX_H*IM_RESIZE+LEG_Y_SHIFT*IM_RESIZE,lg2.Position(3),lg2.Position(4)]);
                    lg2.ItemTokenSize(1) = LEG_TOKEN_SIZE;
                end
                %## AX EDITS
                if meas_i == 1
                    ylabel(ax,'10*log_{10}(Flattened PSD)');
                    annotation('textbox',[0.5-(0.1/2),y_shift-(0.1/2)+AX_H*IM_RESIZE,.1,.1],...
                        'String',sprintf('%s',g_chars{grp_i}),'HorizontalAlignment','left',...
                        'VerticalAlignment','top','LineStyle','none','FontName',AX_FONT_NAME,...
                        'FontSize',12,'FontWeight','Bold','Units','normalized');
                else
                    ylabel(ax,'')
                end
                if grp_i == 3
                    xlabel(ax,KINNAMES_LABS{kin_i});
                else
                    xlabel(ax,'')
                end
                if grp_i == 1
                    title(ax,EEG_MEASURES_LABS{meas_i});
                end
                set(ax,'FontWeight','bold','FontSize',AX_TXT_SIZE);
                ylim(ax,y_lim_calc);   
                xlim(ax,x_lim_calc);
                set(ax,'OuterPosition',[0,0,1,1]);
                set(ax,'Position',[x_shift,y_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
                %## AX SHIFT
                if x_cnt < X_DIM
                    x_shift = x_shift + AX_X_SHIFT*IM_RESIZE*AX_W;
                else
                    y_shift = y_shift - AX_Y_SHIFT*IM_RESIZE*AX_H;
                    x_shift = AX_INIT_X;
                    x_cnt = 0;
                end
                x_cnt = x_cnt + 1;
            end
        end
        hold off;
        %##
        exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_kin%s_plot_group.tiff',string(clusters(cl_i)),KINNAMES{kin_i})],'Resolution',300)
        close(fig)
    end
end
%% (ALL SUBJS MODEL) =================================================== %%
%## ALPHA & BETA
IM_RESIZE = 0.5;
X_DIM = 4; % number of measures
EEG_MEASURES = {'avg_beta_post_swing1','avg_beta_post_stance2', ...
    'avg_alpha_post_swing1','avg_alpha_post_stance2'};
EEG_MEASURES_LABS = {'Mean \beta swing1','Mean \beta stance2', ...
    'Mean \alpha swing1','Mean \alpha stance2'};
MEASURE_SYMBOLS = {'\beta','\beta', ...
    '\alpha','\alpha'};
tmp_savedir = [save_dir filesep 'ersp_kin_beta_alpha'];
%## THETA
% IM_RESIZE = 0.7;
% X_DIM = 3; % number of measures
% EEG_MEASURES = {'avg_theta_post_swing2','avg_theta_post_stance2','avg_theta_post_swing1', ...
%     };
% EEG_MEASURES_LABS = {'Mean \theta swing2','Mean \theta stance2','Mean \theta swing1', ...
%     };
% MEASURE_SYMBOLS = {'\theta','\theta','\theta', ...
%     };
% tmp_savedir = [save_dir filesep 'ersp_kin_theta'];
%-
mkdir(tmp_savedir);
grp_i = 1;
for cl_i = 1:length(clusters)
    %##
    for kin_i = 1:length(KINNAMES)
        %##
        atlas_name = cluster_titles{double(string(clusters(cl_i)))};
        fig = figure('color','white');
        sgtitle(atlas_name,'FontName',AX_FONT_NAME,'FontSize',TITLE_TXT_SIZE,'FontWeight','bold','Interpreter','none');
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
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
        for meas_i = 1:length(EEG_MEASURES)
            %##
            cond_plot_store = [];
            group_plot_store = [];
            if ~isempty(grp_i)
                inds = KIN_TABLE.model_n == des_i &... %designs(des_i) &...
                    KIN_TABLE.cluster_n == clusters(cl_i) &...
                    KIN_TABLE.group_n == groups(grp_i);
            else
                inds = KIN_TABLE.model_n == des_i &... %designs(des_i) &...
                    KIN_TABLE.cluster_n == clusters(cl_i);
            end
            %## CUT-OUT KINEMATIC VALUES >/< 3 STD'S FROM MEAN
            muc = mean(KIN_TABLE.(KINNAMES{kin_i}));
            stdc = std(KIN_TABLE.(KINNAMES{kin_i}));
            inds_keep = KIN_TABLE.(KINNAMES{kin_i}) > muc-3*stdc & KIN_TABLE.(KINNAMES{kin_i}) < muc+3*stdc;                
            tmp_tbl = KIN_TABLE(inds&inds_keep,:);
            %- determine axes limits
            y_lim_calc = [min(tmp_tbl.(EEG_MEASURES{meas_i}))-std(tmp_tbl.(EEG_MEASURES{meas_i})),...
                max(tmp_tbl.(EEG_MEASURES{meas_i}))+std(tmp_tbl.(EEG_MEASURES{meas_i}))];
            x_lim_calc = [min(tmp_tbl.(KINNAMES{kin_i}))-std(tmp_tbl.(KINNAMES{kin_i})),...
                max(tmp_tbl.(KINNAMES{kin_i}))+std(tmp_tbl.(KINNAMES{kin_i}))];
            %- get speeds & subjects
            speed_ns = unique(tmp_tbl.speed_n);
            subjects = unique(tmp_tbl.subj_char);
            
            %## ASSIGN STATS
            tmp_stats = RSTATS_IMPORT.cluster_num==double(string(clusters(cl_i))) &...
                strcmp(RSTATS_IMPORT.model_char,'all_interaction') &...
                strcmp(RSTATS_IMPORT.freq_band_char,EEG_MEASURES{meas_i}) &...
                strcmp(RSTATS_IMPORT.kinematic_char,KINNAMES{kin_i}) &...
                strcmp(RSTATS_IMPORT.group_char,'all');
            tmp_stats = RSTATS_IMPORT(tmp_stats,:);
            %- eta/f2
            F2_sp = tmp_stats.fsq_sp;
            F2_kin = tmp_stats.fsq_kin;
            % F2 = tmp_stats.etasq_sp;
            % F2 = tmp_stats.etasq_kin;
            %- intercepts
            % R2 = tmp_stats.r2_m_int;
            % F2 = tmp_stats.f2_m_int;
            R2 = tmp_stats.r2_c_int;
            % F2 = tmp_stats.f2_c_int;
            %-
            anova_p_kin = tmp_stats.pval_kin;
            anova_p_sp = tmp_stats.pval_sp;
            anova_p_intact = tmp_stats.pval_intact;
            slope_sp = tmp_stats.coeff_sp;
            slope_kin = tmp_stats.coeff_kin;
            slope_intact = tmp_stats.coeff_intact;
            inter_mn = tmp_stats.coeff_int;
            %## SCATTER
            %- plot
            ax = axes();
            hold on;
            for cond_i = 1:length(speed_ns)
                inds = tmp_tbl.speed_n==speed_ns(cond_i); %double(string(cond_chars(cond_i)));
                data = tmp_tbl(inds,:);
                [vals,inds] = sort(data.(KINNAMES{kin_i}));
                data = data(inds,:);
                ss = scatter(data,KINNAMES{kin_i},EEG_MEASURES{meas_i},'DisplayName',sprintf('%0.2f',speed_ns(cond_i)));
                ss.CData = COLORS(cond_i,:);
                ss.SizeData = 1;
                % ss.MarkerFaceAlpha = 'flat';
                % ss.AlphaData = repmat(0.3,[size(data,1),1]);
                ss.MarkerEdgeAlpha = 0.6; %0.6;
                if meas_i == 1
                    cond_plot_store = [cond_plot_store, ss];
                end
            end
            %## LINEAR MODEL FIT
            hold on;
            for cond_i = 1:length(speed_ns)
                inds = tmp_tbl.speed_n==double(string(speed_ns(cond_i)));
                data = tmp_tbl(inds,:);
                [vals,inds] = sort(data.(KINNAMES{kin_i}));
                data = data(inds,:);
                if anova_p_intact < ALPHA
                    %## PLOT RANDOM EFFECTS
                    % if DO_PLOT_RANDOM_VARS
                    %     rnd_colors = linspecer(size(bretable,1));
                    %     for subj_i = 1:size(bretable,1)
                    %         x = unique(data.(KINNAMES{var_i}));
                    %         xs = double(string(speed_ns(cond_i)));
                    %         y = [];
                    %         for i = 1:length(x)
                    %             y(i) = x(i)*slope_kin + xs*slope_sp + xs*x(i)*slope_intact + bretable.Estimate(subj_i) + inter_mn ;
                    %         end
                    %         pps = plot(x,y,...
                    %             'DisplayName',sprintf('subj. %s',bretable.Level{subj_i}),...
                    %             'LineWidth',LINE_WIDTH_REFF);
                    %         % pp.Color = rnd_colors(subj_i)*0.8;
                    %         pps.Color = [COLORS(cond_i,:)*0.5,0.35];
                    %         pps.LineStyle = ':';
                    %     end
                    % end
                    %## PLOT MAIN EFFECTS
                    xx = sort([unique(data.(KINNAMES{kin_i}));x_lim_calc']);
                    yy = double(string(speed_ns(cond_i)));
                    zz = [];
                    for i = 1:length(xx)
                        zz(i) = xx(i)*slope_kin + yy*slope_sp + yy*xx(i)*slope_intact + inter_mn;
                    end
                    pp = plot(xx,zz,...
                        'DisplayName',sprintf('LME Fit %0.2f (m/s)',speed_ns(cond_i)),...
                        'LineWidth',2);
                    pp.Color = [COLORS(cond_i,:)*0.8,0.6];
                    %## GIVE CORRELATION COEFFICIENT
                    if DO_PLOT_R2
                        chkl = find((anova_p_intact <= [0.05,0.01,0.001]),1,'last');
                        if ~isempty(chkl)
                            ast = repmat('*',[1,chkl]);
                            % str = sprintf('%s R^2=%0.2g\nm_{speed}=%0.2g\nm_{kin}=%0.2g\nm_{sp:kin}=%0.2g',ast,R2,slope_sp,slope_kin,slope_intact);
                            % str = sprintf('%s R^2=%0.2g  f^2=%0.2g\ny=%0.1g*sp+%0.1g*kin+%0.1g*sp:kin+%0.1g',ast,R2,F2,slope_sp,slope_kin,slope_intact,inter_mn);
                            % str = sprintf('%s R^2=%0.2g  f^2=%0.2g\ny=(%0.1g)sp+(%0.1g)kin+(%0.1g)sp*kin+(%0.1g)',ast,R2,F2,slope_sp,slope_kin,slope_intact,inter_mn);
                            str = sprintf('%s R^2=%0.2g\nf_{speed}^2=%0.2g  f_{kin}^2=%0.2g',ast,R2,F2_sp,F2_kin);
                            
                            r2_pos = [x_shift+REG_X_SHIFT*IM_RESIZE-(0.3/4),y_shift+REG_Y_SHIFT*IM_RESIZE-(0.3/4),0.3,0.3];
                            annotation('textbox',r2_pos,...
                                    'String',str,'HorizontalAlignment','center',...
                                    'VerticalAlignment','middle','LineStyle','none','FontName',AX_FONT_NAME,...
                                    'FontSize',REG_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                    'BackgroundColor','none');
                        end
                    end
                end
            end
            
            %## LEGEND
            if meas_i == 1
                %- lg2                
                legend(gca,cond_plot_store);
                [lg2,icons,plots,txt]  = legend('boxoff');
                tmp = get(lg2,'String');
                cnt = 1;
                for i = 1:length(cond_plot_store)
                    tmp{i} = sprintf('%0.2g',double(string(speed_ns(cnt))));
                    cnt = cnt + 1;
                end
                set(lg2,'String',tmp,'FontName',AX_FONT_NAME,'FontSize',LEG_TXT_SIZE)
                set(lg2,'Orientation','horizontal')
                set(lg2,'Units','normalized')
                set(lg2,'Position',[AX_INIT_X+LEG_X_SHIFT*IM_RESIZE,...
                    AX_INIT_Y+AX_H*IM_RESIZE+LEG_Y_SHIFT*IM_RESIZE,lg2.Position(3),lg2.Position(4)]);
                lg2.ItemTokenSize(1) = LEG_TOKEN_SIZE;
            end
            %## AX EDITS
            if meas_i == 1
                ylabel(ax,'10*log_{10}(Flattened PSD)');
            else
                ylabel(ax,'')
            end            
            title(ax,EEG_MEASURES_LABS{meas_i});
            xlabel(ax,KINNAMES_LABS{kin_i});
            set(ax,'FontWeight','bold','FontSize',AX_TXT_SIZE);
            ylim(ax,y_lim_calc);   
            xlim(ax,x_lim_calc);
            set(ax,'OuterPosition',[0,0,1,1]);
            set(ax,'Position',[x_shift,y_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
            %## AX SHIFT
            if x_cnt < X_DIM
                x_shift = x_shift + AX_X_SHIFT*IM_RESIZE*AX_W;
            else
                y_shift = y_shift - AX_Y_SHIFT*IM_RESIZE*AX_H;
                x_shift = AX_INIT_X;
                x_cnt = 0;
            end
            x_cnt = x_cnt + 1;
        end
        hold off;
        %##
        % exportgraphics(fig,[tmp_savedir filesep sprintf('intacteff_cl%s_kin%s_plot_allsubj.pdf',string(clusters(cl_i)),KINNAMES{kin_i})],'ContentType','vector')
        exportgraphics(fig,[tmp_savedir filesep sprintf('intacteff_cl%s_kin%s_plot_allsubj.tiff',string(clusters(cl_i)),KINNAMES{kin_i})],'Resolution',900)
        close(fig)
    end
end
%% ===================================================================== %%
%##
PLOT_STRUCT = struct('color_map',[],...
    'cond_labels',{{'0.25','0.50','0.75','1.0'}},...
    'cond_offsets',[-0.35,-0.1,0.15,0.40],...
    'group_labels',{{'YA','OHFA','OLFA'}},...
    'group_offsets',[0.125,0.475,0.812],...
    'group_lab_yoffset',-0.23,...
    'group_lab_fontweight','normal',...
    'group_lab_fontsize',10,...
    'y_label','',...
    'y_label_fontsize',10,...
    'y_label_fontweight','bold',...
    'ylim',[],...
    'x_label','Speed (m/s)',...
    'x_label_fontsize',10,...
    'x_label_fontweight','bold',...
    'x_label_yoffset',-0.155,...
    'z_label','eeg',...
    'z_label_fontsize',10,...
    'z_label_fontweight','bold',...
    'z_label_yoffset',-0.155,...
    'xlim',[],...
    'title',{{''}},...
    'title_fontsize',10,...
    'title_fontweight','bold',...
    'font_size',9,...
    'font_name','Arial',...
    'do_combine_groups',false,...
    'regresslab_txt_size',7,...
    'ax_position',[0,0,1,1],...
    'ax_line_width',1,...
    'xtick_angle',75);
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],...
        'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};

%## PARAMS
%-
COLORS = cmap_speed;
IM_RESIZE = 0.7;
TITLE_TXT_SIZE = 14;
ALPHA = 0.05;
DO_PLOT_RANDOM_VARS = false;
%-
AX_W = 0.35;
AX_H = 0.25;
AX_TXT_SIZE = 10;
AX_FONT_NAME = 'Arial';
AX_X_SHIFT = 1.2;
AX_Y_SHIFT = 1.55;
AX_INIT_X = 0.1;
AX_INIT_Y = 0.66;
%-
LINE_WIDTH_REFF = 1;
LINE_WIDTH_MEFF = 3;
%-
DO_PLOT_R2 = true;
REG_TXT_SIZE = 6;
REG_X_SHIFT = 0.08;
REG_Y_SHIFT = 0.12;
%-
LEG_X_SHIFT = 0;
LEG_Y_SHIFT =  0.1;
LEG_TXT_SIZE = 9;
LEG_TOKEN_SIZE = 20;
%-
VIEW_3D = [45,20];
GRID_N = 5; 
%-
X_DIM = 3;
tmp_savedir = [save_dir filesep '3dplot_intact'];
mkdir(tmp_savedir);
%%
for cl_i = 1:length(clusters)
    %##
    for kin_i = 1:length(KINNAMES)
        %##
        atlas_name = cluster_titles{double(string(clusters(cl_i)))};
        fig = figure('color','white');
        sgtitle(atlas_name,'FontName',AX_FONT_NAME,'FontSize',TITLE_TXT_SIZE,'FontWeight','bold','Interpreter','none');
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        hold on;
        set(gca,AXES_DEFAULT_PROPS{:})
        %-
        x_shift = AX_INIT_X;
        x_cnt = 1;
        y_shift = AX_INIT_Y;
        %##
        for grp_i = 1:length(groups)
            vert_shift = 0;
            horiz_shift = 0;
            stats_store = [];
            for meas_i = 1:length(EEG_MEASURES)
                %##
                cond_plot_store = [];
                group_plot_store = [];
                if ~isempty(grp_i)
                    inds = KIN_TABLE.model_n == des_i &... %designs(des_i) &...
                        KIN_TABLE.cluster_n == clusters(cl_i) &...
                        KIN_TABLE.group_n == groups(grp_i);
                else
                    inds = KIN_TABLE.model_n == des_i &... %designs(des_i) &...
                        KIN_TABLE.cluster_n == clusters(cl_i);
                end
                %## CUT-OUT KINEMATIC VALUES >/< 3 STD'S FROM MEAN
                muc = mean(KIN_TABLE.(KINNAMES{kin_i}));
                stdc = std(KIN_TABLE.(KINNAMES{kin_i}));
                inds_keep = KIN_TABLE.(KINNAMES{kin_i}) > muc-3*stdc & KIN_TABLE.(KINNAMES{kin_i}) < muc+3*stdc;                
                tmp_tbl = KIN_TABLE(inds & inds_keep,:);
                %- determine axes limits
                z_lim_calc = [min(tmp_tbl.(EEG_MEASURES{meas_i}))-std(tmp_tbl.(EEG_MEASURES{meas_i})),...
                    max(tmp_tbl.(EEG_MEASURES{meas_i}))+std(tmp_tbl.(EEG_MEASURES{meas_i}))];
                x_lim_calc = [min(tmp_tbl.(KINNAMES{kin_i}))-std(tmp_tbl.(KINNAMES{kin_i})),...
                    max(tmp_tbl.(KINNAMES{kin_i}))+std(tmp_tbl.(KINNAMES{kin_i}))];
                %-
                speed_ns = unique(tmp_tbl.speed_n);
                subjects = unique(tmp_tbl.subj_char);

                %## ASSIGN STATS
                tmp_stats = RSTATS_IMPORT.cluster_num==double(string(clusters(cl_i))) &...
                    strcmp(RSTATS_IMPORT.model_char,'interaction') &...
                    strcmp(RSTATS_IMPORT.freq_band_char,EEG_MEASURES{meas_i}) &...
                    strcmp(RSTATS_IMPORT.kinematic_char,KINNAMES{kin_i}) &...
                    strcmp(RSTATS_IMPORT.group_char,string(group_chars(grp_i)));
                tmp_stats = RSTATS_IMPORT(tmp_stats,:);
                % R2 = tmp_stats.r2_m_int;
                % F2 = tmp_stats.f2_m_int;
                R2 = tmp_stats.r2_c_int;
                F2 = tmp_stats.f2_c_int;
                anova_p_kin = tmp_stats.pval_kin;
                anova_p_sp = tmp_stats.pval_sp;
                anova_p_intact = tmp_stats.pval_intact;
                slope_sp = tmp_stats.coeff_sp;
                slope_kin = tmp_stats.coeff_kin;
                slope_intact = tmp_stats.coeff_intact;
                inter_mn = tmp_stats.coeff_int;
                %## SCATTER
                %- plot
                ax = axes();
                hold on;
                for cond_i = 1:length(speed_ns)
                    inds = tmp_tbl.speed_n==double(string(speed_ns(cond_i)));
                    data = tmp_tbl(inds,:);
                    [vals,inds] = sort(data.(KINNAMES{kin_i}));
                    xx = vals;
                    yy = data.speed_n(inds);
                    zz = data.(EEG_MEASURES{meas_i})(inds);
                    ss = scatter3(xx,yy,zz);
                    ss.CData = COLORS(cond_i,:);
                    ss.SizeData = 15;
                    ss.MarkerEdgeAlpha = 0.8;
                    if meas_i == 1
                        cond_plot_store = [cond_plot_store, ss];
                    end
                end
                %## LINEAR MODEL FIT
                hold on;
                if anova_p_intact < ALPHA
                    str = sprintf('%s R^2=%0.2g  f^2=%0.2g\ny=%0.1g*sp+%0.1g*kin+%0.1g*sp:kin+%0.1g',ast,R2,F2,slope_sp,slope_kin,slope_intact,inter_mn);
                    %## GRID FITS
                    data = tmp_tbl;
                    xx = linspace(min(data.(KINNAMES{kin_i})),max(data.(KINNAMES{kin_i})),GRID_N);
                    yy = linspace(min(data.speed_n),max(data.speed_n),GRID_N);
                    zz = xx*slope_kin + yy'*slope_sp + yy'*xx*slope_intact + inter_mn;
                    pp = surf(xx,yy,zz,'DisplayName',str);
                    pp.CData = zz;
                    pp.MarkerFaceColor = COLORS(cond_i,:);
                    pp.FaceAlpha = 0;
                    pp.EdgeAlpha = 0.5;
                    % colormap(ax,linspecer)
                    %## LINE FITS
                    % for cond_i = 1:length(speed_ns)
                    %     inds = tmp_tbl.speed_n==double(string(speed_ns(cond_i)));
                    %     data = tmp_tbl(inds,:);
                    %     [vals,inds] = sort(data.(KINNAMES{kin_i}));
                    %     data = data(inds,:);
                    %     %## PLOT MAIN EFFECTS
                    %     xx = min(data.(KINNAMES{kin_i})):0.5:max(data.(KINNAMES{kin_i}));
                    %     yy = linspace(min(data.speed_n),max(data.speed_n),length(xx));
                    %     zz = xx*slope_kin + yy*slope_sp + yy.*xx*slope_intact + inter_mn ;
                    %     % pp = surf(xx,yy,zz,'DisplayName',str);
                    %     pp = plot3(xx,yy,zz);
                    %     pp.Color = COLORS(cond_i,:)*0.8;
                    %     pp.LineWidth = 3;
                    % end
                    %## GIVE CORRELATION COEFFICIENT
                    if DO_PLOT_R2
                        chkl = find((anova_p_intact <= [0.05,0.01,0.001]),1,'last');
                        if ~isempty(chkl)
                            ast = repmat('*',[1,chkl]);
                            % str = sprintf('%s R^2=%0.2g\nm_{speed}=%0.2g\nm_{kin}=%0.2g\nm_{sp:kin}=%0.2g',ast,R2,slope_sp,slope_kin,slope_intact);
                            % str = sprintf('%s R^2=%0.2g  f^2=%0.2g\ny=%0.1g*sp+%0.1g*kin+%0.1g*sp:kin+%0.1g',ast,R2,F2,slope_sp,slope_kin,slope_intact,inter_mn);
                            r2_pos = [x_shift+REG_X_SHIFT*IM_RESIZE,y_shift+REG_Y_SHIFT*IM_RESIZE,0.2,0.2];
                            annotation('textbox',r2_pos,...
                                    'String',str,'HorizontalAlignment','center',...
                                    'VerticalAlignment','middle','LineStyle','none','FontName',AX_FONT_NAME,...
                                    'FontSize',REG_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                    'BackgroundColor','none');
                        end
                    end
                end
                
                %## LEGEND
                if meas_i == 1 && grp_i == 1
                    %- lg2
                    legend(gca,cond_plot_store);
                    [lg2,icons,plots,txt]  = legend('boxoff');
                    tmp = get(lg2,'String');
                    cnt = 1;
                    for i = 1:length(cond_plot_store)
                        tmp{i} = sprintf('%0.2g',double(string(speed_ns(cnt))));
                        cnt = cnt + 1;
                    end
                    set(lg2,'String',tmp,'FontName',AX_FONT_NAME,'FontSize',LEG_TXT_SIZE)
                    set(lg2,'Orientation','horizontal')
                    set(lg2,'Units','normalized')
                    set(lg2,'Position',[AX_INIT_X+LEG_X_SHIFT*IM_RESIZE,...
                        AX_INIT_Y+AX_H*IM_RESIZE+LEG_Y_SHIFT*IM_RESIZE,lg2.Position(3),lg2.Position(4)]);
                    lg2.ItemTokenSize(1) = LEG_TOKEN_SIZE;
                end
                %## AX EDITS
                set(ax,'View',VIEW_3D)
                if meas_i == 1
                    zlabel(ax,'10*log_{10}(Flattened PSD)',...
                        'FontSize',PLOT_STRUCT.z_label_fontsize,...
                        'FontWeight',PLOT_STRUCT.z_label_fontweight)
                    %- group label
                    annotation('textbox',[0.5-(0.1/2),y_shift-(0.1/2)+AX_H*IM_RESIZE,.1,.1],...
                        'String',sprintf('%s',g_chars{grp_i}),'HorizontalAlignment','left',...
                        'VerticalAlignment','top','LineStyle','none','FontName',AX_FONT_NAME,...
                        'FontSize',12,'FontWeight','Bold','Units','normalized');
                else
                    zlabel(ax,'')
                end
                if grp_i == 3 && meas_i == 1
                    ylabel(ax,'Speed (m/s',...
                        'FontSize',PLOT_STRUCT.y_label_fontsize, ...
                        'FontWeight',PLOT_STRUCT.y_label_fontweight);
                    xlabel(ax,KINNAMES_LABS{kin_i},...
                        'FontSize',PLOT_STRUCT.x_label_fontsize,...
                        'FontWeight',PLOT_STRUCT.x_label_fontweight);
                else
                    xlabel(ax,'')
                end
                if grp_i == 1
                    title(ax,EEG_MEASURES_LABS{meas_i});
                end
                set(ax,'FontWeight','bold','FontSize',AX_TXT_SIZE);
                ylim(ax,[0,1.25]);   
                xlim(ax,x_lim_calc);
                zlim(ax,z_lim_calc);  
                set(ax,'OuterPosition',[0,0,1,1]);
                set(ax,'Position',[x_shift,y_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
                set(ax,'projection', 'perspective', 'box', 'on')
                %-
                % h = rotate3d;
                % set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
                % set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
                % set(gcf, 'ResizeFcn', @align_axislabel)
                align_axislabel([], ax)
                %## AX SHIFT
                if x_cnt < X_DIM
                    x_shift = x_shift + AX_X_SHIFT*IM_RESIZE*AX_W;
                else
                    y_shift = y_shift - AX_Y_SHIFT*IM_RESIZE*AX_H;
                    x_shift = AX_INIT_X;
                    x_cnt = 0;
                end
                x_cnt = x_cnt + 1;
            end
        end
        hold off;
        %##
        exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_kin%s_plot_3dplot.tiff',string(clusters(cl_i)),KINNAMES{kin_i})],'Resolution',300)
        close(fig)
    end
end