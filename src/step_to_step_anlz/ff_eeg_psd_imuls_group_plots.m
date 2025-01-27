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
global ADD_CLEANING_SUBMODS STUDY_DIR SCRIPT_DIR %#ok<GVMIS>
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
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
%- study name
STUDY_DNAME = '10172024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej_speed';
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
% cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb3'];
cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb5'];
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
%## RAW STEP-BY-STEP DATA
% KIN_TABLE = par_load(save_dir,'sbs_eeg_psd_meandesignb.mat');
% %- chk
% subj_chars = unique(KIN_TABLE.subj_char);
% strides = zeros(length(subj_chars),1);
% for i = 1:length(subj_chars)
%     tmp = KIN_TABLE(strcmp(subj_chars{i},KIN_TABLE.subj_char),:);
%     strides(i) = max(tmp.stride_n);
% end
% fprintf('mean number of strides: %0.2f\n',mean(strides));
% fprintf('std number of strides: %0.2f\n',std(strides));
% %- r-stats
% RSTATS_IMPORT = readtable([r_stats_dir filesep 'lme_eeg_kin_raw_stats_meandesignb.xlsx'], ...
%     "FileType","spreadsheet","UseExcel",true);
% X_DIM = 1;

%## MEAN & SD TABLE
KIN_TABLE = readtable([r_stats_dir filesep 'lme_eeg_kin_mean_sd_tbl_meandesignb.xlsx'], ...
    "FileType","spreadsheet", ...
    "UseExcel",true);
RSTATS_IMPORT = readtable([r_stats_dir filesep 'lme_eeg_kin_mean_sd_stats_meandesignb.xlsx'], ...
    "FileType","spreadsheet","UseExcel",true);
X_DIM = 2;
%% MEASURES TO ANALYZE ================================================= %%
%## MEAN SD MEASURES
% EEG_MEASURES = {'avg_theta_post_fn1','avg_theta_post_fn2', ...
%     'avg_alpha_post_fn1','avg_alpha_post_fn2', ...
%     'avg_theta_post_fn1','avg_theta_post_fn2'};
% EEG_MEASURES_LABS = {'Mean \theta','Std. Dev. \theta', ...
%     'Mean \alpha','Std. Dev. \alpha', ...
%     'Mean \beta','Std. Dev. \beta'};
% MEASURE_SYMBOLS = {'\mu \theta','\sigma \theta', ...
%     '\mu \alpha','\sigma \alpha', ...
%     '\mu \beta','\sigma \beta'};
EEG_MEASURES = {'avg_theta_post_fn1','avg_theta_post_fn2'};
EEG_MEASURES_LABS = {'Mean \theta','Std. Dev. \theta'};
MEASURE_SYMBOLS = {'\mu \theta','\sigma \theta'};
%-
KINNAMES = {'ml_exc_mm_gc_fn1','ml_exc_mm_gc_fn2'};
KINNAMES_LABS = {'Mean ML Excursion (m)','Std. ML Excursion (m)'};
tmp_savedir = [save_dir filesep 'mean_std_group_kin_eeg'];

%## RAW STEP MEASURES
% EEG_MEASURES = {'avg_theta_post','avg_alpha_post','avg_beta_post'};
% EEG_MEASURES_LABS = {'Mean \theta_{i}','Mean \alpha_{i}','Mean \beta_{i}'};
% MEASURE_SYMBOLS = {'\theta','\alpha','\beta'};
% %-
% KINNAMES = {'ml_exc_mm_gc'};
% KINNAMES_LABS = {'ML Excursion (m)'};
% tmp_savedir = [save_dir filesep 'raw_group_kin_eeg'];
%% ===================================================================== %%
%## SPEED MANUSCRIPT GROUP PLOT

%## (01/16/2025) __iclabel_kmeansalt_rb10
% cluster_titles = {'','','Right Sensorimotor','Mid Cingulate','Left Temporal','Right Occipital',...
%     'Right Premotor','Left Occipital','Left Sensorimotor','Right Posterior Parietal',...
%     'Left Posterior Parietal',...
%     'Left Supplementary Motor','Right Temporal'};
%## (01/22/2025) __iclabel_kmeasalt_rb5
cluster_titles = {'','','Anterior Cingulate','Right Sensorimotor','Left Occipital','Left Supplementary Motor',...
    'Mid Cingulate','Precuneus','Left Sensorimotor','Left Posterior Parietal',...
    'Right Temporal',...
    'Right Occipital','Right Posterior Parietal'};
%-
% desdes = cat(1,STUDY.design.variable);
% c_chars_d = desdes(strcmp({desdes.label},'cond'));
% c_chars_d = {c_chars_d.value};
% g_chars_d = desdes(strcmp({desdes.label},'group'));
% g_chars_d = {g_chars_d.value};

%## PLOT PARAMS
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],...
        'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
%-
COLORS = cmap_speed;
ALPHA = 0.05;
DO_PLOT_RANDOM_VARS = false;
%-
% X_DIM = 2;
IM_RESIZE = 0.675;
AX_W = 0.45;
AX_H = 0.3;
AX_TXT_SIZE = 10; % 0.7
AX_FONT_NAME = 'Arial';
AX_X_SHIFT = 1.5;
AX_Y_SHIFT = 1.55;
AX_INIT_X = 0.125;
AX_INIT_Y = 0.7;
%-
TITLE_FONT_SIZE = 14;
TITLE_XSHIFT = 0.4;
TITLE_YSHIFT = 0.975;
TITLE_BOX_SZ = [0.4,0.4];
%-
LINE_WIDTH_REFF = 1; %#ok<NASGU>
LINE_WIDTH_MEFF = 2;
%-
DO_PLOT_STAT_STR = true;
REG_TXT_SIZE = 8; % 7
REG_X_SHIFT = 0.325; % 0.08
REG_Y_SHIFT = 0.1; % 0.1k
%-
LEG_X_SHIFT = -0.125; %-0.1
LEG_Y_SHIFT =  -0.4; %-0.38
LEG_TXT_SIZE = 9;
LEG_TOKEN_SIZE = 15;


%## DATA PARAMS
des_i = 2;
% grp_i = [];
s_chars = {STUDY.datasetinfo.subject};
g_chars = {'Younger Adults','Older High Function','Older Low Function'}; %unique(KIN_TABLE.group_name); %KIN_TABLE.group_char;
designs = unique(KIN_TABLE.model_n);
% clusters = unique(KIN_TABLE.cluster_n);
group_chars = unique(KIN_TABLE.group_char);
cond_chars = unique(KIN_TABLE.cond_char);
speed_ns = unique(KIN_TABLE.speed_n);
r_clusters = unique(RSTATS_IMPORT.cluster_num);
r_groups = unique(RSTATS_IMPORT.group_char);
% models = unique(RSTATS_IMPORT.model_char);
%% (ALL SUBJS MODEL) =================================================== %%
%-
mkdir(tmp_savedir);
% grp_i = 1;
KIN_CUTOFF = 0.005;
for cl_i = 1:length(r_clusters)    
    %##
    for kin_i = 1:length(KINNAMES)
        %## INITIATE FIGURE
        atlas_name = cluster_titles{double(string(r_clusters(cl_i)))};
        fig = figure('color','white');
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        annotation('textbox',[TITLE_XSHIFT-(TITLE_BOX_SZ(1)/2/2),TITLE_YSHIFT-(TITLE_BOX_SZ(2)/2),TITLE_BOX_SZ(1),TITLE_BOX_SZ(2)],...
            'String',atlas_name, ...
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

        %## LOOP THROUGH GROUPS AND MEASURES
        vert_shift = 0;
        for grp_i = 1:3 %length(group_chars)               
            horiz_shift = 0;
            stats_store = [];                 
            for meas_i = 1:length(EEG_MEASURES) 
                %## REMOVE BAD STEPS
                tmp_kint = KIN_TABLE(KIN_TABLE.(KINNAMES{kin_i}) >= KIN_CUTOFF,:);
                %## ORIGINAL DATA EXTRACTION
                inds = strcmp(tmp_kint.group_char,group_chars{grp_i}) &...
                    tmp_kint.model_n == des_i &... 
                    tmp_kint.cluster_n == r_clusters(cl_i);
    
                %## MEAD SD TABLE EXTRACTION
                % if ~isempty(grp_i)
                %     inds = tmp_kint.cluster_n == clusters(cl_i) &...
                %         strcmp(tmp_kint.group_char,group_chars{grp_i});
                % else
                %     inds = tmp_kint.cluster_n == clusters(cl_i);
                % end
    
                %## CUT-OUT KINEMATIC VALUES >/< 3 STD'S FROM MEAN
                tmp_kint = tmp_kint(inds,:);
                % inds_keep = (tmp_kint.(KINNAMES{kin_i}) > KIN_CUTOFF);
                % tmp_kint = tmp_kint(inds_keep,:);
                muc = mean(tmp_kint.(KINNAMES{kin_i}));
                stdc = std(tmp_kint.(KINNAMES{kin_i}));
                inds_keep = tmp_kint.(KINNAMES{kin_i}) > muc-3*stdc & tmp_kint.(KINNAMES{kin_i}) < muc+3*stdc;                
                tmp_tbl = tmp_kint(inds_keep,:);       
                %##
                cond_plot_store = [];
                group_plot_store = [];
                % %- determine axes limits
                % y_lim_calc = [min(tmp_tbl.(EEG_MEASURES{meas_i}))-0.25*std(tmp_tbl.(EEG_MEASURES{meas_i})),...
                %     max(tmp_tbl.(EEG_MEASURES{meas_i}))+0.25*std(tmp_tbl.(EEG_MEASURES{meas_i}))];
                % % y_lim_calc = [mean(tmp_tbl.(EEG_MEASURES{meas_i}))-3.5*std(tmp_tbl.(EEG_MEASURES{meas_i})),...
                % %     mean(tmp_tbl.(EEG_MEASURES{meas_i}))+3.5*std(tmp_tbl.(EEG_MEASURES{meas_i}))];
                % x_lim_calc = [min(tmp_tbl.(KINNAMES{kin_i}))-std(tmp_tbl.(KINNAMES{kin_i})),...
                %     max(tmp_tbl.(KINNAMES{kin_i}))+std(tmp_tbl.(KINNAMES{kin_i}))];
                %- get speeds & subjects
                speed_ns = unique(tmp_tbl.speed_n);
                subjects = unique(tmp_tbl.subj_char);
                
                %## ASSIGN STATS
                tmp_stats = RSTATS_IMPORT.cluster_num==double(string(r_clusters(cl_i))) &...
                    strcmp(RSTATS_IMPORT.model_char,'interaction_group') &...
                    strcmp(RSTATS_IMPORT.freq_band_char,EEG_MEASURES{meas_i}) &...
                    strcmp(RSTATS_IMPORT.kinematic_char,KINNAMES{kin_i}) &...
                    strcmp(RSTATS_IMPORT.group_char,r_groups{grp_i});
                tmp_stats = RSTATS_IMPORT(tmp_stats,:);
                % tmp_stats = tmp_stats(4,:);
                %-
                ran_effs_char = strsplit(tmp_stats.ran_effs_char{1},',');
                ran_effs_n = cellfun(@(x) double(string(x)),strsplit(tmp_stats.ran_effs_n{1},','));
    
                %## NORMALIZE DATA USING SUBJECT INTERCEPTS
                subj_chars = unique(tmp_tbl.subj_char);
                for subj_i = 1:length(subj_chars)
                    ind = strcmp(subj_chars{subj_i},ran_effs_char);
                    int = ran_effs_n(ind);
                    ind = strcmp(subj_chars{subj_i},tmp_tbl.subj_char);
                    tmp_tbl(ind,EEG_MEASURES{meas_i}) = tmp_tbl(ind,EEG_MEASURES{meas_i})-int;
                end
    
                %## DETERMINE AXES LIMITS & EXTRACT STATS
                % y_lim_calc = [mean(tmp_tbl.(EEG_MEASURES{meas_i}))-3.5*std(tmp_tbl.(EEG_MEASURES{meas_i})),...
                %     mean(tmp_tbl.(EEG_MEASURES{meas_i}))+3.5*std(tmp_tbl.(EEG_MEASURES{meas_i}))];
                y_lim_calc = [min(tmp_tbl.(EEG_MEASURES{meas_i}))-1*std(tmp_tbl.(EEG_MEASURES{meas_i})),...
                    max(tmp_tbl.(EEG_MEASURES{meas_i}))+1*std(tmp_tbl.(EEG_MEASURES{meas_i}))];
                
                x_lim_calc = [min(tmp_tbl.(KINNAMES{kin_i}))-std(tmp_tbl.(KINNAMES{kin_i})),...
                    max(tmp_tbl.(KINNAMES{kin_i}))+std(tmp_tbl.(KINNAMES{kin_i}))];
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
                    % ss.SizeData = 2;
                    ss.SizeData = 25;
                    % ss.MarkerFaceAlpha = 'flat';
                    % ss.AlphaData = repmat(0.3,[size(data,1),1]);
                    ss.MarkerFaceColor = "flat";
                    ss.MarkerFaceAlpha = 0.4;
                    ss.MarkerEdgeAlpha = 0.6; %0.6;
                    if meas_i == length(EEG_MEASURES) && grp_i == length(groups)
                        cond_plot_store = [cond_plot_store, ss];
                    end
                end
                %## LINEAR MODEL FIT
                hold on;
                for cond_i = 1:length(speed_ns)
                    inds = tmp_tbl.speed_n==speed_ns(cond_i);
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
                        xx = sort([max(data.(KINNAMES{kin_i})),min(data.(KINNAMES{kin_i})),x_lim_calc]);
                        % xx = sort([unique(data.(KINNAMES{kin_i}));x_lim_calc']);
                        yy = double(string(speed_ns(cond_i)));
                        zz = zeros(length(xx),1);
                        for i = 1:length(xx)
                            zz(i,1) = xx(i)*slope_kin + yy*slope_sp + yy*xx(i)*slope_intact + inter_mn;
                        end
                        pp = plot(xx,zz,...
                            'DisplayName',sprintf('LME Fit %0.2f (m/s)',speed_ns(cond_i)),...
                            'LineWidth',LINE_WIDTH_MEFF);
                        pp.Color = [COLORS(cond_i,:)*0.8,0.8];
                        % if meas_i == 1
                        %     cond_plot_store = [cond_plot_store, pp];
                        % end
                        %## GIVE CORRELATION COEFFICIENT
                        if DO_PLOT_STAT_STR
                            chkl = find((anova_p_intact <= [0.05,0.01,0.001]),1,'last');
                            if ~isempty(chkl) & cond_i == 1
                                ast = repmat('*',[1,chkl]);
                                % str = sprintf('%s R^2=%0.2g\nm_{speed}=%0.2g\nm_{kin}=%0.2g\nm_{sp:kin}=%0.2g',ast,R2,slope_sp,slope_kin,slope_intact);
                                % str = sprintf('%s R^2=%0.2g  f^2=%0.2g\ny=%0.1g*sp+%0.1g*kin+%0.1g*sp:kin+%0.1g',ast,R2,F2,slope_sp,slope_kin,slope_intact,inter_mn);
                                % str = sprintf('%s R^2=%0.2g  f^2=%0.2g\ny=(%0.1g)sp+(%0.1g)kin+(%0.1g)sp*kin+(%0.1g)',ast,R2,F2,slope_sp,slope_kin,slope_intact,inter_mn);
                                % str = sprintf('%s R^2=%0.2g\nf_{speed}^2=%0.2g  f_{kin}^2=%0.2g',ast,R2,F2_sp,F2_kin);
                                fprintf('CL%i) m_{kin} = %0.2g\n',r_clusters(cl_i),slope_kin);
                                fprintf('CL%i) m_{sp} = %0.2g\n',r_clusters(cl_i),slope_sp);
                                fprintf('CL%i) m_{intact} = %0.2g\n',r_clusters(cl_i),slope_intact);
                                fprintf('CL%i) b_{0} = %0.2g\n\n',r_clusters(cl_i),inter_mn);
                                m_25 = (slope_intact*0.25+slope_kin); 
                                b_25 = slope_sp*0.25+inter_mn;
                                m_1 = (slope_intact*1+slope_kin);
                                b_1 = slope_sp*1+inter_mn;
                                str = sprintf('%s R^2=%0.2g\nf^2=%0.2g\n\ny_{0.25}=(%0.2g)x+(%0.2g)\ny_{1.0}=(%0.2g)x+(%0.2g)',ast,R2,fsq_intact,m_25,b_25,m_1,b_1);
                                
                                r2_pos = [x_shift+REG_X_SHIFT*IM_RESIZE-(0.35/4),y_shift+REG_Y_SHIFT*IM_RESIZE-(0.35/4),0.35,0.35];
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
                if meas_i == length(EEG_MEASURES) && grp_i == length(groups)
                    %- lg2                
                    legend(gca,cond_plot_store);
                    [lg2,icons,plots,txt]  = legend('boxoff');
                    tmp = get(lg2,'String');
                    cnt = 1;
                    for i = 1:length(cond_plot_store)
                        tmp{i} = sprintf('%0.2g m/s',double(string(speed_ns(cnt))));
                        cnt = cnt + 1;
                    end
                    set(lg2,'String',tmp,'FontName',AX_FONT_NAME,'FontSize',LEG_TXT_SIZE)
                    set(lg2,'Orientation','horizontal')
                    set(lg2,'Units','normalized')
                    set(lg2,'Position',[AX_INIT_X+LEG_X_SHIFT*IM_RESIZE,...
                        y_shift+AX_H*IM_RESIZE+LEG_Y_SHIFT*IM_RESIZE,lg2.Position(3),lg2.Position(4)]);
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
                if meas_i == 1 && grp_i == 1
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
        end
        hold off;
        %##
        % exportgraphics(fig,[tmp_savedir filesep sprintf('intacteff_cl%s_kin%s_plot_allsubj.pdf',string(clusters(cl_i)),KINNAMES{kin_i})],'ContentType','vector')
        exportgraphics(fig,[tmp_savedir filesep sprintf('intacteff_cl%s_kin%s_plot_groups.tiff',string(r_clusters(cl_i)),KINNAMES{kin_i})],...
            'Resolution',600)
        % close(fig)
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
for cl_i = 1:length(r_clusters)
    %##
    for kin_i = 1:length(KINNAMES)
        %##
        atlas_name = cluster_titles{double(string(r_clusters(cl_i)))};
        fig = figure('color','white');
        sgtitle(atlas_name,'FontName',AX_FONT_NAME,'FontSize',TITLE_TXT_SIZE,'FontWeight','bold','Interpreter','none');
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        hold on;
        set(gca,AXES_DEFAULT_PROPS{:})
        %-
        x_shift = AX_INIT_X;
        hz = 1;
        y_shift = AX_INIT_Y;
        %##
        for grp_i = 1:length(r_groups)
            x_shift = 0;
            horiz_shift = 0;
            stats_store = [];
            for meas_i = 1:length(EEG_MEASURES)
                %##
                cond_plot_store = [];
                group_plot_store = [];
                if ~isempty(grp_i)
                    inds = KIN_TABLE.model_n == des_i &... %designs(des_i) &...
                        KIN_TABLE.cluster_n == r_clusters(cl_i) &...
                        KIN_TABLE.group_n == r_groups(grp_i);
                else
                    inds = KIN_TABLE.model_n == des_i &... %designs(des_i) &...
                        KIN_TABLE.cluster_n == r_clusters(cl_i);
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
                tmp_stats = RSTATS_IMPORT.cluster_num==double(string(r_clusters(cl_i))) &...
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
                if hz < X_DIM
                    x_shift = x_shift + AX_X_SHIFT*IM_RESIZE*AX_W;
                else
                    y_shift = y_shift - AX_Y_SHIFT*IM_RESIZE*AX_H;
                    x_shift = AX_INIT_X;
                    hz = 0;
                end
                hz = hz + 1;
            end
        end
        hold off;
        %##
        exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_kin%s_plot_3dplot.tiff',string(r_clusters(cl_i)),KINNAMES{kin_i})],'Resolution',300)
        close(fig)
    end
end