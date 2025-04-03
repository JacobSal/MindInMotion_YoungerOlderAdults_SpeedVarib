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
STUDY_DESI_PARAMS = {{'subjselect',{},...
            'variable2','cond','values2',{'flat','low','med','high'},...
            'variable1','group','values1',{'H1000','H2000','H3000'}},...
            {'subjselect',{},...
            'variable2','cond','values2',{'0p25','0p5','0p75','1p0'},...
            'variable1','group','values1',{'H1000','H2000','H3000'}}};
STUDY.cache = [];
for des_i = 1:length(STUDY_DESI_PARAMS)
    [STUDY] = std_makedesign(STUDY,ALLEEG,des_i,STUDY_DESI_PARAMS{des_i}{:});
end
%## LOAD CLUSTER DATA
cl_struct = par_load(cluster_k_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
[comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
CLUSTER_PICS = main_cl_inds;
%% ===================================================================== %%
%-
des_i = 2;
cl_n = 3;
s_chars = {STUDY.datasetinfo(STUDY.cluster(cl_n).sets).subject};
desdes = cat(1,STUDY.design.variable);
% c_chars = desdes(strcmp({desdes.label},'cond'));
% c_chars = {c_chars.value};
g_chars = desdes(strcmp({desdes.label},'group'));
g_chars = {g_chars.value};
g_chars = g_chars{1};
%--
cluster_colors = {[1 1 1],...        % White
            [1 1 0]...             % Yellow
            [221,52,151]/255,...    % Pink
            [1 0 0],...             % Red
            [250 140 0]/255,...     % oragne
            [210 173 255]/255,...   % purple
            [0.5 0.5 0],...         % Olive
            [0.5 0 0.5],...         % Purple
            [0.5 0 0],...           % Maroon
            [0 1 1],...             % Aqua
            [0 1 0],...             % Lime
            [0 0.5 0.5],...         % Teal
            [0 0.5 0],...           % Green
            [0 0 1],...             % Blue
            [0 0 0.5],...           % Navy
            [0.8 0.8 0.8]};          % Gray
cluster_colors = cluster_colors([4 11 14 2 13 10 5 6 15 16 1 7 9 3]);
%-- rb3
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
%##
%--
% dip_dir = [cluster_k_dir filesep 'topo_dip_inf'];
dip_dir = [cluster_k_dir filesep 'topo_dip_inf' filesep 'valid_clusts'];
%--
cluster_inds_plot = [3,4,5,6,7,8,9,11,12];
clusters = main_cl_inds;
%% (SUBJECT DIPOLES PLOTS) ============================================= %%
%-
AX_INIT_HORIZ_DIP = 0.075;
AX_INIT_VERT_DIP = 0.6;
AX_INIT_VERT_DIP2 = 0.775;
%## FIGURE IMPLEMENTATION
%-
% AX1_RESIZE = 0.25;
% AX2_RESIZE = 0.39;
% AX3_RESIZE = 0.302;
IM_RESIZE = 1.5;
AX2_ADJ = 1.52;
AX3_ADJ = 1.255;
%-
AX2_X_O = 0.19;
AX2_Y_O = -0.17;
AX3_X_O = 0.24;
AX3_Y_O = -0.1;
%-
FONT_NAME = 'Arial';
%-
AXES_DEFAULT_PROPS = {'box','off', ...
    'xtick',[], ...
    'ytick',[], ...
    'ztick',[], ...
    'xcolor',[1,1,1], ...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
    'ycolor',[1,1,1]};
%% PLOT ============================================================= %%
fig = figure('color','white', ...
    'Renderer','painters');
set(fig,'Units','inches', ...
    'Position',[0.5,0.5,6,9.5], ...
    'PaperUnits','inches', ...
    'PaperSize',[1 1], ...
    'PaperPosition',[0 0 1 1])
p_sz = get(fig,'Position');
set(gca,AXES_DEFAULT_PROPS{:})
hold on;

%## NEXT ROW =========================================================== %%
tmp = openfig([dip_dir filesep 'dipplot_alldipspc_top.fig']);
cnt = length(tmp.Children(end).Children);
while cnt > 3
    tmp.Children(end).Children(cnt).LineWidth = 0.02; %'none';
    tmp.Children(end).Children(cnt).CData = [0,0,0]; %tmp.Children(end).Children(cnt-1).MarkerFaceColor;
    tmp.Children(end).Children(cnt).MarkerFaceAlpha = 0.85;
    tmp.Children(end).Children(cnt).SizeData = 6; %fig.Children(end).Children(cnt+3+1).MarkerSize*2.6;
    cnt = cnt - 1;
end
ax = get(tmp,'CurrentAxes');
view(ax,[90,0])
daspect = get(ax,'DataAspectRatio');
ac = get(ax,'Children');
XLIM = [ac(1).XData(1,1),ac(1).XData(1,2)];
YLIM = [ac(1).YData(1,1),ac(1).YData(2,1)];
ZLIM = [ac(2).ZData(1,1),ac(2).ZData(2,1)];
dx = abs(ac(1).XData(1,2) - ac(1).XData(1,1));
dy = abs(ac(1).YData(1,1) - ac(1).YData(2,1));
dz = abs(ac(2).ZData(1,1) - ac(2).ZData(2,1));

%## AXES 1
ax1 = axes('Parent',fig,'DataAspectRatio',daspect,'Units','normalized');
copyobj(ac,ax1);
set(ax1,AXES_DEFAULT_PROPS{:});
ax1.XRuler.FirstCrossoverValue = 0;
ax1.YRuler.FirstCrossoverValue = 0;
ax1.ZRuler.FirstCrossoverValue = 0;
view([0,90])
camzoom(ax1,1.1^2)
set(ax1,'YLim',YLIM)
set(ax1,'XLim',XLIM)
drawnow;
pp = get(ax1,'Position');
ax1_x = pp(3)*IM_RESIZE/p_sz(3);
ax1_y = pp(4)*IM_RESIZE/p_sz(4);
set(ax1,'OuterPosition',[0,0,1,1],'Position',[AX_INIT_HORIZ_DIP,AX_INIT_VERT_DIP,ax1_x,ax1_y]);

%## AXES 2
ax2 = axes('Parent',fig,'DataAspectRatio',daspect,'Units','normalized');
copyobj(ac,ax2);
%-
set(ax2,AXES_DEFAULT_PROPS{:});
ax2.XRuler.FirstCrossoverValue = 0;
ax2.YRuler.FirstCrossoverValue = 0;
ax2.ZRuler.FirstCrossoverValue = 0;
set(ax2,'View',[90,0]);
camzoom(ax2,1.1^2)
set(ax2,'YLim',YLIM)
set(ax2,'ZLim',ZLIM)
drawnow;
pp2 = get(ax2,'Position');
ax2_x = pp2(3)*IM_RESIZE*AX2_ADJ/p_sz(3);
ax2_y = pp2(4)*IM_RESIZE*AX2_ADJ/p_sz(4);
set(ax2,'OuterPosition',[0,0,1,1],'Position',[AX_INIT_HORIZ_DIP+ax1_x+(ax2_x*AX2_X_O),AX_INIT_VERT_DIP+(ax2_y*AX2_Y_O),ax2_x,ax2_y]);

%## AXES 3
ax3 = axes('Parent',fig,'DataAspectRatio',daspect,'Units','normalized');
copyobj(ac,ax3);
%-
set(ax3,AXES_DEFAULT_PROPS{:});
ax3.XRuler.FirstCrossoverValue = 0;
ax3.YRuler.FirstCrossoverValue = 0;
ax3.ZRuler.FirstCrossoverValue = 0;
set(ax3,'view',[0,0])
set(ax3,'ZLim',ZLIM)
set(ax3,'XLim',XLIM)
camzoom(ax3,1.1^2)
pp3 = get(ax3,'Position');

ax3_x = pp3(3)*IM_RESIZE*AX3_ADJ/p_sz(3);
ax3_y = pp3(4)*IM_RESIZE*AX3_ADJ/p_sz(4);
% AX3_W_O = 0.24;
% AX3_H_O = -0.1;
% ax3.SortMethod='ChildOrder';
set(ax3,'OuterPosition',[0,0,1,1],'Position',[AX_INIT_HORIZ_DIP+ax1_x+ax2_x+(ax2_x*AX2_X_O)+(ax3_x*AX3_X_O),AX_INIT_VERT_DIP+(ax3_y*AX3_Y_O),ax3_x,ax3_y]);
close(tmp)

%## FIRST ROW ========================================================== %%
tmp = openfig([dip_dir filesep 'dipplot_avgdipspc_top.fig']);
cnt = length(tmp.Children(end).Children);
surf = tmp.Children(end).Children(1);
while cnt > 3
    hold on;
    tmp.Children(end).Children(cnt).Marker = 'o'; 
    tmp.Children(end).Children(cnt).LineWidth = 0.05; %'none';
    tmp.Children(end).Children(cnt).MarkerFaceAlpha = 0.85;
    tmp.Children(end).Children(cnt).SizeData = 100; %fig.Children(end).Children(cnt+3+1).MarkerSize*2.6;
    cnt = cnt - 1;
end
ax = get(tmp,'CurrentAxes');
daspect = get(ax,'DataAspectRatio');
ac = get(ax,'Children');
XLIM = [ac(1).XData(1,1),ac(1).XData(1,2)];
YLIM = [ac(1).YData(1,1),ac(1).YData(2,1)];
ZLIM = [ac(2).ZData(1,1),ac(2).ZData(2,1)];
dx = abs(ac(1).XData(1,2) - ac(1).XData(1,1));
dy = abs(ac(1).YData(1,1) - ac(1).YData(2,1));
dz = abs(ac(2).ZData(1,1) - ac(2).ZData(2,1));
%-
ax1 = axes('Parent',fig,'DataAspectRatio',daspect,'Units','normalized');
copyobj(ac,ax1);
%-
set(ax1,AXES_DEFAULT_PROPS{:});
ax1.XRuler.FirstCrossoverValue = 0;
ax1.YRuler.FirstCrossoverValue = 0;
ax1.ZRuler.FirstCrossoverValue = 0;
view([0,90])
camzoom(ax1,1.1^2)
set(ax1,'YLim',YLIM)
set(ax1,'XLim',XLIM)
drawnow;
pp = get(ax1,'Position');
ax1_x = pp(3)*IM_RESIZE/p_sz(3);
ax1_y = pp(4)*IM_RESIZE/p_sz(4);
set(ax1,'OuterPosition',[0,0,1,1],'Position',[AX_INIT_HORIZ_DIP,AX_INIT_VERT_DIP2,ax1_x,ax1_y]);

%## AXES 2
ax2 = axes('Parent',fig,'DataAspectRatio',daspect,'Units','normalized');
copyobj(ac,ax2);
%-
set(ax2,);
ax2.XRuler.FirstCrossoverValue = 0;
ax2.YRuler.FirstCrossoverValue = 0;
ax2.ZRuler.FirstCrossoverValue = 0;
set(ax2,'View',[90,0]);
camzoom(ax2,1.1^2)
set(ax2,'YLim',YLIM)
set(ax2,'ZLim',ZLIM)
drawnow;
pp2 = get(ax2,'Position');
ax2_x = pp2(3)*IM_RESIZE*AX2_ADJ/p_sz(3); 
ax2_y = pp2(4)*IM_RESIZE*AX2_ADJ/p_sz(4);
set(ax2,'OuterPosition',[0,0,1,1],'Position',[AX_INIT_HORIZ_DIP+ax1_x+(ax2_x*AX2_X_O),AX_INIT_VERT_DIP2+(ax2_y*AX2_Y_O),ax2_x,ax2_y]);

%## AXES 3
ax3 = axes('Parent',fig,'DataAspectRatio',daspect,'Units','normalized');
copyobj(ac,ax3);
%-
set(ax3,AXES_DEFAULT_PROPS{:});
ax3.XRuler.FirstCrossoverValue = 0;
ax3.YRuler.FirstCrossoverValue = 0;
ax3.ZRuler.FirstCrossoverValue = 0;
set(ax3,'view',[0,0])
set(ax3,'ZLim',ZLIM)
set(ax3,'XLim',XLIM)
camzoom(ax3,1.1^2)
pp3 = get(ax3,'Position');
ax3_x = pp3(3)*IM_RESIZE*AX3_ADJ/p_sz(3);
ax3_y = pp3(4)*IM_RESIZE*AX3_ADJ/p_sz(4);
set(ax3,'OuterPosition',[0,0,1,1],'Position',[AX_INIT_HORIZ_DIP+ax1_x+ax2_x+(ax2_x*AX2_X_O)+(ax3_x*AX3_X_O),AX_INIT_VERT_DIP2+(ax3_y*AX3_Y_O),ax3_x,ax3_y]);
close(tmp)

%## ANNOTATIONS ======================================================== %%
%-
ANN_H = 0.875;
annotation('textbox',[0.15,ANN_H,.1,.1],...
    'String','Transverse','HorizontalAlignment','center',...
    'VerticalAlignment','top','LineStyle','none','FontName',FONT_NAME,...
    'FontSize',12,'FontWeight','Bold','Units','normalized');
%-
annotation('textbox',[0.43,ANN_H,.1,.1],...
    'String','Sagittal','HorizontalAlignment','center',...
    'VerticalAlignment','top','LineStyle','none','FontName',FONT_NAME,...
    'FontSize',12,'FontWeight','Bold','Units','normalized');
%-
annotation('textbox',[0.745,ANN_H,.1,.1],...
    'String','Coronal','HorizontalAlignment','center',...
    'VerticalAlignment','top','LineStyle','none','FontName',FONT_NAME,...
    'FontSize',12,'FontWeight','Bold','Units','normalized');
hold on;
bc_shift_x = -0.035;
bc_shift_y = -0.07;
p_sz = get(fig,'Position');
circle_dim = 0.2;
circle_dim_x = circle_dim/(p_sz(3));
circle_dim_y = circle_dim/(p_sz(4));
txt_dim_y = 1/p_sz(3);
txt_dim_x = 2.25/p_sz(4);
txt_shift_x = 0.35/p_sz(4);
txt_shift_y = -0.43/p_sz(3);
HZ_DIM = 3;
AX_INIT_HORIZ = AX_INIT_HORIZ_DIP;
HORIZ_SHIFT_VAL = 13.5/p_sz(4);
VERT_SHIFT_VAL = -2/p_sz(3);
hz = 0;
x_shift = AX_INIT_HORIZ;
y_shift = AX_INIT_VERT_DIP;
aa_bc = [];
aa_txt = [];
for cl_i = 1:length(cluster_inds_plot)
    %-- initiate params
    cl_ii = find(cluster_inds_plot(cl_i) == double(string(clusters)));
    cl_n = double(string(clusters(cl_ii)));
    atlas_name = cluster_titles{cl_ii};
    %--
    hz = hz + 1;
    ANN_H = 0.875;
    tmp = annotation(fig,'ellipse',[bc_shift_x+x_shift,bc_shift_y+y_shift,circle_dim_x,circle_dim_y],...
        'Units','normalized',...
        'FaceColor',cluster_colors{cluster_inds_plot(cl_i)});
    aa_bc = [aa_bc, tmp];
    tmp = annotation('textbox',[bc_shift_x+x_shift+txt_shift_x,bc_shift_y+y_shift+txt_shift_y,txt_dim_x,txt_dim_y],...
        'String',sprintf('%s',atlas_name), ...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle', ...
        'LineStyle','none', ...
        'FontName',FONT_NAME,...
        'FontSize',12, ...
        'FontWeight','Bold', ...
        'Units','normalized');
    aa_txt = [aa_txt, tmp];
    if hz < HZ_DIM
        x_shift = x_shift + HORIZ_SHIFT_VAL*txt_dim_x;
    else
        y_shift = y_shift + VERT_SHIFT_VAL*txt_dim_y;
        x_shift = AX_INIT_HORIZ;
        hz = 0;
    end    
end
%-- edits (03/03/2025)
aa_bc(end).Position(1) = aa_bc(end).Position(1) + HORIZ_SHIFT_VAL*txt_dim_x;
aa_txt(end).Position(1) = aa_txt(end).Position(1) + HORIZ_SHIFT_VAL*txt_dim_x;
hold off;
%## SAVE
exportgraphics(fig,[save_dir filesep 'figure_group_dipoleplots.pdf'],'ContentType','vector')
exportgraphics(fig,[save_dir filesep 'figure_group_dipoleplots.tiff'],'ContentType','image','Resolution',600)

%% (BLANK MRI SLICES) ================================================== %%
HIRES_FNAME = 'mni_icbm152_t1_tal_nlin_sym_09a.nii';
HIRES_FPATH = [PATHS.data_dir filesep '_resources' filesep 'mni_icbm152_nlin_sym_09a'];
%- assign hires_template default
fname = strsplit(HIRES_FNAME,'.');
hires_mri = [HIRES_FPATH filesep fname{1} '_dipplotmri.mat'];
mri = load(hires_mri);
mri = mri.mri;

%## PLOIT
% FIGURE_POSITION = [1,1,6.5,9];
% fig = figure('color','white', ...
%         'Renderer','Painters');
% set(fig,'Units','inches', ...
%         'Position',FIGURE_POSITION, ...
%         'PaperUnits','inches', ...
%         'PaperSize',[1 1], ...
%         'PaperPosition',[0 0 1 1])
%-- fieldtrip check
% cfg = struct('method','ortho', ...
%     'anaparamter','anatomy', ...
%     'location','center', ...
%     'crosshair','no', ...
%     'axis','off');
% ft_sourceplot(cfg,mri)

%## GET MRI DIMS
SCALE_FACTOR = 0.9;
CROP_THRESH = 0;
%--
sz = size(mri.anatomy);
sz1 = floor(sz(1)/2);
sz2 = floor(sz(2)/2);
sz3 = floor(sz(3)/2);
axs1 = sz(1)/max(sz)*SCALE_FACTOR;
axs2 = sz(2)/max(sz)*SCALE_FACTOR;
axs3 = sz(3)/max(sz)*SCALE_FACTOR;
%-- sagittal plane
tmpc = (squeeze(mri.anatomy(sz1,:,sz(3):-1:1)));
tmpc(tmpc<CROP_THRESH) = 255;
fig = figure('color','white', ...
        'Renderer','Painters');
ax = axes();
% ss = surf(ax,tmpc, ...
%     'FaceColor','interp', ...
%     'EdgeColor','texturemap', ...
%     'FaceLighting','flat', ...
%     'EdgeLighting','flat');
% ss = surf(ax,tmpc, ...
%     'FaceColor','texturemap', ...
%     'EdgeColor','interp', ...
%     'FaceLighting','flat', ...
%     'EdgeLighting','flat');
ss = surf(ax,tmpc, ...
    'FaceColor','flat', ...
    'EdgeColor','texturemap', ...
    'FaceLighting','none', ...
    'EdgeLighting','none');
%(03/21/2025) JS, seems to be most true to the MRI
set(ax,'box','off', ...
    'xtick',[], ...
    'ytick',[],...
    'ztick',[], ...
    'xcolor',[1,1,1], ...
    'ycolor',[1,1,1], ...
    'zcolor',[1,1,1], ...
    'Position',[0.05,0.05,axs2,axs3]);
colormap('gray')
view([90,90])
% exportgraphics(fig,[save_dir filesep 'mri_sagittal_slice.pdf'],'ContentType','Vector'); 
%(03/21/2025) vector save takes 5-ever
exportgraphics(fig,[save_dir filesep 'mri_sagittal_slice.tif'],'Resolution',1000);
%-- coronal plane
tmpc = (squeeze(mri.anatomy(:,sz2,sz(3):-1:1)));
fig = figure('color','white', ...
        'Renderer','Painters');
ax = axes();
% ss = surf(ax,tmpc, ...
%     'FaceColor','interp', ...
%     'EdgeColor','interp', ...
%     'FaceLighting','gouraud');
ss = surf(ax,tmpc, ...
    'FaceColor','flat', ...
    'EdgeColor','texturemap', ...
    'FaceLighting','flat', ...
    'EdgeLighting','flat');
set(ax,'box','off', ...
    'xtick',[], ...
    'ytick',[],...
    'ztick',[], ...
    'xcolor',[1,1,1], ...
    'ycolor',[1,1,1], ...
    'zcolor',[1,1,1], ...
    'Position',[0.05,0.05,axs1,axs3]);
colormap('gray')
view([90,90])
% exportgraphics(fig,[save_dir filesep 'mri_coronal_slice.pdf'],'ContentType','Vector');
exportgraphics(fig,[save_dir filesep 'mri_coronal_slice.tif'],'Resolution',1000);
%-- top plane
tmpc = (squeeze(mri.anatomy(:,:,sz3)));
fig = figure('color','white', ...
        'Renderer','Painters');
ax = axes();
% ss = surf(ax,tmpc, ...
%     'FaceColor','interp', ...
%     'EdgeColor','interp', ...
%     'FaceLighting','gouraud');
ss = surf(ax,tmpc, ...
    'FaceColor','flat', ...
    'EdgeColor','texturemap', ...
    'FaceLighting','flat', ...
    'EdgeLighting','flat');
set(ax,'box','off', ...
    'xtick',[], ...
    'ytick',[],...
    'ztick',[], ...
    'xcolor',[1,1,1], ...
    'ycolor',[1,1,1], ...
    'zcolor',[1,1,1], ...
    'Position',[0.05,0.05,axs1,axs2]);
colormap('gray')
view([-90,90])
% exportgraphics(fig,[save_dir filesep 'mri_top_slice.pdf'],'ContentType','Vector');
exportgraphics(fig,[save_dir filesep 'mri_top_slice.tif'],'Resolution',1000);