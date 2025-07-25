%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/run_f_cluster_info_plots.sh

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
        STUDY_DIR = SCRIPT_DIR; % change this if in sub folder
        SRC_DIR = STUDY_DIR;
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
    STUDY_DIR = SCRIPT_DIR; % change this if in sub folder
    SRC_DIR = STUDY_DIR;
end
%## Add Study, Src, & Script Paths
addpath(SRC_DIR);
addpath(STUDY_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca_speed');
%% (PATHS) ============================================================= %%
%## DATASET
DATA_SET = 'MIM_dataset';
%## STUDY INFO
% STUDY_DNAME = '10172024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej_speed';
STUDY_DNAME = '02202025_mim_yaoa_powpow0p3_crit_speed';
% STUDY_FNAME = 'epoch_study';
%## soft define
studies_dir = [PATHS.data_dir filesep DATA_SET filesep '_studies'];
%## CLUSTER
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
% cluster_fpath = [studies_dir filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb3'];
cluster_fpath = [studies_dir filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_allcond_rb3'];

cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
cluster_k_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
%- save dir
save_dir = [cluster_k_dir filesep 'topo_dip_inf' filesep 'extra_clusts'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% =============================================label_i===================== %%
%## LOAD STUDY
if ~ispc
    [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '_UNIX.study'],'filepath',cluster_study_fpath);
else
    [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '.study'],'filepath',cluster_study_fpath);
end
%## LOAD STUDY
% if ~ispc
%     tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_NAME)]);
%     STUDY = tmp.STUDY;
% else
%     tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s.study',CLUSTER_STUDY_NAME)]);
%     STUDY = tmp.STUDY;
% end
cl_struct = par_load([cluster_study_fpath filesep sprintf('%i',CLUSTER_K)],sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
[comps_out,main_cl_inds,outlier_cl_inds,valid_clusters] = eeglab_get_cluster_comps(STUDY);
% cluster_inds = main_cl_inds(1:end);
% cluster_inds = [3,4,5,6,7,8,9,11,12];
cluster_inds = [10,13];
% cluster_inds = valid_clusters;
% cluster_inds = [3,4,5,7,10,12,13]; % 01192025_mim_yaoa_nopowpow_crit_speed, rb3

% cluster_inds = [3,4,6,7,8,9,10,13]; % 01192025_mim_yaoa_nopowpow_crit_speed, rb3
%- save dir
% save_dir = [save_dir filesep 'figure_gen']; %sprintf('cl%i-cl%i',min(cluster_inds),max(cluster_inds))];
% if ~exist(save_dir,'dir')
%     mkdir(save_dir);
% end
%% (ANATOMY) =========================================================== %%
%{
%## AAL3 PATH
if ~ispc
    tp = strsplit(path,':');
else
    tp = strsplit(path,';');
end
b1 = contains(tp,'AAL3','IgnoreCase',true);
b2 = tp(b1);
try
    ind = regexp(b2(1),'AAL3','end');
    path_aal3 = b2{1}(1:ind{1});
    fprintf('ALL3 path: %s\n',path_aal3);
catch ME
    switch ME.identifier
        case 'MATLAB:badsubscript'
            fprintf('AAL3 path not found.\n');
    end
end
addpath(path_aal3);
%## ANATOMY CALCULATION
% [path_aal3 filesep 'ROI_MNI_V7_1mm.nii']
ANATOMY_STRUCT = struct('atlas_fpath',{{[path_aal3 filesep 'AAL3v1_1mm.nii.gz'], ...
        [path_aal3 filesep 'ROI_MNI_V7_1mm.nii']}},...
    'group_chars',{unique({STUDY.datasetinfo.group})},...
    'cluster_inds',3:length(STUDY.cluster),...
    'anatomy_calcs',{{'all aggregate','all centroid'}},... % ('all calcs','group centroid','all centroid','group aggregate','all aggregate')
    'save_dir',save_dir,...
    'save_inf',true);
%(01/16/2025) JS, 'all centroid' option has a bug where indexing
%doesn't seem to want to work on line 535 in poi_box
%(03/05/2025) JS, this seems to be fixed after some exception
%handling changes
[STUDY,anat_struct,~,~,txt_out] = eeglab_get_anatomy(STUDY, ...
    'ALLEEG',ALLEEG, ...
    'ANATOMY_STRUCT',ANATOMY_STRUCT);
%}
%% TOPO ================================================================ %%
groups = unique({STUDY.datasetinfo.group});
AX_HORIZ_SHIFT = 0.4;
AX_VERT_SHIFT = 0.05;
DIP_IM_DPI = 900;
AX_INIT_HORIZ_TOPO = 0.085;
AX_INIT_VERT_TOPO = 0.75;
FIGURE_POSITION =[1,1,6.5,9];
PG_SIZE = [6.5,9];
FONT_NAME = 'Arial';
TOPO_FONTSIZE = 6;
HZ_DIM = 3;
g_chars = unique({STUDY.datasetinfo.group});
g_chars_topo = {'Young Adults','Older High Mobility Adults','Older Low Mobility Adults'};
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],...
        'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
% AX_H = 0.25;
% AX_W = 0.225;
% IM_RESIZE = 0.7;
IM_RESIZE = 0.225;
horiz_shift = 0;
vert_shift = 0;
hz = 1;
%## TOPO PLOT
fig = figure('color','white');
set(fig,'Units','inches','Position',FIGURE_POSITION)
set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
set(gca,AXES_DEFAULT_PROPS{:});
hold on;
for k_i = 1:length(cluster_inds)
    cl_i = double(string(cluster_inds(k_i)));
    %## TOPO PLOT
    [~,h] = std_topoplot_CL(STUDY,cl_i,'together');
    set(h,'color','w')
    set(h,'Units','normalized');    
    ax = get(h,'CurrentAxes');
    colormap(h,linspecer);   
    pos = get(ax,'Position');
    opos = get(ax,'OuterPosition');
    ac = get(ax,'Children');
    ax1 = axes('Parent',fig,...
        'DataAspectRatio',get(ax,'DataAspectRatio'));
    copyobj(ac,ax1);
    set(ax1,AXES_DEFAULT_PROPS{:});
    colormap(ax1,linspecer)
    g_counts = cell(length(groups),1);
    for g_i = 1:length(groups)
        g_inds = cellfun(@(x) strcmp(x,g_chars{g_i}),{STUDY.datasetinfo(STUDY.cluster(cl_i).sets).group});
        if length(g_chars_topo{g_i}) == 1 || ischar(g_chars_topo{g_i})
            g_counts{g_i} =sprintf('%s N=%i',g_chars_topo{g_i},sum(g_inds));
        else
            g_counts{g_i} =sprintf('%s\n%s N=%i',g_chars_topo{g_i}{1},g_chars_topo{g_i}{2},sum(g_inds));
        end
    end
    ax1.Title.String = g_counts; %sprintf('N=%i',length(STUDY.cluster(cl_i).sets));
    ax1.Title.Interpreter = 'none';
    ax1.FontSize = TOPO_FONTSIZE; %PLOT_STRUCT.font_size;
    ax1.FontName = FONT_NAME;
    ax1.FontWeight = 'bold';
    ax1.OuterPosition = [0,0,1,1];
    ax1.Units = 'Normalized';
    pp = get(ax1,'Position');    
    ax1.Position = [AX_INIT_HORIZ_TOPO+horiz_shift,AX_INIT_VERT_TOPO+vert_shift,pp(3)*IM_RESIZE,pp(4)*IM_RESIZE];  %[left bottom width height]
    close(h);
    % std_topoplot_CL(STUDY,cl_i,'together');
    % colormap(linspecer); %colormap_ersp)
    % fig_i = get(groot,'CurrentFigure');
    % set(fig_i,'color','w')
    % g_counts = cell(length(g_chars),1);
    % for group_i = 1:length(g_chars)
    %     g_inds = cellfun(@(x) strcmp(x,g_chars{group_i}),{STUDY.datasetinfo(STUDY.cluster(cl_i).sets).group});
    %     if length(g_chars_topo{group_i}) == 1 || ischar(g_chars_topo{group_i})
    %         g_counts{group_i} =sprintf('%s N=%i',g_chars_topo{group_i},sum(g_inds));
    %     else
    %         g_counts{group_i} =sprintf('%s\n%s N=%i',g_chars_topo{group_i}{1},g_chars_topo{group_i}{2},sum(g_inds));
    %     end
    % end
    % obj_i = findobj(fig_i,'type','Axes');
    % obj_i(1).Title.String = [{sprintf('CL%i',cl_i)};g_counts]; %sprintf('N=%i',length(STUDY.cluster(cl_i).sets));
    % obj_i(1).Title.Interpreter = 'none';
    % obj_i(1).FontSize = TOPO_FONTSIZE; %PLOT_STRUCT.font_size;
    % obj_i(1).FontName = FONT_NAME;
    % obj_i(1).FontWeight = 'bold';
    % obj_i(1).OuterPosition = [0,0,1,1];
    % obj_i(1).Units = 'Normalized';
    % obj_i(1).Position = [AX_INIT_HORIZ_TOPO+horiz_shift,AX_INIT_VERT_TOPO+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE];  %[left bottom width height]
    %##
    if hz < HZ_DIM
        horiz_shift = horiz_shift + pp(3)*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
    else
        vert_shift = vert_shift - (pp(4)*IM_RESIZE + AX_VERT_SHIFT*IM_RESIZE);
        horiz_shift = 0;
        hz = 0;
    end
    hz = hz + 1;
end
hold off;
drawnow;
exportgraphics(fig,[save_dir filesep 'topo_plots.pdf'],'ContentType','vector')
exportgraphics(fig,[save_dir filesep 'topo_plots.tif'],'ContentType','image','Resolution',600)
savefig(fig,[save_dir filesep 'topo_plots.fig'])
% close(fig);

%% NEW DIPOLE IMPLEMENTATION
%--
HIRES_FNAME = 'mni_icbm152_t1_tal_nlin_sym_09a.nii';
HIRES_FPATH = [PATHS.data_dir filesep '_resources' filesep 'mni_icbm152_nlin_sym_09a'];
%- assign hires_template default
fname = strsplit(HIRES_FNAME,'.');
hires_mesh = [HIRES_FPATH filesep fname{1} '_dipplotvol.mat'];
hires_mri = [HIRES_FPATH filesep fname{1} '_dipplotmri.mat'];
mri = load(hires_mri);
mri = mri.mri;
vol = hires_mesh;
%## EEGLAB DEFAULTS
%{
tmp = strsplit(path,';');
% tmp = strsplit(path,';');
b1 = regexp(tmp,'eeglab','end');
b2 = tmp(~cellfun(@isempty,b1));
PATH_EEGLAB = b2{1}; %(1:b1{1});
fprintf('EEGLAB path: %s\n',PATH_EEGLAB);
PATH_EEGLAB_BEM  = [PATH_EEGLAB filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM' filesep];
MNI_VOL = [PATH_EEGLAB_BEM filesep 'standard_vol.mat'];
MNI_MRI = [PATH_EEGLAB_BEM filesep 'standard_mri.mat'];
mri = MNI_MRI;
mri = load(mri);
mri = mri.mri;
vol = MNI_VOL;
%}
%-
% transform = [1,0,0,-99;...
%                  0,1,0,-135;...
%                  0,0,1,-73;...
%                  0,0,0,1];
% transform = [1,0,0,1;...
%                  0,1,0,1;...
%                  0,0,1,1;...
%                  0,0,0,1];
DIPPLOT_STRUCT = struct('rvrange',[0,30],... % this is a value from 0 to 100 (e.g., rv = 0.15 is 15)
        'summary','off',...
        'mri',mri,...
        'coordformat','MNI',...
        'transform',[],...
        'image','mri',...
        'plot','on',...
        'color',{{[0,0,0]}},...
        'view',[.5,.5,.5],...
        'mesh','off',...
        'meshdata',vol,...
        'axistight','off',... % display the closest MRI slice to distribution
        'gui','on',...
        'num','off',...
        'cornermri','on',...
        'drawedges','off',...
        'projimg','off',...
        'projlines','off',...
        'projwidth',1,...
        'projcol',{{[0,0,0]}},...
        'projalpha',0.1, ...
        'dipolesize',1,...
        'dipolelength',0,...
        'pointout','off',...
        'sphere',1,...
        'spheres','off',...
        'normlen','off',...
        'dipnames',{{}},...
        'holdon','on',...
        'camera','auto',...
        'density','off');

%% ALL DIPOLES FOR CLUSTERS ============================================ %%
DIPPLOT_STRUCT.dipolesize = 1.25;
DIPPLOT_STRUCT.axistight = 'off';
[fig] = eeglab_dipplot(STUDY,ALLEEG,cluster_inds,...
    'PLOT_TYPE','all_nogroup',...
    'DIPPLOT_STRUCT',DIPPLOT_STRUCT);
% %-- include the brain mesh volume if desired
% hold on;
% vol = load(DIPPLOT_STRUCT.meshdata);
% mesh = vol.vol.ft_vol;
% tmp = mesh.bnd(1);
% tcmap = linspecer(40);
% ft_plot_mesh(tmp, ...
%     'facecolor',tcmap(30), ...
%     'edgecolor','none', ...
%     'facealpha',0.3, ...
%     'edgealpha',0.3);
% drawnow;
%--
pause(2);
camzoom(1.1^2);
%## TOP
view([0,90])
%-- save
exportgraphics(fig,[save_dir filesep sprintf('dipplot_alldipspc_top.tiff')], ...
    'Resolution',DIP_IM_DPI);
exportgraphics(fig,[save_dir filesep sprintf('dipplot_alldipspc_top.pdf')], ...
    'ContentType','vector', ...
    'Resolution',DIP_IM_DPI);
saveas(fig,[save_dir filesep sprintf('dipplot_alldipspc_top.fig')]);
%## SAGGITAL
view([90,0])
%-- save
exportgraphics(fig,[save_dir filesep sprintf('dipplot_alldipspc_sagittal.tiff')], ...
    'Resolution',DIP_IM_DPI);
exportgraphics(fig,[save_dir filesep sprintf('dipplot_alldipspc_sagittal.pdf')], ...
    'ContentType','vector', ...
    'Resolution',DIP_IM_DPI);
%## CORONAL
view([0,0])
%-- save
exportgraphics(fig,[save_dir filesep sprintf('dipplot_alldipspc_coronal.tiff')], ...
    'Resolution',DIP_IM_DPI);
exportgraphics(fig,[save_dir filesep sprintf('dipplot_alldipspc_coronal.pdf')], ...
    'ContentType','vector', ...
    'Resolution',DIP_IM_DPI);
pause(2);
close(fig);

%% AVERAGE DIPOLE FOR CLUSTERS ========================================= %%
DIPPLOT_STRUCT.dipolesize = 1.25;
DIPPLOT_STRUCT.axistight = 'off';
[fig] = eeglab_dipplot(STUDY,ALLEEG,cluster_inds,...
    'PLOT_TYPE','average_nogroup',...
    'DIPPLOT_STRUCT',DIPPLOT_STRUCT);
% %-- include the brain mesh volume if desired
% hold on;
% vol = load(DIPPLOT_STRUCT.meshdata);
% mesh = vol.vol.ft_vol;
% tmp = mesh.bnd(1);
% tcmap = linspecer(40);
% ft_plot_mesh(tmp, ...
%     'facecolor',tcmap(20,:), ...
%     'edgecolor','none', ...
%     'facealpha',0.3, ...
%     'edgealpha',0.3);
% drawnow;
%--
pause(2);
camzoom(1.1^2);
%## TOP
view([0,90])
%-- save
exportgraphics(fig,[save_dir filesep sprintf('dipplot_avgdipspc_top.tiff')], ...
    'Resolution',DIP_IM_DPI);
exportgraphics(fig,[save_dir filesep sprintf('dipplot_avgdipspc_top.pdf')], ...
    'ContentType','vector', ...
    'Resolution',DIP_IM_DPI);
saveas(fig,[save_dir filesep sprintf('dipplot_avgdipspc_top.fig')]);
%## SAGGITAL
view([90,0])
%-- save
exportgraphics(fig,[save_dir filesep sprintf('dipplot_avgdipspc_sagittal.tiff')], ...
    'Resolution',DIP_IM_DPI);
exportgraphics(fig,[save_dir filesep sprintf('dipplot_avgdipspc_sagittal.pdf')], ...
    'ContentType','vector', ...
    'Resolution',DIP_IM_DPI);
%## CORONAL
view([0,0])
%-- save
exportgraphics(fig,[save_dir filesep sprintf('dipplot_avgdipspc_coronal.tiff')], ...
    'Resolution',DIP_IM_DPI);
exportgraphics(fig,[save_dir filesep sprintf('dipplot_avgdipspc_coronal.pdf')], ...
    'ContentType','vector', ...
    'Resolution',DIP_IM_DPI);
pause(2);
close(fig);

%% INDIVIDUAL CLUSTERS ================================================= %%
DIPPLOT_STRUCT.dipolesize = 1.25;
DIPPLOT_STRUCT.axistight = 'on';
DIP_IM_DPI = 300;
%--
for i = 1:length(cluster_inds)
    cluster_i = cluster_inds(i);
    %## ALL ========================================================== %%
    %{
    [fig] = eeglab_dipplot(STUDY,ALLEEG,cluster_i,...
        'PLOT_TYPE','all_nogroup',...
        'DIPPLOT_STRUCT',DIPPLOT_STRUCT);
    camzoom(1.1^2);
    %## TOP
    view([0,90])
    %-- save
    % exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_nogroup_top.tiff',cluster_i)], ...
    %     'Resolution',DIP_IM_DPI);
    % exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_nogroup_top.pdf',cluster_i)], ...
    %     'ContentType','vector', ...
    %     'Resolution',DIP_IM_DPI);
    saveas(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_nogroup_top.fig',cluster_i)]);
    %## SAGGITAL
    % view([90,0])
    %-- save
    % exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_nogroup_sagittal.tiff',cluster_i)], ...
    %     'Resolution',DIP_IM_DPI);
    % exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_nogroup_sagittal.pdf',cluster_i)], ...
    %     'ContentType','vector', ...
    %     'Resolution',DIP_IM_DPI);
    %## CORONAL
    % view([0,0])
    %-- save
    % exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_nogroup_coronal.tiff',cluster_i)], ...
    %     'Resolution',DIP_IM_DPI);
    % exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_nogroup_coronal.pdf',cluster_i)], ...
    %     'ContentType','vector', ...
    %     'Resolution',DIP_IM_DPI);
    close(fig);
    %}

    %## ANGLE ALL WITH GROUP SYMBOLS =================================== %%
    DIPPLOT_STRUCT.projlines = 'on';
    DIPPLOT_STRUCT.projcol = {[0,0,0]};
    DIPPLOT_STRUCT.projalpha = 0.2;
    [fig] = eeglab_dipplot(STUDY,ALLEEG,cluster_i,...
        'PLOT_TYPE','all_group',...
        'DIPPLOT_STRUCT',DIPPLOT_STRUCT);
    % camzoom(1.21);
    camzoom(1.01);
    %## ANGULAR
    view([45,30])    
    exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_angle.png',cluster_i)], ...
        'Resolution',DIP_IM_DPI);
    % exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_angle.pdf',cluster_i)], ...
    %     'ContentType','vector', ...
    %     'Resolution',DIP_IM_DPI);
    saveas(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_angle.fig',cluster_i)]);

    %## NON ANGLE PLOT ===============================================
    DIPPLOT_STRUCT.projlines = 'off';
    [fig] = eeglab_dipplot(STUDY,ALLEEG,cluster_i,...
        'PLOT_TYPE','all_group',...
        'DIPPLOT_STRUCT',DIPPLOT_STRUCT);
    camzoom(1.1^2);
    %## TOP
    view([0,90])    
    exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_top.png',cluster_i)], ...
        'Resolution',DIP_IM_DPI);
    % exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_top.pdf',cluster_i)], ...
    %     'ContentType','vector', ...
    %     'Resolution',DIP_IM_DPI);
    saveas(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_top.fig',cluster_i)]);
    %## SAGGITAL
    view([90,0])
    %-- save
    exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_sagittal.png',cluster_i)], ...
        'Resolution',DIP_IM_DPI);
    % exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_sagittal.pdf',cluster_i)], ...
    %     'ContentType','vector', ...
    %     'Resolution',DIP_IM_DPI);
    %## CORONAL
    view([0,0])
    %-- save
    exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_coronal.png',cluster_i)], ...
        'Resolution',DIP_IM_DPI);
    % exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_coronal.pdf',cluster_i)], ...
    %     'ContentType','vector', ...
    %     'Resolution',DIP_IM_DPI);
    pause(2);
    close(fig);
end