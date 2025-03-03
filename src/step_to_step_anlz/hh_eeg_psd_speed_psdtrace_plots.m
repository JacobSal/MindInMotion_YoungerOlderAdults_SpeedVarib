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
%-
save_dir = [cluster_k_dir filesep ANALYSIS_DNAME];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% (LOAD & CONGREGATE PSD DATA) ===================================== %%
% tmp = 'perstridefb_nfslidingb16';
% fextr = 'perstridefb_nfslidingb26';
fextr = 'slidingb6';
% fextr = 'slidingb12';
% fextr = 'slidingb18';
% fextr = 'slidingb24';
f_range = [3, 40];

%## FNAMES
fnames = {STUDY.datasetinfo.filename};
fpaths = {STUDY.datasetinfo.filepath};
subj_chars = {STUDY.datasetinfo.subject};

%## GET FREQ DIMS
dat = par_load(fpaths{1},sprintf('psd_output_%s.mat',fextr));
psd_mean = dat.psd_mean;
%--
psd_dat_out = nan(length(subj_chars),size(psd_mean,1),250,length(CLUSTER_PICS));
cond_dat_out = cell(length(subj_chars),250,length(CLUSTER_PICS));
subj_dat_out = nan(length(subj_chars),250,length(CLUSTER_PICS));
group_dat_out = cell(length(subj_chars),250,length(CLUSTER_PICS));
conds_extract = {'0p25','0p5','0p75','1p0'};

%##
for subj_i = 1:length(subj_chars)
    %-- initiate params
    tmp_study = STUDY;
    tmp_cl_study = CL_STUDY;
    tmp_subjs = {tmp_cl_study.datasetinfo.subject};
    %## LOAD SET FIlE
    EEG = load([fpaths{subj_i} filesep fnames{subj_i}],'-mat');
    
    %## LOAD PSD DATA
    dat = par_load(fpaths{subj_i},sprintf('psd_output_%s.mat',fextr));
    %--
    psd_mean = dat.psd_mean; %[frequency, epoch/splice, channel/component]
    fooof_freqs = dat.freqs;
    cond_struct = dat.cond_struct;
    fprintf('conds: %s\n',strjoin(unique({cond_struct.cond}),','));
    comp_arr = dat.indx_to_comp;
    psd_mean = psd_mean(:,[cond_struct.ind],:);

    %## GET CLUSTER DATA
    fprintf('Getting cluster information...\n');
    [~,subs] = intersect(comp_arr(3,:),main_cl_inds);
    tmp_arr = comp_arr(:,subs);
    
    %## EXTRACT DATA TO CLUSTERED ARRAYS
    fprintf('Extract data and assigning to big array...\n');    
    psd_dat_out(subj_i,:,1:size(psd_mean,2),tmp_arr(1,:)) = psd_mean(:,:,tmp_arr(1,:));
    subj_dat_out(subj_i,1:size(psd_mean,2),tmp_arr(1,:)) = repmat(subj_i,[1,size(psd_mean,2),length(tmp_arr)]);
    tmp = {cond_struct.cond};
    cond_dat_out(subj_i,1:length(cond_struct),tmp_arr(1,:)) = repmat(tmp',[1,length(tmp_arr)]);
    group_dat_out(subj_i,1:length(cond_struct),tmp_arr(1,:)) = repmat({EEG.group},[1,size(psd_mean,2),length(tmp_arr)]);
    fprintf('psd len: %i\ncond len: %i\n',size(psd_mean,2),length(cond_struct));
    % psd_dat_out(subj_i) = psd_mean(:,:,tmp_arr(1,:));
end
dat_out_struct = struct('psd_dat',psd_dat_out, ...
    'cond_dat',{cond_dat_out}, ...
    'subj_dat',subj_dat_out, ...
    'group_dat',{group_dat_out});
par_save(dat_out_struct,[cluster_k_dir filesep 'kin_eeg_step_to_step' filesep sprintf('raw_psd_dat_%s.mat',fextr)]);
dat_out_struct = struct.empty;
%## LOAD
%{
dat_out_struct = par_load([cluster_k_dir filesep 'kin_eeg_step_to_step' filesep sprintf('raw_psd_dat_%s.mat',fextr)]);
psd_dat_out = dat_out_struct.psd_dat;
cond_dat_out = dat_out_struct.cond_dat;
subj_dat_out = dat_out_struct.subj_dat;
group_dat_out = dat_out_struct.group_dat;
%}
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
xtick_label_c = {'0.25 m/s','0.50 m/s','0.75 m/s','1.0 m/s'};
cluster_inds_plot = [3,4,5,7,10,12,13];
% [~,~,cluster_inds] = intersect(cluster_inds_plot,CLUSTER_PICS);
%% ===================================================================== %%
%## PARAMETERS
%-- colors
cmaps_speed = linspecer(4*3);
cmaps_speed = [cmaps_speed(1,:);cmaps_speed(2,:);cmaps_speed(3,:);cmaps_speed(4,:)];
%--
SAVE_RES = 300;
TITLE_TXT_SIZE = 14;
% IM_RESIZE = 0.9;
% AX_W = 0.3;
% AX_H = 0.25;
AX_W = 0.325;
AX_H = 0.25;
IM_RESIZE = 0.7;
AX_FONT_NAME = 'Arial';
AX_X_SHIFT = 1.7;
AX_Y_SHIFT = 1.4;
AX_INIT_X = 0.13;
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

%## PSD PLOTS
PLOT_STRUCT = struct('y_label',{'10*log_{10}(PSD)'},...
    'y_label_fontsize',12,...
    'y_label_fontweight','bold',...
    'ylim',[],...
    'x_label',{'Frequency (Hz)'},...
    'x_label_fontsize',12,...
    'x_label_fontweight','bold',...
    'x_label_yoffset',0,...
    'xtick_labs',{{}}, ...
    'xticks',[], ...
    'xlim',[],...
    'title',{{''}},...
    'title_fontsize',12,...
    'title_fontweight','normal',...
    'font_size',12,...
    'font_name','Arial',...
    'ax_position',[0,0,1,1],...
    'ax_line_width',1,...
    'xtick_angle',45);
LINE_STRUCT = struct('line_width',2, ...
    'line_style','-', ...
    'line_alpha',0.6, ...
    'line_color',[1,1,1], ...
    'line_label',{'label'}, ...
    'line_avg_fcn',@mean, ...
    'line_avg_fcn_params',{{2}}, ...
    'do_err_shading',true, ...
    'err_alpha',0.6, ...
    'err_color',[0.5,0.5,0.5], ...
    'err_edge_color','none', ...
    'err_bnd_fcn',@prctile, ...
    'err_upr_bnd_params',{{75,2}}, ...
    'err_lwr_bnd_params',{{25,2}}, ...
    'err_line_style',':', ...
    'err_line_width',3);
%% (ALL SUBJS MODEL) =================================================== %%
tmp_savedir = [save_dir filesep fextr];
mkdir(tmp_savedir);
%#%--
x_shift = AX_INIT_X;
x_cnt = 1;
y_shift = AX_INIT_Y;
vert_shift = 0;
horiz_shift = 0;
stats_store = [];
for cl_i = 1:length(cluster_inds_plot)
    %%
    %## INITIATE FIGURE
    fig = figure('color','white');
    set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    annotation('textbox',[TITLE_XSHIFT-(TITLE_BOX_SZ(1)/2/2),TITLE_YSHIFT-(TITLE_BOX_SZ(2)/2),TITLE_BOX_SZ(1),TITLE_BOX_SZ(2)],...
        'String',cluster_titles{double(string(cluster_inds_plot(cl_i)))}, ...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle', ...
        'LineStyle','none', ...
        'FontName',AX_FONT_NAME,...
        'FontSize',TITLE_FONT_SIZE, ...
        'FontWeight','Bold', ...
        'Units','normalized');        
    hold on;
    set(gca,AXES_DEFAULT_PROPS{:})
    
    %## EXTRACT PSD DATA
    tmp_dat = squeeze(psd_dat_out(:,:,:,cluster_inds_plot(cl_i))); %[subject, frequency, epoch/splice, channel/component]
    tmp_dat = reshape(permute(tmp_dat,[3,1,2]),size(tmp_dat,1)*size(tmp_dat,3),size(tmp_dat,2)); %[subject x epoch/splice, frequency];
    chk = all(~isnan(tmp_dat),2);
    tmp_dat = tmp_dat(chk,:);
    %--
    tmp_cond = squeeze(cond_dat_out(:,:,cluster_inds_plot(cl_i))); %[subject, frequency, epoch/splice, channel/component]    
    tmp_cond = reshape(permute(tmp_cond,[2,1]),[size(cond_dat_out,1)*size(cond_dat_out,2),1]);
    chk = cellfun(@isempty,tmp_cond);
    tmp_cond = tmp_cond(~chk);
    conds = unique(tmp_cond);
    % %--
    % tmp_cond = cat(2,cond_dat_out{:});reshape(cond_dat_out,[size(cond_dat_out,1)*size(cond_dat_out,2),1]);
    % chk = cellfun(@isempty,tmp_cond);
    % tmp_cond = tmp_cond(~chk);
    % conds = unique(tmp_cond);
    %--
    tmp_subj = squeeze(subj_dat_out(:,:,cluster_inds_plot(cl_i))); %[subject, frequency, epoch/splice, channel/component]    
    tmp_subj = reshape(permute(tmp_subj,[2,1]),[size(subj_dat_out,1)*size(subj_dat_out,2),1]);
    chk = ~all(isnan(tmp_subj),2);
    tmp_subj = tmp_subj(chk,:);
    subjs = unique(tmp_subj);
    %## STATISTICS
    %-- Ho : all samples come from the same distribution
    %-- Ha : all samples come from different distributions
    tmp_psd_in = cell(length(conds),1);
    for c_i = 1:length(conds)
        indc = strcmp(tmp_cond,conds{c_i});
        tmp = zeros(size(tmp_dat,2),length(subjs));
        for s_i = 1:length(subjs)
            inds = tmp_subj == subjs(s_i);
            chk = indc & inds;
            tmp(:,s_i) = mean(tmp_dat(chk,:),1);
        end
        tmp = tmp(:,all(tmp ~= 0,1));
        tmp_psd_in{c_i,1} = tmp; %tmp_dat(inds,:);
    end
    [temp_pcond, temp_pgroup, temp_pinter, temp_statcond, temp_statgroup, temp_statinter] = ...
        std_stat(tmp_psd_in, stats);    
    temp_pcond=temp_pcond{1} < 0.05;
    
    %## PLOT
    ax = axes();
    IM_RESIZE = 1.1;
    % AX_W = 0.8;
    % AX_H = 0.7;
    for c_i = 1:length(conds)
        inds = strcmp(tmp_cond,conds{c_i});
        %--
        PLOT_STRUCT.title = {sprintf('Cluster %i',cluster_inds_plot(cl_i))};
        PLOT_STRUCT.ax_position = [x_shift,y_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE];
        PLOT_STRUCT.xlim = [3,40];
        PLOT_STRUCT.ylim = sort([prctile([tmp_psd_in{:}],98,'all'),prctile([tmp_psd_in{:}],2,'all')]);
        %--
        LINE_STRUCT.line_color = cmaps_speed(c_i,:);
        LINE_STRUCT.err_color = cmaps_speed(c_i,:)+0.15;
        LINE_STRUCT.line_alpha = 0.7;
        LINE_STRUCT.err_alpha = 0.3;
        LINE_STRUCT.line_label = xtick_label_c{c_i};
        LINE_STRUCT.do_err_shading = true;
        LINE_STRUCT.err_upr_bnd_params = {97.5,2};
        LINE_STRUCT.err_lwr_bnd_params = {2.5,2};
        %--
        plot_psd(ax,tmp_psd_in{c_i,1},fooof_freqs, ...
            'LINE_STRUCT',LINE_STRUCT, ...
            'PLOT_STRUCT',PLOT_STRUCT);
        hold on;
        % for s_i = 1:size(tmp_psd_in{c_i,1},2)
        %     LINE_STRUCT.do_err_shading = false;
        %     plot_psd(ax,tmp_psd_in{c_i,1},fooof_freqs, ...
        %     'LINE_STRUCT',LINE_STRUCT, ...
        %     'PLOT_STRUCT',PLOT_STRUCT);
        %     hold on;
        % end
    end
    [axsignif,Pa] = plot_psd_stats(ax,fooof_freqs,temp_pcond, ...
        'background','Frequency (Hz)');    

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
    hold off;
    %##
    % exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_std_allsubj.tiff',string(clusters(c_i)))],...
    %     'Resolution',300)
    exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_%s_psd_trace.tiff',string(clusters(c_i)),meas_ext)],...
        'Resolution',SAVE_RES)
    % close(fig)
end