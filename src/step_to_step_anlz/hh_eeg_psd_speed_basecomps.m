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
cmap_terrain = linspecer(4);
custom_yellow = [254,223,0]/255;
cmap_terrain = [cmap_terrain(3,:);custom_yellow;cmap_terrain(4,:);cmap_terrain(2,:)];
cmap_speed = linspecer(4*3);
cmap_speed = [cmap_speed(1,:);cmap_speed(2,:);cmap_speed(3,:);cmap_speed(4,:)];
%% (PATHS) ============================================================= %%
%- datset name
DATA_SET = 'MIM_dataset';
%- study name
% STUDY_DNAME = '10172024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej_speed';
STUDY_DNAME = '02202025_mim_yaoa_powpow0p3_crit_speed';
STUDY_FNAME = 'kin_eeg_epoch_study';
ANALYSIS_DNAME = 'kin_eeg_step_to_step';
%-
studies_fpath = [PATHS.data_dir filesep DATA_SET filesep '_studies'];
%- load cluster
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
% cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb10'];
% cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb5'];
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

%%
fpaths = {STUDY.datasetinfo.filepath};
fextr = 'slidingb36';
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
% cluster_inds_plot = [3,4,5,7,10,12,13];
cluster_inds_plot = [3,5,12,13,7];
% [~,~,cluster_inds] = intersect(cluster_inds_plot,CLUSTER_PICS);
%% ===================================================================== %%
%## PARAMETERS
%--
SAVE_RES = 600;
TITLE_TXT_SIZE = 14;
% IM_RESIZE = 0.9;
% AX_W = 0.3;
% AX_H = 0.25;
AX_W = 0.325;
AX_H = 0.25;
IM_RESIZE = 0.8;
AX_FONT_NAME = 'Arial';
AX_X_SHIFT = 1.7;
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
LEG_X_SHIFT = -0.01; %-0.1
LEG_Y_SHIFT =  -1.6; %-0.38
LEG_TXT_SIZE = 9;
LEG_TOKEN_SIZE = 15;
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
    'line_alpha',1, ...
    'line_color',[1,1,1], ...
    'line_label',{'label'}, ...
    'line_avg_fcn',@mean, ...
    'line_avg_fcn_params',{{2}}, ...
    'do_err_shading',true, ...
    'err_alpha',0.6, ...
    'err_color',[0.5,0.5,0.5], ...
    'err_edge_color',[], ...
    'err_upr_bnd_fcn',@(x,p) 1*std(x,p{:}), ...
    'err_lwr_bnd_fcn',@(x,p) -1*std(x,p{:}), ...
    'err_upr_bnd_params',{{[],2}}, ...
    'err_lwr_bnd_params',{{[],2}}, ...
    'err_line_style',':', ...
    'err_line_width',3);
%% (CLUSTER FIGURE FOR BASELINE CONDS) ============================== %%
meas_ext = 'clusterwise_baseline_comps';
tmp_savedir = [save_dir filesep meas_ext];
mkdir(tmp_savedir);
leg_chars = {'0.25 m/s','0.50 m/s','0.75 m/s','1.0 m/s'}; %{'b6','b12','b18','b24'};

% leg_store = [];
%##
%-- colors
cmaps_speed = linspecer(4*3);
cmaps_speed = [cmaps_speed(1,:);cmaps_speed(2,:);cmaps_speed(3,:);cmaps_speed(4,:)];
%## EXTRACT PSD DATA
% for d_i = 2
g_i = 1;
for d_i = 1:length(dat_out_structs)
    %## INITIATE FIGURE
    x_shift = AX_INIT_X;
    x_cnt = 1;
    y_shift = AX_INIT_Y;
    leg_store = [];
    %--
    fig = figure('color','white');
    set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    annotation('textbox',[TITLE_XSHIFT-(TITLE_BOX_SZ(1)/2/2),TITLE_YSHIFT-(TITLE_BOX_SZ(2)/2),TITLE_BOX_SZ(1),TITLE_BOX_SZ(2)],...
        'String','Variation Comparisons', ...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle', ...
        'LineStyle','none', ...
        'FontName',AX_FONT_NAME,...
        'FontSize',TITLE_FONT_SIZE, ...
        'FontWeight','Bold', ...
        'Units','normalized');        
    hold on;
    set(gca,AXES_DEFAULT_PROPS{:})
    ax_store = [];
    psd_avg_char = [];
    y_lim_store = zeros(length(cluster_inds_plot),2);
    ycnt = 1;
    for cl_i = 1:length(cluster_inds_plot)        
        % psd_dat_out1 = dat_out_structs{d_i}.psd_dat;        
        % psd_dat_out2 = dat_out_structs{d_i}.psd_std_dat;
        % psd_dat_out = psd_dat_out2./psd_dat_out1;
        %--
        psd_dat_out = dat_out_structs{d_i}.psd_dat;        
        % psd_dat_out = dat_out_structs{d_i}.psd_std_dat;
        if cl_i == 1
            psd_avg_char = [psd_avg_char,'mean'];
            % psd_avg_char = [psd_avg_char,'std'];
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
        conds = unique(tmp_group);
        groups = unique(tmp_group);
        indsg = strcmp(tmp_group,groups{g_i});
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
            % LINE_STRUCT.line_avg_fcn = @(x,p) std(x,p{:});
            % LINE_STRUCT.line_avg_fcn_params = {[],2};
            LINE_STRUCT.line_avg_fcn = @(x,p) mean(x,p{:});
            LINE_STRUCT.line_avg_fcn_params = {2};
            LINE_STRUCT.line_color = cmaps_speed(c_i,:);
            LINE_STRUCT.line_alpha = 0.7;
            LINE_STRUCT.line_label = xtick_label_c{c_i};
            %--
            LINE_STRUCT.do_err_shading = true;
            LINE_STRUCT.err_color = cmaps_speed(c_i,:)+0.15;
            LINE_STRUCT.err_alpha = 0.3;
            LINE_STRUCT.err_upr_bnd_fcn = @(x,p) mean(x,2) + std(x,p{:});
            LINE_STRUCT.err_lwr_bnd_fcn = @(x,p) mean(x,2) - std(x,p{:});
            LINE_STRUCT.err_upr_bnd_params = {[],2};
            LINE_STRUCT.err_lwr_bnd_params = {[],2};
            % LINE_STRUCT.err_upr_bnd_fcn = @(x,p) prctile(x,p{:});
            % LINE_STRUCT.err_lwr_bnd_fcn = @(x,p) prctile(x,p{:});
            % LINE_STRUCT.err_upr_bnd_params = {25,2};
            % LINE_STRUCT.err_lwr_bnd_params = {75,2};
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

    %--
    y_lim_store = [min(y_lim_store,[],'all'),max(y_lim_store,[],'all')];
    for aa = 1:length(ax_store)
        set(ax_store(aa),'YLim',y_lim_store)
    end
    %## LEGEND
    %- lg2                
    legend(gca,leg_store);
    [lg2,icons,plots,txt]  = legend('boxoff');
    tmp = get(lg2,'String');
    cnt = 1;
    for i = 1:length(leg_store)
        tmp{i} = sprintf('%s',leg_chars{cnt});
        cnt = cnt + 1;
    end
    set(lg2,'String',tmp,'FontName',AX_FONT_NAME,'FontSize',LEG_TXT_SIZE)
    set(lg2,'Orientation','horizontal')
    set(lg2,'Units','normalized')
    set(lg2,'Position',[AX_INIT_X+LEG_X_SHIFT*IM_RESIZE*AX_W,...
        y_shift+AX_H*IM_RESIZE+LEG_Y_SHIFT*IM_RESIZE*AX_H,lg2.Position(3),lg2.Position(4)]);
    lg2.ItemTokenSize(1) = LEG_TOKEN_SIZE;
    hold off;
    %##
    % exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_std_allsubj.tiff',string(clusters(c_i)))],...
    %     'Resolution',300)
    exportgraphics(fig,[tmp_savedir filesep sprintf('allclusts_base_%s_%s.tiff',psd_avg_char,basel_chars{d_i})],...
        'Resolution',SAVE_RES)
    % close(fig)
end
%% (PER GROUP CLUSTER FIGURE FOR BASELINE CONDS) ===================== %%
meas_ext = 'clusterwise_baseline_comps';
tmp_savedir = [save_dir filesep meas_ext];
mkdir(tmp_savedir);
leg_chars = {'0.25 m/s','0.50 m/s','0.75 m/s','1.0 m/s'}; %{'b6','b12','b18','b24'};
groups_proc = {'H1000','H2000','H3000'};
% leg_store = [];
%##
%-- colors
cmaps_speed = linspecer(4*3);
cmaps_speed = [cmaps_speed(1,:);cmaps_speed(2,:);cmaps_speed(3,:);cmaps_speed(4,:)];
%## EXTRACT PSD DATA
% for d_i = 2
for g_i = 1:length(groups_proc)
    for d_i = 1:length(dat_out_structs)
        %## INITIATE FIGURE
        x_shift = AX_INIT_X;
        x_cnt = 1;
        y_shift = AX_INIT_Y;
        leg_store = [];
        %--
        fig = figure('color','white');
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        annotation('textbox',[TITLE_XSHIFT-(TITLE_BOX_SZ(1)/2/2),TITLE_YSHIFT-(TITLE_BOX_SZ(2)/2),TITLE_BOX_SZ(1),TITLE_BOX_SZ(2)],...
            'String','Variation Comparisons', ...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle', ...
            'LineStyle','none', ...
            'FontName',AX_FONT_NAME,...
            'FontSize',TITLE_FONT_SIZE, ...
            'FontWeight','Bold', ...
            'Units','normalized');        
        hold on;
        set(gca,AXES_DEFAULT_PROPS{:})
        ax_store = [];
        psd_avg_char = [];
        y_lim_store = zeros(length(cluster_inds_plot),2);
        ycnt = 1;
        for cl_i = 1:length(cluster_inds_plot)        
            % psd_dat_out1 = dat_out_structs{d_i}.psd_dat;        
            % psd_dat_out2 = dat_out_structs{d_i}.psd_std_dat;
            % psd_dat_out = psd_dat_out2./psd_dat_out1;
            %--
            psd_dat_out = dat_out_structs{d_i}.psd_dat;        
            % psd_dat_out = dat_out_structs{d_i}.psd_std_dat;
            if cl_i == 1
                psd_avg_char = [psd_avg_char,'mean'];
                % psd_avg_char = [psd_avg_char,'std'];
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
            groups = unique(tmp_group);
            %--
            indsg = strcmp(tmp_group,groups{g_i});
            tmp_group = tmp_group(indsg,:);
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
                    % tmp(:,s_i) = mean(tmp_dat(chk,:),1);
                    tmp(:,s_i) = std(tmp_dat(chk,:),[],1);        
                    % tmp(:,s_i) = prctile(tmp_dat(chk,:),75,1) - prctile(tmp_dat(chk,:),25,1);
                    if c_ii == 1 && s_i == 1 && cl_i == 1
                        % psd_avg_char = [psd_avg_char,'mean'];
                        psd_avg_char = [psd_avg_char,'std'];
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
                % LINE_STRUCT.line_avg_fcn = @(x,p) std(x,p{:});
                % LINE_STRUCT.line_avg_fcn_params = {[],2};
                LINE_STRUCT.line_avg_fcn = @(x,p) mean(x,p{:});
                LINE_STRUCT.line_avg_fcn_params = {2};
                LINE_STRUCT.line_color = cmaps_speed(c_i,:);
                LINE_STRUCT.line_alpha = 0.7;
                LINE_STRUCT.line_label = xtick_label_c{c_i};
                %--
                LINE_STRUCT.do_err_shading = true;
                LINE_STRUCT.err_color = cmaps_speed(c_i,:)+0.15;
                LINE_STRUCT.err_alpha = 0.3;
                LINE_STRUCT.err_upr_bnd_fcn = @(x,p) mean(x,2) + std(x,p{:});
                LINE_STRUCT.err_lwr_bnd_fcn = @(x,p) mean(x,2) - std(x,p{:});
                LINE_STRUCT.err_upr_bnd_params = {[],2};
                LINE_STRUCT.err_lwr_bnd_params = {[],2};
                % LINE_STRUCT.err_upr_bnd_fcn = @(x,p) prctile(x,p{:});
                % LINE_STRUCT.err_lwr_bnd_fcn = @(x,p) prctile(x,p{:});
                % LINE_STRUCT.err_upr_bnd_params = {25,2};
                % LINE_STRUCT.err_lwr_bnd_params = {75,2};
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
    
        %--
        y_lim_store = [min(y_lim_store,[],'all'),max(y_lim_store,[],'all')];
        for aa = 1:length(ax_store)
            set(ax_store(aa),'YLim',y_lim_store)
        end
        %## LEGEND
        %- lg2                
        legend(gca,leg_store);
        [lg2,icons,plots,txt]  = legend('boxoff');
        tmp = get(lg2,'String');
        cnt = 1;
        for i = 1:length(leg_store)
            tmp{i} = sprintf('%s',leg_chars{cnt});
            cnt = cnt + 1;
        end
        set(lg2,'String',tmp,'FontName',AX_FONT_NAME,'FontSize',LEG_TXT_SIZE)
        set(lg2,'Orientation','horizontal')
        set(lg2,'Units','normalized')
        set(lg2,'Position',[AX_INIT_X+LEG_X_SHIFT*IM_RESIZE*AX_W,...
            y_shift+AX_H*IM_RESIZE+LEG_Y_SHIFT*IM_RESIZE*AX_H,lg2.Position(3),lg2.Position(4)]);
        lg2.ItemTokenSize(1) = LEG_TOKEN_SIZE;
        hold off;
        %##
        % exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_std_allsubj.tiff',string(clusters(c_i)))],...
        %     'Resolution',300)
        exportgraphics(fig,[tmp_savedir filesep sprintf('g%s_base_%s_%s.tiff',groups{g_i},psd_avg_char,basel_chars{d_i})],...
            'Resolution',SAVE_RES)
        close(fig)
    end
end

%% (BASELINE COMPS PER CONDITION) ====================================== %%
meas_ext = 'baseline_comps_percond';
tmp_savedir = [save_dir filesep meas_ext];
mkdir(tmp_savedir);
leg_chars = basel_chars; %{'b6','b12','b18','b24'};
leg_store = [];
%##
cmaps_speed = linspecer(length(dat_out_structs));
x_shift = AX_INIT_X;
x_cnt = 1;
y_shift = AX_INIT_Y;
cl_i = 2;

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
for c_i = 1:length(conds)
    %## PLOT
    ax = axes();
    tmp_psd_in = cell(length(dat_out_structs),1);
    hold on;
    for d_i = 1:length(dat_out_structs)
        psd_dat_out = dat_out_structs{d_i}.psd_dat;
        cond_dat_out = dat_out_structs{d_i}.cond_dat;
        subj_dat_out = dat_out_structs{d_i}.subj_dat;
        % group_dat_out = dat_out_structs{d_i}.group_dat;
        tmp_dat = squeeze(psd_dat_out(:,:,:,cluster_inds_plot(cl_i))); %[subject, frequency, epoch/splice, channel/component]
        tmp_dat = reshape(permute(tmp_dat,[3,1,2]),size(tmp_dat,1)*size(tmp_dat,3),size(tmp_dat,2)); %[subject x epoch/splice, frequency];
        chk = all(~isnan(tmp_dat),2);
        tmp_dat = tmp_dat(chk,:);
        sum(chk)
        %--
        tmp_cond = squeeze(cond_dat_out(:,:,cluster_inds_plot(cl_i))); %[subject, frequency, epoch/splice, channel/component]    
        tmp_cond = reshape(permute(tmp_cond,[2,1]),[size(cond_dat_out,1)*size(cond_dat_out,2),1]);
        chk = cellfun(@isempty,tmp_cond);
        tmp_cond = tmp_cond(~chk);
        conds = unique(tmp_cond);
        %--
        tmp_subj = squeeze(subj_dat_out(:,:,cluster_inds_plot(cl_i))); %[subject, frequency, epoch/splice, channel/component]    
        tmp_subj = reshape(permute(tmp_subj,[2,1]),[size(subj_dat_out,1)*size(subj_dat_out,2),1]);
        chk = ~all(isnan(tmp_subj),2);
        tmp_subj = tmp_subj(chk,:);
        subjs = unique(tmp_subj);
        %##
        indc = strcmp(tmp_cond,conds{c_i});
        tmp = zeros(size(tmp_dat,2),length(subjs));
        for s_i = 1:length(subjs)
            inds = tmp_subj == subjs(s_i);
            chk = indc & inds;
            tmp(:,s_i) = mean(tmp_dat(chk,:),1);
        end
        tmp = tmp(:,all(tmp ~= 0,1));
        tmp_psd_in{d_i,1} = tmp; %tmp_dat(inds,:);
        
    end
    %## STATISTICS
    %-- Ho : all samples come from the same distribution
    %-- Ha : all samples come from different distributions
    % tmp_psd_in = cell(length(conds),1);
    % for c_ii = 1:length(conds)
    %     indc = strcmp(tmp_cond,conds{c_ii});
    %     tmp = zeros(size(tmp_dat,2),length(subjs));
    %     for s_i = 1:length(subjs)
    %         inds = tmp_subj == subjs(s_i);
    %         chk = indc & inds;
    %         tmp(:,s_i) = mean(tmp_dat(chk,:),1);
    %     end
    %     tmp = tmp(:,all(tmp ~= 0,1));
    %     tmp_psd_in{c_ii,1} = tmp; %tmp_dat(inds,:);
    % end
    % [pcond, pgroup, pinter, statcond, statgroup, statinter] = ...
    %     std_stat(tmp_psd_in, stats);    
    % pcond=pcond{1} < 0.05;    

    %## PLOT
    for d_i = 1:length(dat_out_structs)
        PLOT_STRUCT.title = {sprintf('Condition %s for Cluster %i',conds{c_i},cluster_inds_plot(cl_i))};
        PLOT_STRUCT.ax_position = [x_shift,y_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE];
        PLOT_STRUCT.xlim = [3,40];
        PLOT_STRUCT.ylim = [-1.5,3]; %sort([prctile([tmp_psd_in{:}],99,'all'),prctile([tmp_psd_in{:}],1,'all')]);
        disp(PLOT_STRUCT.ylim);
        %--
        LINE_STRUCT.line_color
        % LINE_STRUCT.line_color = cmaps_speed(c_i,:);
        % LINE_STRUCT.err_color = cmaps_speed(c_i,:);
        LINE_STRUCT.line_color = cmaps_speed(d_i,:);
        LINE_STRUCT.err_color = cmaps_speed(d_i,:);
        LINE_STRUCT.line_alpha = 0.7;
        LINE_STRUCT.err_alpha = 0.3;
        LINE_STRUCT.line_label = xtick_label_c{c_i};
        LINE_STRUCT.do_err_shading = true;
        LINE_STRUCT.err_bnd_fcn = @std;
        LINE_STRUCT.err_upr_bnd_params = {[],2};
        LINE_STRUCT.err_lwr_bnd_params = {[],2};
        % LINE_STRUCT.err_edge_color = cmaps_speed(c_i,:);
        LINE_STRUCT.err_edge_color = cmaps_speed(d_i,:);
        LINE_STRUCT.err_line_width = 1;
        %--
        [ax,Pa,Li] = plot_psd(ax,tmp_psd_in{d_i,1},fooof_freqs, ...
            'LINE_STRUCT',LINE_STRUCT, ...
            'PLOT_STRUCT',PLOT_STRUCT);
        if c_i == 1
            leg_store = [leg_store, Li];
        end
        hold on;
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
end
%## LEGEND
%- lg2                
legend(gca,leg_store);
[lg2,icons,plots,txt]  = legend('boxoff');
tmp = get(lg2,'String');
cnt = 1;
for i = 1:length(leg_store)
    tmp{i} = sprintf('%s',leg_chars{cnt});
    cnt = cnt + 1;
end
set(lg2,'String',tmp,'FontName',AX_FONT_NAME,'FontSize',LEG_TXT_SIZE)
set(lg2,'Orientation','horizontal')
set(lg2,'Units','normalized')
set(lg2,'Position',[AX_INIT_X+LEG_X_SHIFT*IM_RESIZE,...
    y_shift+AX_H*IM_RESIZE+LEG_Y_SHIFT*IM_RESIZE,lg2.Position(3),lg2.Position(4)]);
lg2.ItemTokenSize(1) = LEG_TOKEN_SIZE;
hold off;
%##
% exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_std_allsubj.tiff',string(clusters(c_i)))],...
%     'Resolution',300)
exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_%s_cond_comps.tiff',string(cluster_inds_plot(cl_i)),meas_ext)],...
    'Resolution',SAVE_RES)
% close(fig)
%% (SUBJ. ACROSS TIME PLOTS) =========================================== %%
meas_ext = 'baseline_comps';
tmp_savedir = [save_dir filesep meas_ext];
mkdir(tmp_savedir);
leg_chars = {'b6','b12','b18','b24'};
%##
cmaps_speed = linspecer(length(dat_out_structs));
x_shift = AX_INIT_X;
x_cnt = 1;
y_shift = AX_INIT_Y;
vert_shift = 0;
horiz_shift = 0;
stats_store = [];
cl_i = 2;

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
c_i = 1;
d_i = 1;
s_i = 1;
%## PLOT
ax = axes();
leg_store = [];
tmp_psd_in = cell(length(dat_out_structs),1);
psd_dat_out = dat_out_structs{d_i}.psd_dat;
cond_dat_out = dat_out_structs{d_i}.cond_dat;
subj_dat_out = dat_out_structs{d_i}.subj_dat;
group_dat_out = dat_out_structs{d_i}.group_dat;
tmp_dat = squeeze(psd_dat_out(:,:,:,cluster_inds_plot(cl_i))); %[subject, frequency, epoch/splice, channel/component]
tmp_dat = reshape(permute(tmp_dat,[3,1,2]),size(tmp_dat,1)*size(tmp_dat,3),size(tmp_dat,2)); %[subject x epoch/splice, frequency];
chk = all(~isnan(tmp_dat),2);
tmp_dat = tmp_dat(chk,:);
sum(chk)
%--
tmp_cond = squeeze(cond_dat_out(:,:,cluster_inds_plot(cl_i))); %[subject, frequency, epoch/splice, channel/component]    
tmp_cond = reshape(permute(tmp_cond,[2,1]),[size(cond_dat_out,1)*size(cond_dat_out,2),1]);
chk = cellfun(@isempty,tmp_cond);
tmp_cond = tmp_cond(~chk);
conds = unique(tmp_cond);
%--
tmp_subj = squeeze(subj_dat_out(:,:,cluster_inds_plot(cl_i))); %[subject, frequency, epoch/splice, channel/component]    
tmp_subj = reshape(permute(tmp_subj,[2,1]),[size(subj_dat_out,1)*size(subj_dat_out,2),1]);
chk = ~all(isnan(tmp_subj),2);
tmp_subj = tmp_subj(chk,:);
subjs = unique(tmp_subj);
%--
indc = strcmp(tmp_cond,conds{c_i});
inds = tmp_subj == subjs(s_i);
chk = indc & inds;
tmp_psd_in = tmp_dat(chk,:);

%## STATISTICS
% [pcond, pgroup, pinter, statcond, statgroup, statinter] = std_stat(tmp_psd_in, stats);    
% pcond=pcond{1} < 0.05;

%## PLOT
PLOT_STRUCT.title = {sprintf('Cluster %i',cluster_inds_plot(cl_i))};
PLOT_STRUCT.ax_position = [x_shift,y_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE];
PLOT_STRUCT.xlim = [3,40];
PLOT_STRUCT.ylim = [-1.5,3]; %sort([prctile([tmp_psd_in{:}],99,'all'),prctile([tmp_psd_in{:}],1,'all')]);
disp(PLOT_STRUCT.ylim);
%--
LINE_STRUCT.line_color
% LINE_STRUCT.line_color = cmaps_speed(c_i,:);
% LINE_STRUCT.err_color = cmaps_speed(c_i,:);
LINE_STRUCT.line_color = cmaps_speed(d_i,:);
LINE_STRUCT.err_color = cmaps_speed(d_i,:);
LINE_STRUCT.line_alpha = 0.7;
LINE_STRUCT.err_alpha = 0.3;
LINE_STRUCT.line_label = xtick_label_c{c_i};
LINE_STRUCT.do_err_shading = true;
LINE_STRUCT.err_upr_bnd_params = {75,2};
LINE_STRUCT.err_lwr_bnd_params = {25,2};
% LINE_STRUCT.err_edge_color = cmaps_speed(c_i,:);
LINE_STRUCT.err_edge_color = cmaps_speed(d_i,:);
LINE_STRUCT.err_line_width = 1;
%--
[ax] = plot_psd(ax,tmp_psd_in{c_i,1},fooof_freqs, ...
    'LINE_STRUCT',LINE_STRUCT, ...
    'PLOT_STRUCT',PLOT_STRUCT);
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
%- lg2                
legend(gca,leg_store);
[lg2,icons,plots,txt]  = legend('boxoff');
tmp = get(lg2,'String');
cnt = 1;
for i = 1:length(leg_store)
    tmp{i} = sprintf('%0.2g m/s',double(string(leg_chars(cnt))));
    cnt = cnt + 1;
end
set(lg2,'String',tmp,'FontName',AX_FONT_NAME,'FontSize',LEG_TXT_SIZE)
set(lg2,'Orientation','horizontal')
set(lg2,'Units','normalized')
set(lg2,'Position',[AX_INIT_X+LEG_X_SHIFT*IM_RESIZE,...
    y_shift+AX_H*IM_RESIZE+LEG_Y_SHIFT*IM_RESIZE,lg2.Position(3),lg2.Position(4)]);
lg2.ItemTokenSize(1) = LEG_TOKEN_SIZE;
hold off;
%##
% exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_std_allsubj.tiff',string(clusters(c_i)))],...
%     'Resolution',300)
exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_%s_cond_comps.tiff',string(cluster_inds_plot(cl_i)),meas_ext)],...
    'Resolution',SAVE_RES)
% close(fig)
