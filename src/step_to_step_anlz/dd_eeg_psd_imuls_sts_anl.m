%   Project Title: MIM OA & YA SPEED & KINETICS ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_oa_speed_eeg_out/run_a_epoch_process.sh

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
% global ADD_CLEANING_SUBMODS STUDY_DIR SCRIPT_DIR %#ok<GVMIS>
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
ANALYSIS_DNAME = 'kin_eeg_step_to_step';
%-
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
%- load cluster
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
% cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb10'];
cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb5'];
% cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb3'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
cluster_k_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
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
%% (STUDY DESIGN) ====================================================== %%
ERSP_STAT_PARAMS_COND = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','cluster',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
STUDY_DESI_PARAMS = {{'subjselect',{},...
    'variable1','cond','values1',{'flat','low','med','high'},'pairing','on',...
    'variable2','group','values2',{'H1000''s','H2000''s','H3000''s'}},...
    {'subjselect',{},...
    'variable1','cond','values1',{'0p25','0p5','0p75','1p0'},'pairing','on',...
    'variable2','group','values2',{'H1000''s','H2000''s','H3000''s'}}};
% condition_gait = unique({STUDY.datasetinfo(1).trialinfo.cond}); 
% % (09/18/2024) JS, reorders it weirdly. Manually override
%## ersp plot per cluster per condition
args = eeglab_struct2args(ERSP_STAT_PARAMS_COND);
STUDY = pop_statparams(STUDY,args{:});
args = eeglab_struct2args(ERSP_PARAMS);
STUDY = pop_erspparams(STUDY,args{:});
STUDY.cache = [];
for des_i = 1:length(STUDY_DESI_PARAMS)
    [STUDY] = std_makedesign(STUDY,[],des_i,STUDY_DESI_PARAMS{des_i}{:});
end
%% (SUBJECT LOOP)
%## PARAMETERS
des_i = 2;
speed_n_chars = {'0.25','0.50','0.75','1.0'};
conds = STUDY.design(des_i).variable(1).value; %event.cond});
% conds = conds(1:4);
theta_band_lims = [4, 8];
alpha_band_lims = [8 12];
beta_band_lims  = [12 30];
alpha1_band_lims = [8,10.5];
alpha2_band_lims = [10.5,13];
beta1_band_lims = [13,20];
beta2_band_lims = [20,30];
TIME_FACTOR = 1/1000;
def_step_structs = struct('subj_char',{''}, ...
        'cond_char',{''}, ...
        'speed_char',{''}, ...
        'speed_n',[], ...
        'model_char',{''}, ...
        'model_n',[], ...
        'comp_n',[], ...
        'cluster_n',[], ...
        'comp_n_old',[], ...
        'stride_n',[], ...
        'rhs_1',[], ...
        'lto_1',[], ...
        'lhs_1',[], ...
        'rto_1',[], ...
        'rhs_2',[], ...
        'eeg_psd',[], ...
        'eeg_freqs',[], ...
        'step_dur',[], ...
        'gait_cycle_dur',[], ...
        'stance_dur',[], ...
        'swing_dur',[], ...
        'double_sup_dur',[], ...
        'single_sup_dur',[], ...
        'avg_theta',[], ...
        'avg_alpha',[], ...
        'avg_beta',[]);
cmaps_speed = linspecer(4*3);
cmaps_speed = [cmaps_speed(1,:);cmaps_speed(2,:);cmaps_speed(3,:);cmaps_speed(4,:)];
%##
tmp_alleeg = cell(length(STUDY.datasetinfo));
conds_keep = {'0p25','0p5','0p75','1p0'};
xtick_label_g = {'0.25','0.50','0.75','1.0'}; %{'0.25','0.50','0.75','1.0'};
fname_ext_store = '';
%% ===================================================================== %%
for subj_i = 1:length(STUDY.datasetinfo)
% parfor subj_i = 1:length(STUDY.datasetinfo)
    tmp_study = STUDY;
    tmp_cl_study = CL_STUDY;
    fooof_freqs = [];
    ff = 1;
    % fname_ext = cell(1,5);
    %## LOAD EEG DATA
    try
        EEG = pop_loadset('filepath',tmp_study.datasetinfo(subj_i).filepath,'filename',tmp_study.datasetinfo(subj_i).filename);
        subj_ind = strcmp(EEG.subject,{tmp_cl_study.datasetinfo.subject});
        fprintf('Running subject %s...\n',EEG.subject);             
    catch e
        fprintf('%s',getReport(e));
        continue
    end
    if ~any(subj_ind)
        continue
    end    
    %- sanity check
    %{
    pop_eegplot( EEG, 1, 1, 1);
    %}
    %## LOAD PSD DATA
    % epoched_fPath = strsplit(EEG.filepath,filesep);
    epoched_fPath = strsplit(tmp_cl_study.datasetinfo(subj_ind).filepath,filesep);
    icatimef_f = [strjoin(epoched_fPath,filesep) filesep sprintf('%s.icaspec',EEG.subject)];
    %- load .icatimef load-in parameters
    fprintf('Loading Time-Frequency Data...\n');
    tmpf = load(icatimef_f,'-mat');
    %-
    fn = fieldnames(tmpf);
    inds = find(contains(fieldnames(tmpf),'comp'));
    test = tmpf.(fn{inds(1)});
    %-
    eeg_psd = zeros(size(test,1),size(test,2),length(inds),'double');
    for i = 1:length(inds)
        eeg_psd(:,:,i) = tmpf.(fn{inds(i)}); % freq x epoch x chan
    end
    %-    
    freqs_orig = tmpf.freqs;
    trialinfo = tmpf.trialinfo;
    %- subset data based on condition
    inds_cond = cellfun(@(x) any(strcmp(x,conds_keep)),{trialinfo.cond});
    trialinfo = trialinfo(inds_cond);
    eeg_psd = eeg_psd(:,inds_cond,:);
    nolog_eeg_psd = 10.^(eeg_psd/10);
    %## SANITY CHECK
    %{
    i = 1;
    inds = randi([1,size(eeg_psd,2)],50);
    figure;
    hold on;
    plot(freqs,squeeze(eeg_psd(:,inds,i)),'DisplayName','original avg.');
    ylabel('10*log_{10}(PSD)');
    xlabel('Frequency (Hz)');
    hold off;
    %}

    %## RUN FOOOF
    fprintf('Running FOOOF...\n');
    settings = struct('peak_width_limits',[1,8],...
        'min_peak_height',0.05,...
        'max_n_peaks',5);
    f_range = [3, 40];
    f_ind = find(freqs_orig > f_range(1) & freqs_orig < f_range(2));
    f_ind = sort([f_ind; min(f_ind)-1;max(f_ind)+1]);
    return_model = true;
    tmp = zeros(length(f_ind),size(eeg_psd,2),size(eeg_psd,3));
    %-
    % tmp_conds = unique({trialinfo.trial_num_code});
    % tmp_conds = tmp_conds(~cellfun(@isempty,tmp_conds));
    %-
    tmp_conds = unique({trialinfo.cond});
    %-
    fname_ext = cell(1,5);
    for ct = 1:size(nolog_eeg_psd,3)
        %## GET DATA
        ext_tmp = nolog_eeg_psd(f_ind,:,ct);     

        %## REMOVE FOOOF OF EACH STRIDE
        % for i = 1:size(nolog_eeg_psd,2)
        %     % fr = fooof(freqs,mean(squeeze(nolog_eeg_psd(:,:,j)),2),f_range,settings,return_model);
        %     % fooof_diff = 10*(fr.power_spectrum) - 10*(fr.ap_fit);
        %     spec_in = log10(nolog_eeg_psd(tmp_f,i,j));            
        %     fooof_diff = 10*(spec_in') - 10*(fr.ap_fit);
        %     tmp(:,i,j) = fooof_diff';
        %     fooof_freqs = fr.freqs;
        % end

        %## REMOVE MEAN FOOOF OF COMPONENT
        %- un-log
        ext_tmp = 10.^(ext_tmp/10);
        %- fooof
        fr = fooof(freqs_orig(f_ind),mean(squeeze(ext_tmp),2),f_range,settings,return_model);
        spec_in = 10*log10(squeeze(ext_tmp));
        ext_tmp = spec_in - repmat(10*fr.ap_fit',[1,size(ext_tmp,2)]);
        fooof_freqs = fr.freqs;
        if ct == 1
            fname_ext{ff} = 'meandesignb';
            ff = ff + 1;
        end
        %(01/17/2025) JS, This baselining seems most similar to the group
        %average speed results

        %## REMOVE MEAN FOOOF OF EACH CONDITION | TRIAL
        %- un-log
        % ext_tmp = 10.^(ext_tmp/10);
        % for c_i = 1:length(tmp_conds)
        %     %-
        %     inds_cond = cellfun(@(x) strcmp(x,tmp_conds{c_i}),{trialinfo.cond});
        %     %-
        %     % inds = cellfun(@(x) any(strcmp(x,sprintf('%s',tmp_conds{c_i}))),{trialinfo.trial_num_code});
        %     %-
        %     fr = fooof(freqs_orig(f_ind),mean(squeeze(ext_tmp(:,inds_cond)),2),f_range,settings,return_model);
        %     spec_in = 10*log10(squeeze(ext_tmp(:,inds_cond)));
        %     ext_tmp(:,inds_cond) = spec_in - repmat(10*fr.ap_fit',[1,sum(inds_cond)]);
        % end
        % if ct == 1
        %     fname_ext{ff} = 'meancondb';
        %     ff = ff + 1;
        % end
        %(01/17/2025) JS, removing the condition baseline significantly
        %reduces the power.
       
        %## ASSIGN DATA
        tmp(:,:,ct) = ext_tmp;        
    end
    inds_keep = find(~all(squeeze(tmp(:,:,3))==0,1));
    tmp = tmp(:,inds_keep,:);
    trialinfo = trialinfo(inds_keep);

    %{
    j = 1;
    inds = randi([1,size(tmp,2)],50);
    figure;
    hold on;
    plot(freqs(tmp_f),squeeze(eeg_psd(tmp_f,inds,j)),'DisplayName','original avg.');
    plot(freqs(tmp_f),10*log10(mean(squeeze(nolog_eeg_psd(tmp_f,inds,j)),2)),'k-','LineWidth',5);
    ylabel('10*log_{10}(PSD)');
    xlabel('Frequency (Hz)');
    hold off;
    %-
    j = 1;
    figure;
    hold on;
    plot(freqs(tmp_f),squeeze(tmp(:,inds,j)),'DisplayName','aperiodic');
    plot(freqs(tmp_f),mean(squeeze(tmp(:,inds,j)),2),'k-','DisplayName','aperiodic mean','LineWidth',5);
    ylabel('10*log_{10}(PSD)');
    xlabel('Frequency (Hz)');
    ylim([-17.5,16]);
    hold off;
    %-
    figure;
    hold on;
    plot(log10(mean(squeeze(nolog_eeg_psd(:,:,j)),2)))
    plot(log10(mean(squeeze(nolog_eeg_psd(tmp_f,:,j)),2)));
    % plot(mean(squeeze(eeg_psd(:,:,j)),2))
    plot(fr.power_spectrum);
    plot(fr.ap_fit);
    plot(fr.fooofed_spectrum);
    hold off;
    %}    
    %-
    
    fprintf('considering %i trials, spanning conditions: %s\n',length(trialinfo),strjoin(unique({trialinfo.cond}),','));
    eeg_psd = tmp;
    freqs_orig = fooof_freqs;
    %## GET CLUSTER DATA
    fprintf('Getting cluster information...\n');
    tmp_subjs = {tmp_cl_study.datasetinfo.subject};
    comp_arr = zeros(3,length(CLUSTER_PICS)); %(1, comps; 2, old_comps; 3, cluster
    for i = 1:length(CLUSTER_PICS)
        cl_i = CLUSTER_PICS(i);
        ind = find(strcmp(EEG.subject,tmp_subjs));
        ind = tmp_study.cluster(cl_i).sets == ind;
        if any(ind)
            comp_arr(1,i) = tmp_study.cluster(cl_i).comps(ind);
            comp_arr(3,i) = cl_i;
            comp_arr(2,i) = EEG.etc.urreject.ic_keep(comp_arr(1,i));
        end
    end
    comp_arr = comp_arr(:,all(comp_arr,1));

    %## SANITY CHECK
    %- PLOT
    % j = randi([1,size(tmp,3)],1);
    [~,subs] = intersect(comp_arr(3,:),main_cl_inds);
    tmp_arr = comp_arr(:,subs);
    % ci = randi([1,size(tmp_arr,2)],1);
    % ci = 5;
    for ci = 1:size(tmp_arr,2)
        ct = tmp_arr(1,ci);
        fig = figure;
        hold on;
        %- standard error shading
        LINE_ALPHA_PSDFF = 0.8;
        LINE_WIDTH_PSDFF = 2;    
        color_dark = cmaps_speed;
        color_light = cmaps_speed+0.15;    
        for c_i = 1:length(tmp_conds)
            %-
            inds_cond = cellfun(@(x) strcmp(x,tmp_conds{c_i}),{trialinfo.cond});
            data = squeeze(tmp(:,inds_cond,ct));
            %- std error
            % std_error = std(data,[],2); %std(data,[],2)/sqrt(size(data,2));
            % [Pa,Li] = JackKnife_sung(freqs_orig,mean(data,2),mean(data,2)-std_error,mean(data,2)+std_error,...
            %     color_dark(c_i,:),color_light(c_i,:));
            %- spread
            [Pa,Li] = JackKnife_sung(freqs_orig,mean(data,2),prctile(data,25,2),prctile(data,75,2),...
                color_dark(c_i,:),color_light(c_i,:));
            Pa.EdgeColor = "none";
        end
        for c_i = 1:length(tmp_conds)
            %-
            inds_cond = cellfun(@(x) strcmp(x,tmp_conds{c_i}),{trialinfo.cond});
            data = tmp(:,inds_cond,ct)';
            ax = plot(freqs_orig,mean(data),'color',color_dark(c_i,:),...
                'linewidth',LINE_WIDTH_PSDFF,'LineStyle','-','displayname',sprintf('%s',xtick_label_g{c_i}));
            set(ax,'Color',[color_dark(c_i,:),LINE_ALPHA_PSDFF]);
        end
        ylabel('10*log_{10}(PSD)');
        xlabel('Frequency (Hz)');
        legend();
        title(sprintf("CL%i) IC Clean Index: %i, IC Orig. Index: %i",tmp_arr(3,ci),tmp_arr(1,ci),tmp_arr(2,ci)));
        fname_ext = fname_ext(~cellfun(@isempty,fname_ext));
        exportgraphics(fig,[EEG.filepath filesep sprintf('valid_fig_cl%i_ic%i_%s.png',tmp_arr(3,ci),tmp_arr(1,ci),strjoin(fname_ext,'_'))]);
        close(fig);
    end

    %## EXTRACT IMU DATA
    fprintf('Performing gait kinematic calculations...\n');
    tband = (freqs_orig > theta_band_lims(1) & freqs_orig < theta_band_lims(2));
    aband = (freqs_orig > alpha_band_lims(1) & freqs_orig < alpha_band_lims(2));
    bband = (freqs_orig > beta_band_lims(1) & freqs_orig < beta_band_lims(2));
    body_posx = find(contains({EEG.chanlocs.labels},'body_xpos','IgnoreCase',true));
    body_posy = find(contains({EEG.chanlocs.labels},'body_ypos','IgnoreCase',true));
    body_posz = find(contains({EEG.chanlocs.labels},'body_zpos','IgnoreCase',true));
    % world_posx = find(contains({EEG.chanlocs.labels},'world_xpos','IgnoreCase',true));
    % world_posy = find(contains({EEG.chanlocs.labels},'world_ypos','IgnoreCase',true));
    % world_posz = find(contains({EEG.chanlocs.labels},'world_zpos','IgnoreCase',true));
    EEG = eeg_checkset(EEG,'loaddata');
    %##
    steps_struct = repmat(def_step_structs,[1,length(EEG.epoch)*size(comp_arr,2)]);
    cnt = 1;
    %- loop through each condition
    % errors: {'H1018'   }
    % {'H1019'   }
    % {'H1029'   }
    % {'H1036'   }
    % {'H1038'   }
    % {'H1039'   }
    % {'H2018_FU'}
    % {'H2021'   }
    % {'H2038'   }
    % {'H2052'   }
    % {'H3072'   }
    % {'NH3070'  }
    % {'NH3086'  }
    % {'NH3090'  }
    % {'NH3112'  }
    % {'NH3113'  }
    % {'NH3128'  }
    for cond_i = 1:length(conds)
        %##
        % inds = find(cellfun(@(x) any(strcmp(x,conds{cond_i})),{EEG.epoch.eventcond}));
        inds_cond = find(cellfun(@(x) any(strcmp(x,conds{cond_i})),{trialinfo.cond}));
        %- imu data
        datx = squeeze(EEG.data(body_posx,:,inds_cond)); % AP (anteroposterior)
        daty = squeeze(EEG.data(body_posy,:,inds_cond)); % ML (mediolateral)
        datz = squeeze(EEG.data(body_posz,:,inds_cond)); % UD (Up-Down)
        %- imu data
        % datx = squeeze(EEG.data(world_posx,:,inds_cond)); % AP (anteroposterior)
        % daty = squeeze(EEG.data(world_posy,:,inds_cond)); % ML (mediolateral)
        % datz = squeeze(EEG.data(world_posz,:,inds_cond)); % UD (Up-Down)
        
        %- sanity check
        %{
        AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],...
                'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
        fig = figure('color','white');
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        hold on;
        set(gca,AXES_DEFAULT_PROPS{:})
        % title(conds{cond_i})
        ax = axes();
        hold on;
        for s_i = 10:13
            plot3(ax,datx(:,s_i),daty(:,s_i),datz(:,s_i),'LineWidth',2);
        end
        view(ax,[45,25])
        ylabel('Mediolateral (m)');
        xlabel('Anteroposterior (m)');
        zlabel('Up-down (m)');
        set(ax,'Position',[0.25,0.5,0.5,0.5]);  %[left bottom width height]                
        set(ax,'projection', 'perspective', 'box', 'on')
        align_axislabel([], ax)
        hold off;
        exportgraphics(fig,[save_dir filesep sprintf('example_imu_fig.pdf')],'ContentType','vector')
        exportgraphics(fig,[save_dir filesep sprintf('example_imu_fig.tif')],'Resolution',900)
        %}

        %## 2D PLOT
        %{
        AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],...
                'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
        fig = figure('color','white');
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        hold on;
        set(gca,AXES_DEFAULT_PROPS{:})
        title(conds{cond_i})
        %-
        ax = axes();
        hold on;
        for s_i = 10:13
            plot(ax,EEG.times,daty(:,s_i),'LineWidth',2);
        end
        set(ax,'Position',[0.25,0.75,0.5,0.25]);  %[left bottom width height]
        ylabel(ax,'Mediolateral Position (m)');
        xlabel(ax,'time (ms)');
        %-
        hold on;
        ax = axes();
        for s_i = 10:13
            plot(ax,EEG.times,datx(:,s_i),'LineWidth',2);
        end
        set(ax,'Position',[0.25,0.45,0.5,0.25]);  %[left bottom width height]
        ylabel(ax,'Anteroposterior Position (m)');
        xlabel(ax,'time (ms)');
        %-
        hold on;
        ax = axes();
        for s_i = 10:13
            plot(ax,EEG.times,datz(:,s_i),'LineWidth',2);
        end
        set(ax,'Position',[0.25,0.15,0.5,0.25]);  %[left bottom width height]
        ylabel(ax,'Up-down Position (m)');
        xlabel(ax,'time (ms)');
        hold off;
        exportgraphics(fig,[save_dir filesep sprintf('%s_%s_xyz_imu_fig.tif',subj_chars{subj_i},conds{cond_i})],'Resolution',300)
        %}
        %## ADVANCED EPOCH'ING?
        %- before & after inds
        do_log_s = zeros(1,length(inds_cond));
        tmp_step_struct = repmat(def_step_structs,[1,length(inds_cond)]);        
        for s_i = 2:length(inds_cond)-1
            %## DETERMINE LS FOOT EVENTS for the event and surrounding           
            tmp = EEG.epoch(inds_cond(s_i));
            rhs_1 = nan();
            rhs_2 = nan();
            lhs_1 = nan();
            lto_1 = nan();
            rto_1 = nan();
            %## 
            t_rhi = find(strcmp(tmp.eventtype,'RHS'));
            lhi = find(strcmp(tmp.eventtype,'LHS'));
            rti = find(strcmp(tmp.eventtype,'RTO'));
            lti = find(strcmp(tmp.eventtype,'LTO'));
            %##
            good_s = true;
            ii = 0;
            while good_s
                if length(t_rhi) > 2
                    rhi = [t_rhi(1+ii),t_rhi(2+ii)];
                else
                    rhi = t_rhi;
                end
                rht = [tmp.eventlatency{rhi}]*TIME_FACTOR;
                lht = [tmp.eventlatency{lhi}]*TIME_FACTOR;
                rtt = [tmp.eventlatency{rti}]*TIME_FACTOR;
                ltt = [tmp.eventlatency{lti}]*TIME_FACTOR;
                %-
                rhs_1 = rht(1); %tmp.eventlatency{rhi(1)}*TIME_FACTOR;
                rhs_2 = rht(2); %tmp.eventlatency{rhi(2)}*TIME_FACTOR;
                %- index critera
                % ind = lhi > rhi(1) & lhi < rhi(2);
                % lhs_1 = tmp.eventlatency{lhi(ind)}*TIME_FACTOR;
                % ind = lti > rhi(1) & lti < rhi(2);
                % lto_1 = tmp.eventlatency{lti(ind)}*TIME_FACTOR;
                % ind = rti > rhi(1) & rti < rhi(2);
                % rto_1 = tmp.eventlatency{rti(ind)}*TIME_FACTOR;
                %- latency criteria
                ind = lht > rht(1) & lht < rht(2);
                lhs_1 = lht(ind); %tmp.eventlatency{lhi(ind)}*TIME_FACTOR;
                ind = ltt > rht(1) & ltt < rht(2);
                lto_1 = ltt(ind); %tmp.eventlatency{lti(ind)}*TIME_FACTOR;
                ind = rtt > rht(1) & rtt < rht(2);
                rto_1 = rtt(ind); %tmp.eventlatency{rtt(ind)}*TIME_FACTOR;
                if all(~cellfun(@isempty,{rhs_1,lhs_1,lto_1,rto_1,rhs_2})) && ...
                        length([rhs_1,lhs_1,lto_1,rto_1,rhs_2]) == 5
                    %-
                    log_s = true;
                    good_s = false;
                    fprintf('x');
                else
                    % fprintf('Error. Trying next consecutive strikes...\n');
                    fprintf('.');
                    ii = ii + 1;
                    if ii+2 > length(t_rhi)
                        fprintf('o');
                        good_s = false;
                        log_s = false;
                    end
                end
            end
            %## CALCULATE IMU MEASURES
            %- rhs_1, lto_1, lhs_1, rto_1, rhs_2
            EVENT_ERR = 2.001; % potentially a millisecond discrepency for some reason
            if log_s
                %- rhs_1 ind
                diff = EEG.times - rhs_1/TIME_FACTOR;
                rhsi_1 = find(diff > -EVENT_ERR & diff < EVENT_ERR,1,'first');
                %- rhs_1 ind
                diff = EEG.times - lto_1/TIME_FACTOR;
                ltoi_1 = find(diff > -EVENT_ERR & diff < EVENT_ERR,1,'first');
                %- rhs_1 ind
                diff = EEG.times - lhs_1/TIME_FACTOR;
                lhsi_1 = find(diff > -EVENT_ERR & diff < EVENT_ERR,1,'first');
                %- rhs_1 ind
                diff = EEG.times - rto_1/TIME_FACTOR;
                rtoi_1 = find(diff > -EVENT_ERR & diff < EVENT_ERR,1,'first');
                %- rhs_2 ind
                diff = EEG.times - rhs_2/TIME_FACTOR;
                rhsi_2 = find(diff > -EVENT_ERR & diff < EVENT_ERR,1,'first');
                %- calculate pelvis movement
                %* step width (ml exc)
                vec = [daty(rhsi_1,s_i),daty(lhsi_1,s_i)];
                m1 = min(vec);
                m2 = max(vec);
                ml_exc = sqrt((m2-m1)^2);
                %* ap exc
                vec = [datx(rhsi_1,s_i),datx(lhsi_1,s_i)];
                m1 = min(vec);
                m2 = max(vec);
                ap_exc = sqrt((m2-m1)^2);
                %* ud exc
                vec = [datz(rhsi_1,s_i),datz(lhsi_1,s_i)];
                m1 = min(vec);
                m2 = max(vec);
                ud_exc = sqrt((m2-m1)^2);
                %* min max exc (ml)
                vec1 = [daty(rhsi_1:lhsi_1,s_i)];
                vec2 = [daty(lhsi_1:rhsi_2,s_i)];
                m1 = mean([min(vec1),min(vec2)]);
                m2 = mean([max(vec1),max(vec2)]);
                ml_mm_exc = sqrt((m2-m1)^2);
                %* min max exc (ap)
                vec1 = [datx(rhsi_1:lhsi_1,s_i)];
                vec2 = [datx(lhsi_1:rhsi_2,s_i)];
                m1 = mean([min(vec1),min(vec2)]);
                m2 = mean([max(vec1),max(vec2)]);
                ap_mm_exc = sqrt((m2-m1)^2);
                %* min max exc (ud)
                vec1 = [datz(rhsi_1:lhsi_1,s_i)];
                vec2 = [datz(lhsi_1:rhsi_2,s_i)];
                m1 = mean([min(vec1),min(vec2)]);
                m2 = mean([max(vec1),max(vec2)]);
                ud_mm_exc = sqrt((m2-m1)^2);
                %* min max exc gc (ml)
                vec = [daty(rhsi_1:rhsi_2,s_i)];
                m1 = min(vec);
                m2 = max(vec);
                ml_mm_gc_exc = sqrt((m2-m1)^2);
                %* min max exc gc (ap)
                vec = [datx(rhsi_1:rhsi_2,s_i)];
                m1 = min(vec);
                m2 = max(vec);
                ap_mm_gc_exc = sqrt((m2-m1)^2);
                %* min max exc gc (ud)
                vec = [datz(rhsi_1:rhsi_2,s_i)];
                m1 = min(vec);
                m2 = max(vec);
                ud_mm_gc_exc = sqrt((m2-m1)^2);
                %-
                tmp_step_struct(s_i).stride_n = inds_cond(s_i);
                tmp_step_struct(s_i).rhs_1 = rhs_1;
                tmp_step_struct(s_i).rhs_2 = rhs_2;
                tmp_step_struct(s_i).lhs_1 = lhs_1;
                tmp_step_struct(s_i).lto_1 = lto_1;
                tmp_step_struct(s_i).rto_1 = rto_1;
                tmp_step_struct(s_i).gait_cycle_dur = rhs_2 - rhs_1;
                tmp_step_struct(s_i).stance_dur = rto_1 - rhs_1;
                tmp_step_struct(s_i).swing_dur = rhs_2 - rto_1;
                tmp_step_struct(s_i).double_sup_dur = rto_1 - lhs_1;
                tmp_step_struct(s_i).single_sup_dur = lhs_1 - lto_1;
                tmp_step_struct(s_i).step_dur = rhs_2 - lhs_1;
                tmp_step_struct(s_i).ml_exc = ml_exc;
                tmp_step_struct(s_i).ap_exc = ap_exc;
                tmp_step_struct(s_i).ud_exc = ud_exc;
                % fprintf('ml_exc_mm: %0.2g\n',ml_mm_exc);
                tmp_step_struct(s_i).ml_exc_mm = ml_mm_exc;
                tmp_step_struct(s_i).ap_exc_mm = ap_mm_exc;
                tmp_step_struct(s_i).ud_exc_mm = ud_mm_exc;
                %-
                tmp_step_struct(s_i).ml_exc_mm_gc = ml_mm_gc_exc;
                tmp_step_struct(s_i).ap_exc_mm_gc = ap_mm_gc_exc;
                tmp_step_struct(s_i).ud_exc_mm_gc = ud_mm_gc_exc;
                %-
                do_log_s(s_i) = true;
            end            
        end
        fprintf('\n');
        %## ADD EEG MEASURES
        for s_i = 2:length(inds_cond)-1
            %- struct
            if do_log_s(s_i) %&& all(~isempty([rhs_1,rhs_2,lhs_1,lto_1,rto_1]))
                for i = 1:length(CLUSTER_PICS)
                    cl_i = CLUSTER_PICS(i);
                    comp_n = comp_arr(1,comp_arr(3,:) == cl_i);
                    if ~isempty(comp_n)
                        steps_struct(cnt).subj_char = EEG.subject;
                        steps_struct(cnt).cond_char = conds{cond_i}; %speed_n_chars{cond_i}; %conds{cond_i};
                        steps_struct(cnt).group_char = EEG.group;
                        steps_struct(cnt).speed_char = speed_n_chars{cond_i}; %conds{cond_i};
                        steps_struct(cnt).speed_n = double(string(speed_n_chars{cond_i})); %conds{cond_i};
                        steps_struct(cnt).model_char = 'speed';
                        steps_struct(cnt).model_n = des_i;
                        steps_struct(cnt).cluster_n = cl_i;
                        steps_struct(cnt).comp_n = comp_n; % new comp number
                        steps_struct(cnt).comp_n_old = comp_arr(2,comp_arr(3,:) == cl_i); % old comp number
                        steps_struct(cnt).stride_n = tmp_step_struct(s_i).stride_n;
                        steps_struct(cnt).rhs_1 = tmp_step_struct(s_i).rhs_1;
                        steps_struct(cnt).rhs_2 = tmp_step_struct(s_i).rhs_2;
                        steps_struct(cnt).lhs_1 = tmp_step_struct(s_i).lhs_1;
                        steps_struct(cnt).lto_1 = tmp_step_struct(s_i).lto_1;
                        steps_struct(cnt).rto_1 = tmp_step_struct(s_i).rto_1;
                        steps_struct(cnt).gait_cycle_dur = tmp_step_struct(s_i).gait_cycle_dur;
                        steps_struct(cnt).stance_dur = tmp_step_struct(s_i).stance_dur;
                        steps_struct(cnt).swing_dur = tmp_step_struct(s_i).swing_dur;
                        steps_struct(cnt).double_sup_dur = tmp_step_struct(s_i).double_sup_dur;
                        steps_struct(cnt).single_sup_dur = tmp_step_struct(s_i).single_sup_dur;
                        steps_struct(cnt).step_dur = tmp_step_struct(s_i).step_dur;
                        %-
                        steps_struct(cnt).ml_exc = tmp_step_struct(s_i).ml_exc;
                        steps_struct(cnt).ap_exc = tmp_step_struct(s_i).ap_exc;
                        steps_struct(cnt).ud_exc = tmp_step_struct(s_i).ud_exc;
                        % fprintf('ml_exc_mm: %0.2g\n',ml_mm_exc);
                        steps_struct(cnt).ml_exc_mm = tmp_step_struct(s_i).ml_exc_mm;
                        steps_struct(cnt).ap_exc_mm = tmp_step_struct(s_i).ap_exc_mm;
                        steps_struct(cnt).ud_exc_mm = tmp_step_struct(s_i).ud_exc_mm;
                        %-
                        steps_struct(cnt).ml_exc_mm_gc = tmp_step_struct(s_i).ml_exc_mm_gc;
                        steps_struct(cnt).ap_exc_mm_gc = tmp_step_struct(s_i).ap_exc_mm_gc;
                        steps_struct(cnt).ud_exc_mm_gc = tmp_step_struct(s_i).ud_exc_mm_gc;
                        % step_structs(cnt).eeg_psd = squeeze(eeg_psd(:,cnt,comp_n));
                        % step_structs(cnt).freqs = freqs;
                        steps_struct(cnt).avg_theta_pre = mean(eeg_psd(tband,inds_cond(s_i-1),comp_n));
                        steps_struct(cnt).avg_alpha_pre = mean(eeg_psd(aband,inds_cond(s_i-1),comp_n));
                        steps_struct(cnt).avg_beta_pre = mean(eeg_psd(bband,inds_cond(s_i-1),comp_n));
                        steps_struct(cnt).avg_theta = mean(eeg_psd(tband,inds_cond(s_i),comp_n));
                        steps_struct(cnt).avg_alpha = mean(eeg_psd(aband,inds_cond(s_i),comp_n));
                        steps_struct(cnt).avg_beta = mean(eeg_psd(bband,inds_cond(s_i),comp_n));
                        steps_struct(cnt).avg_theta_post = mean(eeg_psd(tband,inds_cond(s_i+1),comp_n));
                        steps_struct(cnt).avg_alpha_post = mean(eeg_psd(aband,inds_cond(s_i+1),comp_n));
                        steps_struct(cnt).avg_beta_post = mean(eeg_psd(bband,inds_cond(s_i+1),comp_n));
                        % step_structs(cnt).eeg_comp_arr = comp_arr;
                        % step_structs(cnt).eeg_comp_arr_deats = ['row 1 is the component in the reduced eeg set;' ...
                        %     ' row 2 is the component in the old eeg set; row 3 is the cluster number'];
                        cnt = cnt + 1;
                    end
                end
            end
        end
        steps_struct = steps_struct(~cellfun(@isempty,{steps_struct.avg_theta_post}));
        steps_struct = steps_struct(~cellfun(@isempty,{steps_struct.ml_exc}));
        steps_struct = steps_struct(~cellfun(@isempty,{steps_struct.comp_n}));
        %- grab condition     
        tmp_ss = steps_struct(inds_cond);

        % %- per condition measures (step dur)
        % musd = mean([tmp_ss.step_dur]);
        % stdsd = std([tmp_ss.step_dur]);
        % tmp = [tmp_ss.step_dur];
        % %-
        % vals = sqrt((tmp-musd).^2);
        % vals = num2cell(vals);
        % [steps_struct(inds_cond).var_step_dur_1] = deal(vals{:});
        % %-
        % vals = stdsd./sqrt((tmp-musd).^2);
        % vals = num2cell(vals);
        % [steps_struct(inds_cond).var_step_dur_2] = deal(vals{:});
        % %-
        % vals =sqrt((tmp-musd).^2)./stdsd;
        % vals = num2cell(vals);
        % [steps_struct(inds_cond).var_step_dur_3] = deal(vals{:});
        % 
        % %- per condition measures (step width)
        % musd = mean([tmp_ss.ml_exc_mm]);
        % stdsd = std([tmp_ss.ml_exc_mm]);
        % tmp = [tmp_ss.ml_exc_mm];
        % %-
        % vals = sqrt((tmp-musd).^2);
        % vals = num2cell(vals);
        % [steps_struct(inds_cond).var_ml_exc_1] = deal(vals{:});
        % %-
        % vals = stdsd./sqrt((tmp-musd).^2);
        % vals = num2cell(vals);
        % [steps_struct(inds_cond).var_ml_exc_2] = deal(vals{:});
        % fprintf('\n');
    end
    %- assign to struct
    % vals = num2cell(vals)';
    % [steps_struct(inds_cond)] = deal(vals{:});
    %- per design measures
    %## SAVE
    steps_struct = steps_struct(~cellfun(@isempty,{steps_struct.avg_theta_post}));
    steps_struct = steps_struct(~cellfun(@isempty,{steps_struct.ml_exc}));
    steps_struct = steps_struct(~cellfun(@isempty,{steps_struct.comp_n}));
    steps_struct = struct2table(steps_struct);
    fname_ext = fname_ext(~cellfun(@isempty,fname_ext));
    par_save(steps_struct,EEG.filepath,sprintf('kin_eeg_psd_meas_%s.mat',strjoin(fname_ext,'_')));
    fname_ext = fname_ext(~cellfun(@isempty,fname_ext));
    fname_ext_store = strjoin(fname_ext,'_');
    %- turn off when doing parallel processing.
    EEG.data = 'eeg_imu_epochs.fdt';
    tmp_alleeg{subj_i} = EEG;
end
%% (RESAVE STUDY) ====================================================== %%
% %- remove bugged out subjects
% tmp_alleeg = tmp_alleeg(~cellfun(@isempty,tmp_alleeg));
% %- clear data for memory
% % for i = 1:length(tmp_alleeg)
% %     tmp_alleeg{i}.data = 'eeg_imu_epochs.fdt';
% % end
% %## BOOKKEEPING (i.e., ADD fields not similar across EEG structures)
% tmp_alleeg = util_resolve_struct(tmp_alleeg);
% %##
% [STUDY, ALLEEG] = std_editset([],tmp_alleeg,...
%                                 'updatedat','off',...
%                                 'savedat','off',...
%                                 'name','d_eeg_psd_imuls_sts_anl',...
%                                 'filename','d_eeg_psd_imuls_sts_anl',...
%                                 'filepath',save_dir);
% [STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
% [STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
%                                         STUDY.filename,STUDY.filepath,...
%                                         'RESAVE_DATASETS','off');
%% (TEST CHANGES ACROSS SPEEDS) ======================================== %%
data_store = table.empty;
groups = {'H1000','H2000','H3000'};
group_name = {'younger_adults','older_high_function','older_low_function'};
for i = 1:length(STUDY.datasetinfo)
    try
        % tmp = par_load(STUDY.datasetinfo(i).filepath,'ls_eeg_psd_struct.mat');
        % tmp = par_load(STUDY.datasetinfo(i).filepath,'ls_eeg_psd_struct_prepost_gait_imufix_condb.mat');
        % tmp = par_load(fileparts(STUDY.datasetinfo(i).filepath),'kin_eeg_anl.mat'); %(01/16/2025) JS, silly bug from changing par_save need to fix in future
        tmp = par_load(STUDY.datasetinfo(i).filepath,sprintf('kin_eeg_psd_meas_%s.mat',fname_ext_store));
        vals = regexp(tmp.subj_char{1},'[nN]?[hH](\d)\d*','tokens');
        vals = double(string(vals{1}{1}));
        tt = repmat(vals,[1,size(tmp,1)]);
        tmp.group_n = tt';
        tt = repmat(groups(vals),[1,size(tmp,1)]);
        [tmp.group_char] = tt'; %deal(tt{:});
        tt = repmat(group_name(vals),[1,size(tmp,1)]);
        [tmp.group_name] = tt'; %deal(tt{:});
        data_store = [data_store; tmp];
    catch e
        fprintf('Couldn''t load %s...\n%s\n',STUDY.datasetinfo(i).subject,getReport(e));
    end
    % data_store = cat(1,data_store,tmp);
end
% steps_struct = data_store;
% data_store = struct.empty;
% writetable(data_store,[save_dir filesep 'step_by_step_eeg_psd_table.xlsx']);
% par_save(data_store,save_dir,'step_by_step_eeg_psd_table.mat');
writetable(data_store,[save_dir filesep sprintf('sbs_eeg_psd_%s.xlsx',fname_ext_store)]);
par_save(data_store,save_dir,sprintf('sbs_eeg_psd_%s.mat',fname_ext_store));
% writetable(data_store,[save_dir filesep 'step_by_step_eeg_table_indvfooof.xlsx']);
% par_save(data_store,save_dir,'step_by_step_eeg_table_indvfooof.mat');
%% (TEST DATA) ========================================================= %%
%{
%## PARAMSd
% kin_char = 'step_dur';
kin_char = 'var_step_dur_1'; xlab = 'SD_{step}-SD_{mu} (s)';
% kin_char = 'var_step_dur_2'; xlab = 'SD_{std}/(sqrt(SD_{step}-SD_{mu})^2)';
% kin_char = 'var_step_dur_3'; xlab = '(sqrt(SD_{step}-SD_{mu})^2)/SD_{std}';
eeg_char = 'avg_theta';
i = 1;
cl_i = CLUSTER_PICS(i);

%## PICK CLUSTER & MODEL
disp(unique(data_store.cluster_n));
inds = strcmp(data_store.model_char,'speed') & data_store.cluster_n == cl_i;
tmp_tbl = data_store(inds,:);
subj_chars = unique(tmp_tbl.subj_char);
cond_chars = unique(tmp_tbl.c ond_char);
fprintf('CL%i) Subjects: %s\n',cl_i,sprintf('%s,',subj_chars{:}));

%## MODEL DATA
%- (1) LME w/ random intercept for subject : step_dur ~ speed + (1|subj)
mod_out = sprintf('%s ~ %s*speed_n',eeg_char,kin_char);
stats_out = fitlme(tmp_tbl,mod_out);
anv_out = anova(stats_out);
%- test normality
[norm_h,norm_p] = lillietest(stats_out.Residuals.Raw);
%- intercept only model
altmod_out = sprintf('%s ~ 1',eeg_char);
altstats_out = fitlm(tmp_tbl,altmod_out);
%- alternative f2?
R21 = altstats_out.SSR/altstats_out.SST; % coefficient of determination
R22 = stats_out.SSR/stats_out.SST; % coefficient of determination
alt_f2 = (R22-R21)/(1-R22);

%## VISUALIZE DATA
COLOR_MAPS_SPEED = linspecer(4*3);
COLOR_MAPS_SPEED = [COLOR_MAPS_SPEED(1,:);COLOR_MAPS_SPEED(2,:);COLOR_MAPS_SPEED(3,:);COLOR_MAPS_SPEED(4,:)];
color_conds = COLOR_MAPS_SPEED;


%## VISUALIZE DATA
%- condition-wise coloring
figure;
title({'eeg theta variations with step duration','colored by condition'});
hold on;
for cond_i = 1:length(cond_chars)
    %-
    tmp_dat = strcmp(tmp_tbl.cond_char,cond_chars{cond_i});
    tmp_dat = tmp_tbl(tmp_dat,:);
    ss = scatter(tmp_dat,kin_char,eeg_char,'DisplayName',sprintf('%s',cond_chars{cond_i}));
    ss.CData = color_conds(cond_i,:);
    ss.SizeData = 15;
    ss.MarkerEdgeAlpha = 0.6;
    %-
end
xlabel(xlab)
ylabel('10*log_{10}(PSD-Ap. Fit)')


%- subject-wise coloring
figure;
title({'eeg theta variations','colored by subject'});
hold on;
color_subjs = linspecer(length(subj_chars));
for subj_i = 1:length(subj_chars)
    tmp_dat = strcmp(tmp_tbl.subj_char,subj_chars{subj_i});
    tmp_dat = tmp_tbl(tmp_dat,:);
    ss = scatter(tmp_dat,kin_char,eeg_char,'DisplayName',sprintf('%s',subj_chars{subj_i}));
    ss.CData = color_subjs(subj_i,:);
    ss.SizeData = 15;
    % ss.MarkerFaceAlpha = 'flat';
    % ss.AlphaData = repmat(0.3,[size(data,1),1]);
    ss.MarkerEdgeAlpha = 0.6;
end
xlabel(xlab)
ylabel('10*log_{10}(PSD-Ap. Fit)')

%- violin plot
figure;
title('In Each Step: Step Duration Changes w/ Walking Speed')
tmp_dat = tmp_tbl.(kin_char);
tmp_cats = tmp_tbl.cond_char;
violinplot(tmp_dat,tmp_cats);
xlabel('Walking Speed (m/s)')
ylabel(xlab)

%- violin plot
figure;
title('In Each Step: Average Theta Changes w/ Walking Speed')
tmp_dat = tmp_tbl.avg_theta;
tmp_cats = tmp_tbl.cond_char;
violinplot(tmp_dat,tmp_cats);
xlabel('Speed (m/s)')
ylabel('10*log_{10}(PSD-Ap. Fit)')
%}


