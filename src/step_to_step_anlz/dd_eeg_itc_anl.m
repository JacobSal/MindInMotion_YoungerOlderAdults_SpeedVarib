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
%% (SAVE TABLE) ======================================================== %%
%{
% fext = 'phasec_notw_mw';
% fext = 'phasec_notw_mw_based';

fext = 'phasec_notw_mw_based';
itc_so = cell(length(CL_STUDY.datasetinfo),1);
for subj_i = 1:length(CL_STUDY.datasetinfo)
    try
        tmp_sbs_study = SBS_STUDY;
        tmp_cl_study = CL_STUDY;
        subj_char = tmp_sbs_study.datasetinfo(subj_i).subject;
        subj_ind = strcmp(tmp_cl_study.datasetinfo(subj_i).subject,{tmp_sbs_study.datasetinfo.subject});
        %--
        if ~isempty(subj_ind) && any(subj_ind)
            itc_dato = par_load([tmp_cl_study.datasetinfo(subj_i).filepath filesep sprintf('%s_custom_itc_%s.mat',tmp_cl_study.datasetinfo(subj_i).subject,fext)]);
            %--
            fprintf('%s) Adding ITC structure\n',tmp_cl_study.datasetinfo(subj_i).subject);
            tmp = itc_dato.itc_dat_struct;
            tmp = rmfield(tmp,{'itc_freqs','itc_times'});
            tmp = struct2table(tmp);
            itc_so{subj_i} = tmp;
        end
    catch e 
        fprintf('%s) %s\n',tmp_cl_study.datasetinfo(subj_i).subject ,getReport(e));
    end
end
%--
itc_so = util_resolve_table(itc_so);
par_save(itc_so,save_dir,sprintf('itc_table_%s_new.mat',fext));
%(04/30/2025) JS, NH3128 gets an OOM error on the fext='phasec_notw_mw'
%run. Need to include them for the final analysis
%}
%% (SAVE R TABLE) =================================================== %%

% fext = 'phasec_notw_mw';
% fext = 'phasec_notw_mw_based';
fext = 'phasec_notw_mw_based';
itc_so = cell(length(CL_STUDY.datasetinfo),1);
CL_NUM_CUTOFF = 13;
GROUP_CHARS = {'H1000','H2000','H3000'};
for subj_i = 1:length(CL_STUDY.datasetinfo)
    try
        tmp_sbs_study = SBS_STUDY;
        tmp_cl_study = CL_STUDY;
        subj_char = tmp_sbs_study.datasetinfo(subj_i).subject;
        subj_ind = strcmp(tmp_cl_study.datasetinfo(subj_i).subject,{tmp_sbs_study.datasetinfo.subject});
        %--
        if ~isempty(subj_ind) && any(subj_ind)
            itc_dato = par_load([tmp_cl_study.datasetinfo(subj_i).filepath filesep sprintf('%s_custom_itc_%s.mat',tmp_cl_study.datasetinfo(subj_i).subject,fext)]);
            %--
            fprintf('%s) Adding ITC structure\n',tmp_cl_study.datasetinfo(subj_i).subject);
            tmp = itc_dato.itc_dat_struct;
            twp = itc_dato.parameters;
            %--
            FREQ_BOUND = [3,60];
            TIME_BOUND = [twp.timewarpms(1),twp.timewarpms(end)];
            % FREQ_BOUND = [3,60];
            % TIME_BOUND = [-500,3000];
            alltimes = tmp(1).itc_times';
            allfreqs = tmp(1).itc_freqs';
            tinds = alltimes > TIME_BOUND(1) & alltimes < TIME_BOUND(2);
            finds = allfreqs > FREQ_BOUND(1) & allfreqs < FREQ_BOUND(2);
            %--
            tclu = unique([tmp.cluster_n]);
            tcu = unique({tmp.cond_char});
            dats_store = cell(length(tclu)*length(tcu),1);
            cnt = 1;
            tti = tic();
            for cl_i = 1:length(tclu)
                if tclu(cl_i) <= CL_NUM_CUTOFF
                    for c_i = 1:length(tcu)
                        inds = [tmp.cluster_n] == tclu(cl_i) & ...
                            strcmp({tmp.cond_char},tcu{c_i});
                        tt = tmp(inds);
                        sz = size(tt.itc_dat(finds,tinds));
                        rtf = repmat(tt.itc_freqs(finds),[1,sz(2)]);
                        rtt = repmat(tt.itc_times(tinds)',[sz(1),1]);
                        %-- convert & take abs of data
                        % rtd = cellfun(@(x) x(finds,tinds),{tt.itc_dat});
                        % rtt = cellfun(@(x) x(finds,tinds),{tt.itc_dat});

                        rtd = num2cell(abs(reshape(tt.itc_dat(finds,tinds),[1,sz(1)*sz(2)])));
                        rtt = num2cell(reshape(rtt,[1,sz(1)*sz(2)]));
                        rtf = num2cell(reshape(rtf,[1,sz(1)*sz(2)]));
                        % tmp_dats = struct(...
                        %     'subj_char',tt.subj_char, ...
                        %     'group_char',tt.group_char, ...
                        %     'cond_char',tt.cond_char, ...
                        %     'mod_n',tt.mod_n, ...
                        %     'cluster_n',tt.cluster_n, ...
                        %     'itc_dat',0, ...
                        %     'itc_freq',0, ...
                        %     'itc_time',0 ...
                        %     );
                        gg = strcmp(tt.group_char,GROUP_CHARS);
                        cc = strcmp(tt.cond_char,tcu);
                        tmp_dats = struct(...
                            'subj_n',subj_i, ...
                            'group_n',find(gg), ...
                            'cond_n',find(cc), ...
                            'mod_n',tt.mod_n, ...
                            'cluster_n',tt.cluster_n, ...
                            'itc_dat',0, ...
                            'itc_freq',0, ...
                            'itc_time',0 ...
                            );
                        tmp_dats = repmat(tmp_dats,[1,length(rtd)]);
                        [tmp_dats(:).itc_dat] = deal(rtd{:});
                        [tmp_dats(:).itc_freq] = deal(rtf{:});
                        [tmp_dats(:).itc_time] = deal(rtt{:});
                        dats_store{cnt} = tmp_dats;
                        cnt = cnt + 1;
                        %--
                        % figure;
                        % contourf(tt.itc_times,tt.itc_freqs,real(tt.itc_dat),30,'lines','none');
                        % set(gca,'fontsize',14);
                        % ylabel('Frequency (Hz)');
                        % xlabel('Time');
                        % colorbar;
                    end
                end
            end
            tmp = util_resolve_struct(dats_store);
            tmp = struct2table(tmp);
            itc_so{subj_i} = tmp;
            fprintf('%s) loading and conversion done. %0.2f min',tmp_cl_study.datasetinfo(subj_i).subject,toc(tti)/60)
        end
    catch e 
        fprintf('%s) %s\n',tmp_cl_study.datasetinfo(subj_i).subject ,getReport(e));
    end
end
%--
itc_so = util_resolve_table(itc_so);
par_save(itc_so,save_dir,sprintf('itc_rdata_table_%s.mat',fext));
%(04/30/2025) JS, NH3128 gets an OOM error on the fext='phasec_notw_mw'
%run. Need to include them for the final analysis

%## EXPORT DATA TO CSV
itc_so = par_load(save_dir,sprintf('itc_rdata_table_%s.mat',fext));
clu = unique(itc_so.cluster_n);
for cl_i = 1:length(clu)
    inds = itc_so.cluster_n == clu(cl_i);
    tmp = itc_so(inds,:);
    writetable(tmp,[save_dir filesep sprintf('%i_itc_rdata_table_%s.csv',clu(cl_i),fext)]);
end
writetable(itc_so,[save_dir filesep sprintf('itc_rdata_table_%s.csv',fext)]);

%% ================================================================== %%
%## VALIDATION PLOT
%{
tmp_plot_struct = struct( ...
    'alltitles',{''},...
    'xlabel','Gait Cycle Time (ms)',...
    'ylabel','Frequency (Hz)',...
    'xticklabel_times',tmp_ntf_struct.timewarpms,...
    'xticklabel_chars',{{'RHS','RTO','LHS','LTO','RHS'}},...
    'xticklabel_angle',45,...
    'clim',double([min(real(itc_dat),[],'all'),max(real(itc_dat),[],'all')]),...
    'font_size',8,...
    'font_name','Arial',...
    'freq_lims',[3,60],...
    'time_lims',[tmp_ntf_struct.timewarpms(1),tmp_ntf_struct.timewarpms(end)],...    
    'contourf_grain',ceil((500/pi())),...
    'alpha_multiple',0.6,...
    'position',[0.1,0.1,0.5,0.5], ...
    'colorbar_shift',[1,1], ...
    'do_display_bandmarks',true);
%--
allersp_mask = ones(size(itc_dat));
allersp_pcond = zeros(size(itc_dat));
% allersp = real(itc_dat); 
allersp = abs(itc_dat);
alltimes = ttimes;
allfreqs = lfreqs;
tmp_plot_struct.clim = double([min(allersp,[],'all'),max(allersp,[],'all')]);
figure;
ax = axes();
[ax] = plot_contourf_cell(ax,allersp,alltimes,allfreqs,...
    allersp_mask,allersp_pcond, ...
    'PLOT_STRUCT',tmp_plot_struct)
%}
%% (GROUP-CONDITION PLOTS) ========================================== %%
%-- load dat table
% fext = 'phasec_notw_mw';
fext = 'phasec_notw_mw_based';
itc_so = par_load([save_dir filesep sprintf('itc_table_%s_new.mat',fext)]);
%-- tmp load
itc_dato = par_load([CL_STUDY.datasetinfo(1).filepath filesep sprintf('%s_custom_itc_%s.mat',CL_STUDY.datasetinfo(1).subject,fext)]);
twp = itc_dato.parameters;
tws = itc_dato.itc_dat_struct;
%% PARAMS =========================================================== %%
%--
TW_CHARS = {'RHS','RTO','LHS','LTO','RHS'};
ALPHA = 0.05;
FREQ_BOUND = [3,60];
TIME_BOUND = [twp.timewarpms(1),twp.timewarpms(end)];
% FREQ_BOUND = [3,60];
% TIME_BOUND = [-500,3000];
alltimes = tws(1).itc_times';
allfreqs = tws(1).itc_freqs';
tinds = alltimes > TIME_BOUND(1) & alltimes < TIME_BOUND(2);
finds = allfreqs > FREQ_BOUND(1) & allfreqs < FREQ_BOUND(2);
tcrop = alltimes(tinds);
fcrop = allfreqs(finds);
%## PARAMETERS
%##
cluster_inds_plot = [3,4,6,8];
clusters = unique(itc_so.cluster_n);
subjects = unique(itc_so.subj_char);
groups = unique(itc_so.group_char);
conditions = {'0p25','0p5','0p75','1p0'};
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
PLOT_STRUCT = struct( ...
    'alltitles',{''},...
    'xlabel','Gait Cycle Time (ms)',...
    'ylabel','Frequency (Hz)',...
    'xticklabel_times',twp.timewarpms, ...
    'xticklabel_chars',{TW_CHARS},...
    'xticklabel_angle',45,...
    'clim',[],...
    'font_size',8,...
    'font_name','Arial',...
    'freq_lims',[3,60],...
    'time_lims',[twp.timewarpms(1),twp.timewarpms(end)],...    
    'contourf_grain',ceil((500/pi())),...
    'alpha_multiple',0.9,...
    'position',[0.1,0.1,0.5,0.5], ...
    'do_display_bandmarks',true ...
    );
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
FIGURE_POSITION =[1,1,6.5,9];
AX_DIM = [0.13,0.16];
AX_SHIFT = [1.2,-1.25];
AX_INIT_X = 0.2;
AX_INIT_Y = 0.775;
X_DIM = 4;

%##
tmp_save_dir = [save_dir filesep 'itc_plots'];
mkdir(tmp_save_dir);

%% (MORE EFFICIENT PLOTTING) ============================================
%## LOOP
for cl_i = 1:length(cluster_inds_plot)
    %##
    %%
    blext = 'none';
    stat_ext = 'none';    
    %--
    cl_ii = find(cluster_inds_plot(cl_i) == double(string(clusters)));
    cl_n = clusters(cl_i);
    atlas_name = cluster_titles{cl_n};

    %## EXTRACT DATA
    allitc = cell(length(conditions),length(groups));
    cropitc = cell(length(conditions),length(groups));
    %--
    for g_i = 1:length(groups)
        for c_i = 1:length(conditions)
            inds = itc_so.cluster_n == cl_n & ...
                strcmp(itc_so.group_char,groups{g_i}) & ...
                strcmp(itc_so.cond_char,conditions{c_i});
            tt = itc_so(inds,:);
            % itc_dat = cellfun(@real,tt.itc_dat,'UniformOutput',false);
            itc_dat = cellfun(@abs,tt.itc_dat,'UniformOutput',false);
            itc_dat = cat(3,itc_dat{:});
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
    %--
    % itco_c = cellfun(@(x) std(x,[],3)/size(x,3),itco_c,'UniformOutput',false);

    %## TTEST MASKING
    % bs_ersp = cell(size(itco_c,1),size(itco_c,2));
    % bs_masked = cell(size(itco_c,1),size(itco_c,2));
    % bs_pval = cell(size(itco_c,1),size(itco_c,2));
    % %--
    % TTEST_STRUCT = struct(...
    %     'tail','both', ...
    %     'alpha',0.05, ...
    %     'cluster_thresh',300);
    % for g_i = 1:size(itco_c,2)
    %     for c_i = 1:size(itco_c,1)
    %         [tfmu,tfmask,tfpv] = eeglab_time_boot(itco_c{c_i,g_i}, ...
    %             'BOOT_STRUCT',BOOT_STRUCT);
    %         comp_mu = mean(tfmask(tfmask>0));
    %         %-- store
    %         tmp_tts = TTEST_STRUCT;
    %         tmp_tts.tail = 'left';
    %         [tfmu,tfmask,tfpv] = eeglab_tf_ttest(itco_c{c_i,g_i}, comp_mu,...
    %             'TTEST_STRUCT',tmp_tts);
    %         %-- store
    %         bs_ersp{c_i,g_i} = tfmu;
    %         bs_masked{c_i,g_i} = tfmask;
    %         bs_pval{c_i,g_i} = tfpv;
    %     end            
    % end

    %## BOOT STRAP MASKING
    % stats = CL_STUDY.etc.statistics;
    % stats.condstats = 'on';
    % stats.groupstats = 'off';
    % stats.paired = {'on','off'};
    %--
    % stats = CL_STUDY.etc.statistics;
    % stats.condstats = 'off';
    % stats.groupstats = 'on';
    % stats.paired = {'on','off'};
    %--
    % [pcond, pgroup, pinter,scond,sgroup,sinter] = std_stat(itco_c, stats);
    % stat_ext = [stats.condstats,stats.groupstats];
    %-- clims
    clim = cellfun(@(x) [prctile(x,3,'all'),prctile(x,97,'all')],itco_c, ...
        'UniformOutput',false);
    clim = mean(cat(1,clim{:}),1);
    % clim = double([-max(abs(clim)),max(abs(clim))]);
    clim = [0,0.20];

    %## INITIATE FIGURE
    x_shift = AX_INIT_X;
    y_shift = AX_INIT_Y;
    x_cnt = 1;
    y_cnt = 1;
    %-- fig
    fig = figure('color','white');
    set(fig,'Units','inches', ...
        'Position',FIGURE_POSITION, ...
        'PaperUnits','inches', ...
        'PaperSize',[1 1], ...
        'PaperPosition',[0 0 1 1])
    p_sz = get(fig,'Position');
    set(gca,AXES_DEFAULT_PROPS{:});
    %-- axes
    hold on;
    %--
    % rsubj = randi(size(itco_c{1,1},3),1);
    for g_i = 1:length(groups)
        for c_i = 1:length(conditions)            
            % tmpdat = itco_c{c_i,g_i}(:,:,rsubj);
            %--
            % itco_c{c_i,g_i} = 1 - itco_c{c_i,g_i};
            % allersp = mean(itco_c{c_i,g_i},3);
            % allersp_mask = mean(itco_c{c_i,g_i},3);
            % allersp = median(itco_c{c_i,g_i},3);
            % allersp_mask = median(itco_c{c_i,g_i},3);
            % allersp_pcond = zeros(size(allersp));
            %--
            allersp = bs_ersp{c_i,g_i};
            allersp_mask = bs_masked{c_i,g_i};
            allersp_pcond = bs_pval{c_i,g_i};
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
            tmp_plot_struct.alltitles = alltitles{c_i};
            tmp_plot_struct.position = [x_shift,y_shift,AX_DIM(1),AX_DIM(2)];
            tmp_plot_struct.clim = clim; %[0,0.175]; %double(clim); %double([min(allersp,[],'all'),max(allersp,[],'all')]);
            if c_i > 1
                tmp_plot_struct.do_display_bandmarks = false;
                tmp_plot_struct.ylabel = '';
            end
            if c_i < length(conditions)
                tmp_plot_struct.do_display_colorbar = false;
            end
            if g_i < length(groups)
                tmp_plot_struct.xticklabel_chars = {''};
                tmp_plot_struct.xlabel = '';
            elseif c_i > 1
                tmp_plot_struct.xticklabel_chars = {''};
                tmp_plot_struct.xlabel = '';
            end
            %--
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
    exportgraphics(fig,[tmp_save_dir filesep sprintf('cl%i_stdstats_%s_%s_%s.tiff',cl_n,fext,stat_ext,blext)],'Resolution',300);
    % close(fig);
end

%% MOD INDEX =========================================================== %%
%--
ALPHA = 0.05;
mi_pfreq_vec = itc_dato.itc_dat_struct.mi_phase;
mi_afreq_vec = itc_dato.itc_dat_struct.mi_freqs;
%## PARAMETERS
%##
cluster_inds_plot = [3,4,6,8];
clusters = unique(itc_so.cluster_n);
subjects = unique(itc_so.subj_char);
groups = unique(itc_so.group_char);
conditions = {'0p25','0p5','0p75','1p0'};
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
PLOT_STRUCT = struct( ...
    'alltitles',{''},...
    'xlabel','Phase Frequency (Hz)',...
    'ylabel','Amplitude Frequency (Hz)',...
    'xticklabel_times',(min(mi_pfreq_vec):10:max(mi_pfreq_vec)), ...
    'xticklabel_chars',{cellstr(string((min(mi_pfreq_vec):5:max(mi_pfreq_vec))))},...
    'xticklabel_angle',45,...
    'clim',[],...
    'font_size',8,...
    'font_name','Arial',...
    'freq_lims',[min(mi_afreq_vec),max(mi_afreq_vec)],...
    'time_lims',[min(mi_pfreq_vec),max(mi_pfreq_vec)],...    
    'contourf_grain',ceil((500/pi())),...
    'alpha_multiple',0.9,...
    'position',[0.1,0.1,0.5,0.5], ...
    'do_display_bandmarks',false ...
    );
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
FIGURE_POSITION =[1,1,6.5,9];
AX_DIM = [0.13,0.16];
AX_SHIFT = [1.2,-1.25];
AX_INIT_X = 0.2;
AX_INIT_Y = 0.775;
X_DIM = 4;

%##
tmp_save_dir = [save_dir filesep 'itc_plots'];
mkdir(tmp_save_dir);

%% (MORE EFFICIENT PLOTTING) ============================================
%## LOOP
for cl_i = 1:length(cluster_inds_plot)
    %##
    %%
    blext = 'none';
    stat_ext = 'none';    
    %--
    cl_ii = find(cluster_inds_plot(cl_i) == double(string(clusters)));
    cl_n = clusters(cl_i);
    atlas_name = cluster_titles{cl_n};

    %## EXTRACT DATA
    allitc = cell(length(conditions),length(groups));
    cropitc = cell(length(conditions),length(groups));
    %--
    for g_i = 1:length(groups)
        for c_i = 1:length(conditions)
            inds = itc_so.cluster_n == cl_n & ...
                strcmp(itc_so.group_char,groups{g_i}) & ...
                strcmp(itc_so.cond_char,conditions{c_i});
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
    % stats = CL_STUDY.etc.statistics;
    % stats.condstats = 'on';
    % stats.groupstats = 'off';
    % stats.paired = {'on','off'};
    %--
    stats = CL_STUDY.etc.statistics;
    stats.condstats = 'off';
    stats.groupstats = 'on';
    stats.paired = {'on','off'};
    %--
    [pcond,pgroup,pinter,scond,sgroup,sinter] = std_stat(itco_c, stats);
    stat_ext = [stats.condstats,stats.groupstats];
    %##
    %-- clims
    clim = cellfun(@(x) [prctile(x,3,'all'),prctile(x,97,'all')],itco_c, ...
        'UniformOutput',false);
    clim = mean(cat(1,clim{:}),1);
    % clim = double([-max(abs(clim)),max(abs(clim))]);
    clim = [0,5e-5];

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
    %--
    % rsubj = randi(size(itco_c{1,1},3),1);
    for g_i = 1:length(groups)
        for c_i = 1:length(conditions)            
            % tmpdat = itco_c{c_i,g_i}(:,:,rsubj);
            %--
            % allersp = mean(itco_c{c_i,g_i},3);
            % allersp_mask = mean(itco_c{c_i,g_i},3);
            % allersp_pcond = zeros(size(allersp));
            %--
            if ~isempty(pcond)
                tmpp = double(pcond{g_i});
                tmps = double(scond{g_i});
            elseif ~isempty(pgroup)
                tmpp = double(pgroup{c_i});
                tmps = double(sgroup{c_i});
            end
            allersp = mean(itco_c{c_i,g_i},3);
            allersp_mask = allersp.*tmpp;
            allersp_pcond = tmps;
            
            %## PLOT
            tmp_plot_struct = PLOT_STRUCT;
            tmp_plot_struct.alltitles = alltitles{c_i};
            tmp_plot_struct.position = [x_shift,y_shift,AX_DIM(1),AX_DIM(2)];
            tmp_plot_struct.clim = clim; %[0,0.175]; %double(clim); %double([min(allersp,[],'all'),max(allersp,[],'all')]);
            if c_i > 1
                tmp_plot_struct.do_display_bandmarks = false;
                tmp_plot_struct.ylabel = '';
            end
            if c_i < length(conditions)
                tmp_plot_struct.do_display_colorbar = false;
            end
            if g_i < length(groups)
                tmp_plot_struct.xticklabel_chars = {''};
                tmp_plot_struct.xlabel = '';
            elseif c_i > 1
                tmp_plot_struct.xticklabel_chars = {''};
                tmp_plot_struct.xlabel = '';
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
                x_shift = AX_INIT_X;
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


