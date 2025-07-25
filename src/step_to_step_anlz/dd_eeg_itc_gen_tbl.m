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
%## R-STATS LOADING
r_stats_dir = [PATHS.src_dir filesep 'r_scripts' filesep 'sbs_lme_mods'];
%--
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
%-
cl_struct = par_load(cluster_k_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
[comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
CLUSTER_PICS = main_cl_inds;
%% (SAVE TABLE) ======================================================== %%
% fext = 'phasec_notw_mw';
% fext = 'phasec_notw_mw_based';
fext = 'phasec_notw_mw_based';
%## LOOP
itc_so = cell(length(CL_STUDY.datasetinfo),1);
parfor subj_i = 1:length(CL_STUDY.datasetinfo)
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
            %## CREATE AVERAGE ERSP MEASURES
            tt = cellfun(@(x) mean(abs(x),3),{tmp.itc_dat_slide},'UniformOutput',false);
            [tmp(:).itc_dat_slide] = tt{:}; %abs(mean(tmp.itc_dat_slide))
            %--
            tt = cellfun(@(x) mean(abs(x),3),{tmp.ersp_dat_slide},'UniformOutput',false);
            [tmp(:).ersp_dat_slide] = tt{:}; %abs(mean(tmp.itc_dat_slide))
            %--
            tt = cellfun(@(x) mean(abs(x),3),{tmp.erspb_dat_slide},'UniformOutput',false);
            [tmp(:).erspb_dat_slide] = tt{:}; %abs(mean(tmp.itc_dat_slide))

            %## RMV FIELDS & STORE
            % tmp = rmfield(tmp,{'itc_freqs','itc_times','ersp_dat_slide','erspb_dat_slide','itc_dat'});
            % tmp = rmfield(tmp,{'ersp_dat_slide','erspb_dat_slide','itc_dat'});
            tmp = rmfield(tmp,{'itc_dat'});
            tmp = struct2table(tmp);            
            itc_so{subj_i} = tmp;
            
            %## CREATE AVERAGE ERSP MEASURES
            % i = 1;
            % b_i = 1;
            % itc_times = tmp(i).itc_times;
            % itc_freqs = tmp(i).itc_freqs;
            % dat = tmp(i).ersp_dat_slide(:,:,b_i);
            % dat = tmp(i).erspb_dat_slide(:,:,b_i);
            % dat = mean(abs(tmp(i).itc_dat_slide(:,:,b_i)),3);
            % dat = std(abs(tmp(i).itc_dat_slide(:,:,:)),[],3);
            % dat = tmp(i).itc_dat;
            % figure;
            % contourf(itc_times,itc_freqs,abs(dat),30,'lines','none');
            % set(gca,'fontsize',14);
            % ylabel('Frequency (Hz)');
            % xlabel('Time (ms)');
            % colorbar;
        end
    catch e 
        fprintf('%s) %s\n',tmp_cl_study.datasetinfo(subj_i).subject ,getReport(e));
    end
end
%--
itc_so = util_resolve_table(itc_so);
%-- remove all entries of itc_freqs and itc_times except 1
tmp = num2cell(zeros(size(itc_so,1)-1,1));
itc_so.itc_freqs(2:size(itc_so,1)) = tmp(:);
itc_so.itc_times(2:size(itc_so,1)) = tmp(:);
itc_so.mi_freqs(2:size(itc_so,1),:) = zeros(size(itc_so,1)-1,size(itc_so.mi_freqs(1,:),2)); %tmp(:);
itc_so.mi_phase(2:size(itc_so,1),:) = zeros(size(itc_so,1)-1,size(itc_so.mi_phase(1,:),2));
%--
% itc_so = util_resolve_table(itc_so);
par_save(itc_so,save_dir,sprintf('itc_table_%s.mat',fext));
%(04/30/2025) JS, NH3128 gets an OOM error on the fext='phasec_notw_mw'
%run. Need to include them for the final analysis

%% (SAVE R TABLE) =================================================== %%
% fext = 'phasec_notw_mw';
% fext = 'phasec_notw_mw_based';
fext = 'phasec_notw_mw_based';
fext_save = 'crop_phasec_notw_based';
CL_NUM_CUTOFF = 13;
GROUP_CHARS = {'H1000','H2000','H3000'};
COND_CHARS = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
FREQ_BOUND = [3,80];
% FREQ_BOUND = [3,250];
TIME_BOUND = [0,1500];

%## LOOP
% itc_so = cell(11,length(CL_STUDY.datasetinfo));
itc_so = cell(length(CL_STUDY.datasetinfo),1);
parfor subj_i = 1:length(CL_STUDY.datasetinfo)
    tmp_sbs_study = SBS_STUDY;
    tmp_cl_study = CL_STUDY;
    subj_char = tmp_sbs_study.datasetinfo(subj_i).subject;
    subj_ind = strcmp(tmp_cl_study.datasetinfo(subj_i).subject,{tmp_sbs_study.datasetinfo.subject});
    tmp_tbnd = TIME_BOUND;
    tmp_fbnd = FREQ_BOUND;

    %## TRY
    try      
        if ~isempty(subj_ind) && any(subj_ind)
            itc_dato = par_load([tmp_cl_study.datasetinfo(subj_i).filepath filesep sprintf('%s_custom_itc_%s.mat',tmp_cl_study.datasetinfo(subj_i).subject,fext)]);
            %--
            fprintf('%s) Adding ITC structure\n',tmp_cl_study.datasetinfo(subj_i).subject);
            tmp = itc_dato.itc_dat_struct;
            twp = itc_dato.parameters;
            %--
            alltimes = tmp(1).itc_times';
            allfreqs = tmp(1).itc_freqs';
            tinds = alltimes > tmp_tbnd(1) & alltimes < tmp_tbnd(2);
            finds = allfreqs > tmp_fbnd(1) & allfreqs < tmp_fbnd(2);
            %--
            tclu = unique([tmp.cluster_n]);
            % tcu = unique({tmp.cond_char});
            % dats_store = cell(length(tclu)*length(COND_CHARS),11);
            dats_store = cell(length(tclu)*length(COND_CHARS),1);
            cnt = 1;
            tti = tic();
            for cl_i = 1:length(tclu)
                if tclu(cl_i) <= CL_NUM_CUTOFF
                    for c_i = 1:length(COND_CHARS)
                        inds = [tmp.cluster_n] == tclu(cl_i) & ...
                            strcmp({tmp.cond_char},COND_CHARS{c_i});
                        tt = tmp(inds);
                        %-- convert & take abs of data
                        % rtd = abs(tt.itc_dat_slide(finds,tinds,:));
                        % rtd = abs(tt.erspb_dat_slide(finds,tinds,:));
                        % rtd = abs(tt.ersp_dat_slide(finds,tinds,:));
                        % rtd = abs(tt.itc_dat(finds,tinds));
                        %-- one image
                        % rtd = mean(abs(tt.itc_dat_slide(finds,tinds,:)),3);
                        %-- multi image
                        rtd = abs(tt.itc_dat_slide(finds,tinds,:));
                        rtd = mean(rtd,3);
                        %## STORE DATA
                        gg = strcmp(tt.group_char,GROUP_CHARS);
                        cc = strcmp(tt.cond_char,COND_CHARS);
                        if cnt == 1
                            tmp_dats = struct(...
                                'subj_char',tt.subj_char, ...
                                'subj_n',subj_i, ...
                                'group_char',tt.group_char, ...
                                'group_n',find(gg), ...
                                'cond_char',tt.cond_char, ...
                                'cond_n',find(cc), ...
                                'mod_n',tt.mod_n, ...
                                'cluster_n',tt.cluster_n, ...
                                'itc_dat',rtd, ...
                                'itc_freq',tt.itc_freqs(finds), ...
                                'itc_time',tt.itc_times(tinds) ...
                            );
                            % tmp_dats(cnt,:) = {tt.subj_char,subj_i,tt.group_char,find(gg),tt.cond_char, ...
                            %     find(cc),tt.mod_n,tt.cluster_n,rtd,tt.itc_freqs(finds),tt.itc_times(tinds)};
                        else
                            tmp_dats = struct(...
                                'subj_char',tt.subj_char, ...
                                'subj_n',subj_i, ...
                                'group_char',tt.group_char, ...
                                'group_n',find(gg), ...
                                'cond_char',tt.cond_char, ...
                                'cond_n',find(cc), ...
                                'mod_n',tt.mod_n, ...
                                'cluster_n',tt.cluster_n, ...
                                'itc_dat',rtd, ...
                                'itc_freq',0, ...
                                'itc_time',0 ...
                            );
                            % tmp_dats = {tt.subj_char,subj_i,tt.group_char,find(gg),tt.cond_char, ...
                            %     find(cc),tt.mod_n,tt.cluster_n,rtd,0,0};
                        end
                        dats_store{cnt} = tmp_dats;
                        % dats_store(cnt,:) = tmp_dats;
                        cnt = cnt + 1;
                    end
                end
            end
            % tmp = util_resolve_struct(dats_store);
            % tmp = dats_store(all(~cellfun(@isempty,dats_store),2),:);
            tmp = dats_store(~cellfun(@isempty,dats_store));
            itc_so{subj_i} = tmp;
            fprintf('%s) loading and conversion done. %0.2f min\n',tmp_cl_study.datasetinfo(subj_i).subject,toc(tti)/60)
        end
    catch e 
        fprintf('%s) %s\n',tmp_cl_study.datasetinfo(subj_i).subject ,getReport(e));
    end
end
%--
itc_so = itc_so(~cellfun(@isempty,itc_so));
tmp = itc_so;
tmp = cat(1,tmp{:});
%--
save([save_dir filesep sprintf('rdata_struct_ext_%s.mat',fext_save)],'tmp');
% par_save(tmp,save_dir,sprintf('rdata_struct_%s.mat',fext_save));
% tmp = par_load(save_dir,sprintf('rdata_struct_%s.mat',fext));
%% (SAVE R TABLE SPCA ERSP) =================================================== %%
fext = 'ersp_spca';
%- spca dir
SPCA_STUDY_DNAME = '02202025_mim_yaoa_spca_calcs';
STUDY_FNAME_GAIT = 'spca_gait_epoch_study_all';
spca_dir = [studies_fpath filesep sprintf('%s',SPCA_STUDY_DNAME)];
%- gait
if ~ispc
    tmp = load('-mat',[spca_dir filesep sprintf('%s_UNIX.study',STUDY_FNAME_GAIT)]);
    SPCA_STUDY = tmp.STUDY;
else
    tmp = load('-mat',[spca_dir filesep sprintf('%s.study',STUDY_FNAME_GAIT)]);
    SPCA_STUDY = tmp.STUDY;
end

%## GET TIMEWARPING INFORMATION
tmpf = par_load(SPCA_STUDY.datasetinfo(1).filepath,'gait_ersp_spca.mat');
timef_params = tmpf.icatimefopts;
timef_params.timewarpms = tmpf.warptimes;
hardcode_times = tmpf.icatimefopts.times;
hardcode_freqs = tmpf.icatimefopts.freqs;
%-- bounds and crops
time_bound = [timef_params.timewarpms(1),timef_params.timewarpms(end)];
finds = find(hardcode_freqs>=FREQ_BOUND(1) & hardcode_freqs<=FREQ_BOUND(2));
tinds = find(hardcode_times>=time_bound(1) & hardcode_times<=time_bound(2));

SPCA_TABLE = par_load(cluster_k_dir,'spca_cluster_table_ersp.mat');
alltimes = hardcode_times(tinds);
allfreqs = hardcode_freqs(finds);

%## REMOVE UNNECESSARY INFORMATION
% meas_del = 'tf_gpmcorr_c';
meas_del = 'tf_erspcorr_c'; 
outs = table2struct(SPCA_TABLE);
outs = rmfield(outs,{'tf_ersporig_c','tf_gpmorig_c','tf_pc1_c','tf_coeff_c', ...
    'comp_c','cluster_n',meas_del});
outc = cell2struct(struct2cell(outs), ...
    {'subj_char','group_char','cluster_n','group_n','cond_char','itc_dat'});

%## CONVERT TABLE TO RDATA
ter_cs = {'flat','low','med','high'};
sp_cs = {'0p25','0p5','0p75','1p0'};
all_cs = [ter_cs,sp_cs];
%-- add subj_n
subj_n = cell(length(outc),1);
for subj_i = 1:length(CL_STUDY.datasetinfo)
    tmps = CL_STUDY.datasetinfo(subj_i).subject;
    inds = strcmp(tmps,{outc.subj_char});
    subj_n(inds) = repmat({subj_i},[sum(inds),1]);    
end
%-- add mod_n
mod_n = cell(length(outc),1);
indst = cellfun(@(x) strcmp(x,{outc.cond_char}),ter_cs,'UniformOutput',false);
indst = any(cat(1,indst{:}),1);
indss = cellfun(@(x) strcmp(x,{outc.cond_char}),sp_cs,'UniformOutput',false);
indss = any(cat(1,indss{:}),1);
mod_n(indst) = repmat({1},[sum(indst),1]);
mod_n(indss) = repmat({2},[sum(indss),1]);
%-- add cond_n
cond_n = cell(length(outc),1);
for c_i = 1:length(all_cs)
    inds = strcmp(all_cs{c_i},{outc.cond_char});
    cond_n(inds) = repmat({c_i},[sum(inds),1]);
end
%-- crop ersp dat
dd = cellfun(@(x) x(finds,tinds),{outc.itc_dat},'UniformOutput',false);
%-- remake struct
[outc(:).subj_n] = subj_n{:};
[outc(:).mod_n] = mod_n{:};
[outc(:).cond_n] = cond_n{:};
[outc(1).itc_freq] = allfreqs;
[outc(1).itc_time] = alltimes;
titt = repmat({0},[length(outc)-1,1]);
[outc(2:length(outc)).itc_freq] = titt{:};
[outc(2:length(outc)).itc_time] = titt{:};
[outc(:).itc_dat] = dd{:};
%-- put each element into a cell
tmp = cell(length(outc),1);
for o_i = 1:length(outc)
    tmp{o_i} = outc(o_i);
end
outc = tmp;
%--
par_save(outc,save_dir,sprintf('rdata_struct_%s.mat',fext));

%% ===================================================================== %%
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
%% IMPORT R FLASSO ===================================================== %%
% fext = 'itc_rdata_table_phasec_notw_mw_based';
% fext = sprintf('%s_tables_figs',fext);
% fext = 'itc_rdata_table_phasec_notw_mw_based_flasso_results_bsz5';
% fext = 'itc_rdata_table_phasec_notw_mw_based_fl_res_bsz5_nob';
% fext = 'itc_rdata_flasso_out_bsz5_nob_sliding';
% fext = 'itc_rdata_flasso_125f_out_bsz5_nob_sliding';
% fext = 'rdata_extc_phasec_notw_condb';
fext = 'rdata_extc_phasec_notw_nocondb';
COND_CHARS = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
GROUP_CHARS = {'H1000','H2000','H3000'};
clusters = [3,4,6,8,5,9,11];
% conds = [5,6,7,8];
conds = [1,2,3,4];
%## LOOP
stats_struct = struct('itc_dat',[], ...
    'itc_mod',[], ...
    'lamb_val',[], ...
    'itc_freqs',[], ...
    'itc_times',[], ...
    'subj_char',{''}, ...
    'cond_char',{''}, ...
    'cluster_n',[], ...
    'group_char',[]);
stats_struct = repmat(stats_struct,[1,length(COND_CHARS)*length(GROUP_CHARS)*90]);
cnt = 1;
for cl_i = 1:length(clusters)
    
    for c_i = 1:length(conds)
        flasso_fpath = [r_stats_dir filesep fext filesep sprintf('allmat_cl%i-c%i.mat',clusters(cl_i),conds(c_i))];
        tmp = load(flasso_fpath);
        tmp = tmp.flasso_tbl;
        subjs = tmp.svec;
        for s_i = 1:length(subjs)
            tmp_cl_study = CL_STUDY;
            %--
            itc_freqs = tmp.freqs;
            itc_times = tmp.times;
            %--            
            stats_struct(cnt).subj_char = tmp_cl_study.datasetinfo(subjs(s_i)).subject;
            stats_struct(cnt).cond_char = COND_CHARS{conds(c_i)};
            stats_struct(cnt).group_char = tmp_cl_study.datasetinfo(subjs(s_i)).group;
            stats_struct(cnt).cluster_n = clusters(cl_i);
            stats_struct(cnt).itc_dat = reshape(tmp.odat(:,s_i),length(itc_freqs),length(itc_times));
            stats_struct(cnt).itc_mod = reshape(tmp.mmat(:,s_i),length(itc_freqs),length(itc_times));
            stats_struct(cnt).lamb_val = tmp.lvec(s_i);
            stats_struct(cnt).itc_freqs = itc_freqs;
            stats_struct(cnt).itc_times = itc_times;
            %--     
            % dat = reshape(tmp.mmat(:,s_i),length(itc_freqs),length(itc_times));
            % dat = reshape(tmp.odat(:,s_i),length(itc_freqs),length(itc_times));
            % figure;
            % contourf(itc_times,itc_freqs,dat,30,'lines','none');
            % set(gca,'fontsize',14);
            % ylabel('Frequency (Hz)');
            % xlabel('Time (ms)');
            % colorbar;
            cnt = cnt + 1;
        end
    end
end

%## CONVERT TO TABLE
stats_struct = stats_struct(~cellfun(@isempty,{stats_struct.subj_char}));
tmp = struct2table(stats_struct);
par_save(tmp,save_dir,sprintf('itc_rdata_table_%s_cl%s_c%s.mat',fext,strjoin(string(conds),''),strjoin(string(clusters),'')));
