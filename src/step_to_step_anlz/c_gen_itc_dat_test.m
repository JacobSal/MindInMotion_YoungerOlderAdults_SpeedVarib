%   Project Title: MIM OA & YA SPEED & KINETICS ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/step_to_step_anlz/run_sts_c_gen_itc_dat.sh

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
%% (PARAMETERS) ======================================================== %%
FORCE_RECALC_ERSP = true;
ERSP_PARAMS = struct('subbaseline','off',...
    'timerange',[],...
    'ersplim',[-2,2],...
    'freqfac',4,...
    'cycles',[3,0.8],...
    'freqrange',[1,200]);
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
NEWTIMEF_STRUCT = struct(...
    'timewindow',[], ...
    'timewarp',[], ...
    'timewarpms',[], ...
    'ntimesout',[], ...
    'pointrange',[], ...
    'baseline',nan(), ...
    'timelimits',[], ...
    'freqs',[3,250], ...
    'freqscale','log', ...
    'do_plot_ersp',false, ...
    'do_plot_itc',false, ...
    'do_plot_phase',false, ...
    'cycles',[3,0.8], ...
    'padratio',1, ...
    'nfreqs',200, ...
    'alpha',nan(), ...
    'srate',[] ...
    );
%% (PATHS)
%- datset name
DATA_SET = 'MIM_dataset';
%- study name
STUDY_DNAME = '02202025_mim_yaoa_powpow0p3_crit_speed';
STUDY_FNAME = 'kin_eeg_epoch_study';
STUDY_EPOCH_FNAME = 'epoch_study';
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
%## LOAD EPOCH STUDY
% need this for the gait timewarp information
if ~ispc
    tmp = load('-mat',[studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep sprintf('%s_UNIX.study',STUDY_EPOCH_FNAME)]);
    E_STUDY = tmp.STUDY;
else
    tmp = load('-mat',[studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep sprintf('%s.study',STUDY_EPOCH_FNAME)]);
    E_STUDY = tmp.STUDY;
end
%% (STUDY DESIGN) ====================================================== %%
% condition_gait = unique({STUDY.datasetinfo(1).trialinfo.cond}); 
% % (09/18/2024) JS, reorders it weirdly. Manually override
%## ersp plot per cluster per condition
args = eeglab_struct2args(ERSP_STAT_PARAMS_COND);
SBS_STUDY = pop_statparams(SBS_STUDY,args{:});
args = eeglab_struct2args(ERSP_PARAMS);
SBS_STUDY = pop_erspparams(SBS_STUDY,args{:});
SBS_STUDY.cache = [];
for des_i = 1:length(STUDY_DESI_PARAMS)
    [SBS_STUDY] = std_makedesign(SBS_STUDY,[],des_i,STUDY_DESI_PARAMS{des_i}{:});
end

%## REASSIGN CLUSTER
cl_struct = par_load(cluster_k_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
SBS_STUDY.cluster = cl_struct;
[comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(SBS_STUDY);

%## TIME WARP
SRATE = 500;
TIMEWARP_NTIMES = floor(SRATE/pi); % conservative nyquist frequency. making this too big can cause overlap between gait cyles
averaged_warpto_events = E_STUDY.etc.a_epoch_process.avg_warpto_events;
ERSP_CROP_TIMES = [averaged_warpto_events(1), averaged_warpto_events(end)+1];
disp(['Grand average (across all subj) warp to: ',num2str(averaged_warpto_events)]);
%--
CL_STUDY.cluster = cl_struct;
%% ERSP
% parfor s_i = 1:length(SBS_STUDY.datasetinfo)
for s_i = 1:length(CL_STUDY.datasetinfo)
    %## PARAM COPIES
    % tmp_study_sbs = SBS_STUDY;
    tmp_cl_study = CL_STUDY;
    tmp_ersp_params = ERSP_PARAMS;
    % tmp_ntf_struct = NEWTIMEF_STRUCT;
    try
        EEG = pop_loadset('filepath',tmp_cl_study.datasetinfo(s_i).filepath, ...
            'filename',tmp_cl_study.datasetinfo(s_i).filename);
        fprintf('Running Subject %s\n',EEG.subject);
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        end
        ics_calc = EEG.etc.urreject.ic_keep;
        % dat_in = EEG.icaact(ics_calc,:,:);
        % EEG.icaact = dat_in;
        % dat_in = double.empty;
        %(04/28/2025) JS, overriding EEG.icaact due to a bug higher up (or
        %purposeful overlook?) where I use pop_select instead of
        %pop_subcomp... If I try to load icaact it just gives me 0's
        % Just using cluster study for now. May need to fix
        % later.

        %## PARAMETERS IN
        tlimits = [EEG(1).xmin EEG(1).xmax]*1000;
        pointrange1 = round(max((tlimits(1)/1000-EEG(1).xmin)*EEG(1).srate, 1));
        pointrange2 = round(min(((tlimits(2)+1000/EEG(1).srate)/1000-EEG(1).xmin)*EEG(1).srate, EEG(1).pnts));
        pointrange = (pointrange1:pointrange2);
        
        %-- storage
        MOD_CHARS = {{'flat','low','med','high'},{'0p25','0p5','0p75','1p0'}};
        ITC_DAT_STRUCT = struct(...
            'subj_char',{''}, ...
            'group_char',{''}, ...
            'cond_char',{''}, ...
            'mod_n',[], ...
            'comp_n',[], ...
            'cluster_n',[], ...
            'itc_dat',[], ...
            'itc_times',[], ...
            'itc_freqs',[] ...
            );
        
        %-- newtimef func params
        ntf_struct = struct(...
            'timewindow',[], ...
            'timewarp',EEG.timewarp.latencies, ...
            'timewarpms',averaged_warpto_events, ...
            'ntimesout',floor(EEG.srate/pi()), ...
            'pointrange',pointrange, ...
            'baseline',nan(), ...
            'timelimits',[EEG(1).xmin,EEG(1).xmax]*1000, ...
            'freqs',[3,250], ...
            'itctype','phasecoher', ...  % this is the measure AS used
            'freqscale','log', ...
            'do_plot_ersp',false, ...
            'do_plot_itc',false, ...
            'do_plot_phase',false, ...
            'do_single_dat',false, ...
            'do_timewarp',false, ...
            'cycles',[3,0.8], ...
            'padratio',4, ...
            'nfreqs',length((tmp_ersp_params.freqrange(1):tmp_ersp_params.freqrange(2))), ...
            'alpha',nan(), ...
            'srate',EEG.srate ...
        );

        %-- get conditiion info
        cond_chars = unique({EEG.event.cond});
        ind = find(strcmp(EEG.subject,{tmp_cl_study.datasetinfo.subject}));
        trialinfo = std_combtrialinfo(tmp_cl_study.datasetinfo,ind);
        %-- iters
        itersk = zeros(size(EEG.icaact,1)*length(cond_chars),1);
        itersc = zeros(size(EEG.icaact,1)*length(cond_chars),1);
        iterscl = zeros(size(EEG.icaact,1)*length(cond_chars),1);
        %-- data store
        itc_dato = cell(size(EEG.icaact,1)*length(cond_chars),1);
        sdinfo = {tmp_cl_study.datasetinfo.subject};
        scinfo = tmp_cl_study.cluster;
        %-- iter assign
        cnt = 1; 
        for k_i = 1:size(EEG.icaact,1)
            for c_i = 1:length(cond_chars)
                cc = false;
                cl_i = 2;
                while ~any(cc) && cl_i < length(scinfo)
                    cl_i = cl_i + 1;
                    ind = find(strcmp(EEG.subject,sdinfo));
                    ind = scinfo(cl_i).sets == ind;
                    cc = scinfo(cl_i).comps(ind) == k_i;                
                end
                iterscl(cnt) = cl_i;
                itersk(cnt) = k_i;
                itersc(cnt) = c_i;
                cnt = cnt + 1;
            end
        end        
        %--
        fprintf('%s) Running ITC calculation:\n',EEG.subject);
        tt = tic();
        fext = 'phasec_notw_mw_based';
        parfor i = 1:length(itersk)
            %-- init temps
            k_i = itersk(i);
            c_i = itersc(i);
            cl_i = iterscl(i);
            fprintf('%i ',k_i);
            %--
            tmp_cc = cond_chars;
            tmp_ti = trialinfo;
            tmp_eeg = EEG;
            tmp_ntf_struct = ntf_struct;
            tmp_itc_dato = ITC_DAT_STRUCT;
            tmp_ptr = pointrange;
            %--
            % tmp_tls = [-500,1000];
            tmp_tls = tmp_ntf_struct.timelimits;
            tmp_tim = tmp_eeg.times;
            inds = tmp_tim > tmp_tls(1) & tmp_tim < tmp_tls(2);
            tmp_ptr = tmp_ptr(inds);
            lf = floor(pi()*(tmp_eeg.srate/(375)));
            % tmp_tls = [0,sum(abs(tmp_tls))];
            %-- subset conditions
            inds_cond = cellfun(@(x) strcmp(x,tmp_cc{c_i}),{tmp_ti.cond});
            tmp_eeg.icaact = tmp_eeg.icaact(k_i,:,inds_cond);
            % tmp_eeg.icaact = tmp_eeg.icaact(k_i,:,:);
            
            %## SUBSET DATA & RUN
            sp = 25;
            spin = 500*(sp/1000);
            nt = length(tmp_ptr)/spin;
            %--
            % tmp_ntf_struct.timewarp = tmp_eeg.timewarp.latencies(inds_cond,:);
            tmp_ntf_struct.timewarp = tmp_eeg.timewarp.latencies;
            tmp_ntf_struct.do_timewarp = false;
            tmp_ntf_struct.itctype = 'phasecoher'; %'phasecoher';
            tmp_ntf_struct.nfreqs = floor(1/(pi()/tmp_eeg.srate)); %ceil(1*pi()*(2/(length(tmp_ptr)/tmp_eeg.srate)));
            tmp_ntf_struct.padratio =2^2;
            tmp_ntf_struct.cycles = [3,0.8]; %[3,0.8]; %[3,0.8];
            % tmp_ntf_struct.padratio =2^2;
            % tmp_ntf_struct.cycles = 0;
            tmp_ntf_struct.freqs = [lf,tmp_eeg.srate/2];
            tmp_ntf_struct.pointrange = tmp_ptr;
            tmp_ntf_struct.timelimits = tmp_tls;
            tmp_ntf_struct.ntimesout = ceil(nt); %floor(length(tmp_ptr)/(pi()));
            % tmp_ntf_struct.newtimef_params = {...
            %     'ffttaper','hanning', ...
            %     'winsize',round((375/1000)*500,TieBreaker="even"), ... %376 works 
            %     'timesout',ceil(nt)};
            %-- run newtimef
            [e_dat,itc_dat,~,ttimes,lfreqs,~,~,all_e_dat] = ...
                eeglab_newtimef_iter(tmp_eeg,1, ...
                'NEWTIMEF_STRUCT',tmp_ntf_struct);
            % min(ttimes)

            twmed = mean(tmp_eeg.timewarp.latencies(inds_cond,:),1);
            inds = ttimes > twmed(1) & ttimes < twmed(end);
            be = mean(all_e_dat(:,inds,:),2);
            tmpedat = bsxfun(@minus,all_e_dat,be);
            % tmpedat = all_e_dat;
            itcvals = newtimefitc(tmpedat(:,:,:),tmp_ntf_struct.itctype);
            %--
            %(04/29/2025) JS, newtimefitc seems to take advantage of the
            %precomputed time-freq measures from std_precompute, yet still
            %generates the same values as the newtimef "itc" output;
            %however, by default both use "phasecoher" as itctype.
            %## VALIDATION PLOT
            
            tmpd = itc_dat;
            % tmpd = e_dat;
            % tmpd = itcvals;
            %--
            [vals,inds] = unique(lfreqs);
            tmpd = tmpd(inds,:);
            lfreqs = vals;
            %--
            twmed = mean(tmp_eeg.timewarp.latencies(inds_cond,:),1);
            % twmed = [0,200,600,800];
            % xt = twmed;
            % xtl = [twmed(1),twmed(end)];
            % xtt = {'RHS','LTO','LHS','RTO','RHS'}
            %--
            % xt = [ttimes(1),ttimes(end)] 
            % xtl = [ttimes(1),ttimes(end)]
            % xtt = {string(ttimes(1)),string(ttimes(end))};
            %--
            % xt = [ttimes(1),0,ttimes(end)] 
            % xtl = [ttimes(1),ttimes(end)]
            % xtt = {string(ttimes(1)),'RHS',string(ttimes(end))};
            %--
            % xtl = [ttimes(1),ttimes(end)]
            % xt = [ttimes(1),0,200,ttimes(end)] 
            % xtt = {string(ttimes(1)),'RHS','LTO',string(ttimes(end))};
            %--
            % xt = [ttimes(1),500] %ttimes(end)];
            % xtl = [ttimes(1),500] %ttimes(end)];
            % xtt = {'RHS','500 ms'}

            tmp_plot_struct = struct( ...
                'alltitles',{''},...
                'xlabel','Gait Cycle Time (ms)',...
                'ylabel','Frequency (Hz)',...
                'xticklabel_times',xt,...
                'xticklabel_chars',{xtt},...
                'xticklabel_angle',45,...
                'clim',double([min(abs(tmpd),[],'all'),max(abs(tmpd),[],'all')]),...
                'font_size',8,...
                'font_name','Arial',...
                'freq_lims',[3,60],...
                'time_lims',xtl,...    
                'contourf_grain',ceil((500/pi())),...
                'alpha_multiple',0.6,...
                'position',[0.2,0.2,0.5,0.5], ...
                'colorbar_shift',[1,1], ...
                'do_display_bandmarks',true);
            %--
            allersp = abs(tmpd);
            allersp_mask = allersp; %ones(size(tmpd));
            allersp_pcond = zeros(size(allersp));
            % allersp = real(tmpd); 
            alltimes = ttimes;
            allfreqs = lfreqs;
            tmp_plot_struct.clim = double([min(allersp,[],'all'),max(allersp,[],'all')]);
            figure;
            ax = axes();
            [ax] = plot_contourf_cell(ax,allersp,alltimes,allfreqs,...
                allersp_mask,allersp_pcond, ...
                'PLOT_STRUCT',tmp_plot_struct);
            %%
            %## MI INDEX (Tort et al. 2010)
            tmpd = tmp_eeg.icaact(:,tmp_ntf_struct.pointrange,:);
            sigd  = reshape(tmpd,[1,length(tmp_ntf_struct.pointrange)*size(tmpd,3)]);
            [comod,p_vec,f_vec] = mod_index_calc(sigd);

            % figure;
            % contourf(p_vec,f_vec,comod',30,'lines','none');
            % set(gca,'fontsize',14);
            % ylabel('Amplitude Frequency (Hz)');
            % xlabel('Phase Frequency (Hz)');
            % colorbar;

            %## STORE IN STRUCT
            mod_n = cellfun(@(x) any(strcmp(tmp_cc{c_i},x)),MOD_CHARS);   
            twmed = mean(tmp_eeg.timewarp.latencies(inds_cond,:),1);
            %--
            tmp_itc_dato.subj_char = tmp_eeg.subject;
            tmp_itc_dato.group_char = tmp_eeg.group;
            tmp_itc_dato.cond_char = tmp_cc{c_i};
            tmp_itc_dato.mod_n = find(mod_n);
            tmp_itc_dato.subj_char = tmp_eeg.subject;
            tmp_itc_dato.comp_n = k_i;
            tmp_itc_dato.cluster_n = cl_i;
            tmp_itc_dato.itc_dat = single(itc_dat);
            % tmp_itc_dato.e_dat = single(e_dat);
            tmp_itc_dato.itc_times = ttimes';
            tmp_itc_dato.itc_freqs = lfreqs';
            tmp_itc_dato.mi_freqs = f_vec;
            tmp_itc_dato.mi_phase = p_vec;
            tmp_itc_dato.mi_dat =  comod;
            tmp_itc_dato.twc_times = twmed;
            tmp_itc_dato.twg_times = averaged_warpto_events; 
            %-- store
            itc_dato{i} = tmp_itc_dato;

            
            
        end
        fprintf('%s) ITC generation took %0.2f min\n',EEG.subject,toc(tt)/60);
        %## RESHAPE DATA
        % tmp = cat(3,datout{:});
        itc_dato = util_resolve_struct(itc_dato);

        %## ASSIGN DATA & ?SAVE?
        timef_data_struct = struct(...
            'itc_dat_struct',itc_dato, ...
            'datatype','TIMEF', ...
            'datafiles',[EEG.filepath filesep EEG.filename], ...
            'datatrials',ics_calc, ...
            'parameters',ntf_struct, ...
            'trialinfo',trialinfo);
        %-- assign data to individual fields for more memory efficient extraction
        % for k_i = 1:length(g.indices)  % for each (specified) component/channel
        %     timef_data_struct.(sprintf('comp%i',k_i)) = datout{k_i}; 
        % end
        par_save(timef_data_struct,[EEG.filepath filesep sprintf('%s_custom_itc_%s.mat',tmp_cl_study.datasetinfo(s_i).subject,fext)]);
    catch e
        fprintf('%s) Something happened during ERSP generation...\n%s\n',tmp_cl_study.datasetinfo(s_i).subject,getReport(e));
    end
end
