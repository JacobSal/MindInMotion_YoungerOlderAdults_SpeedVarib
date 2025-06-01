function [paramso] = local_psd_plot_inter(ax,psd_dat_in,params,PLOT_STRUCT,LINE_STRUCT,varargin)
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB R2020b
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, liu.chang1@ufl.edu
% cat_logo();
% %-
% PRINT_CHKS = true;
% % alltitles = 'ERSP';
% tmpfg = gcf();
% psr = tmpfg.Position(4)/tmpfg.Position(3); % h/w (e.g., 9/6.5 for letter paper
% %-
% label_struct = struct('FontName','Arial', ...
%     'FontSize',8, ...
%     'FontWeight','bold');
% DEFAULT_PLOT_STRUCT = struct( ...
%     'title',{'ERSP'},...
%     'title_props',label_struct, ...
%     'xlabel','Gait Cycle Time (ms)',...
%     'xlabel_props',label_struct, ...
%     'ylabel','Frequency (Hz)',...
%     'ylabel_props',label_struct, ...
%     'xticklabel_times',[],...
%     'xticklabel_chars',{{}},...
%     'xticklabel_angle',45,...
%     'clim',[],...
%     'font_size',8,...
%     'font_name','Arial',...
%     'freq_lims',[],...
%     'time_lims',[],...    
%     'contourf_grain',ceil((500/pi())),...
%     'alpha_multiple',0.6,...
%     'position',[0.06,0.7,0.13,0.16], ...
%     'do_display_sgtitle',true, ...
%     'sgtitle_char',{''}, ...
%     'sgtitle_shift',[0,0.65], ...
%     'sgtitle_boxsz',[0.1,0.1], ...
%     'sgtitle_props',struct(...
%         'LineStyle','none',...
%         'FontName','Arial', ...
%         'FontSize',12,...
%         'FontWeight','bold',...
%         'HorizontalAlignment','center', ...
%         'VerticalAlignment','top',...
%         'Units','normalized'), ...
%     'do_display_colorbar',true, ...
%     'cbar_shift',[0.8*psr,psr*(1/psr)],...
%     'cbar_label_shift',[1.65*psr,1.27*(1/psr)],...
%     'cbar_ticks',[],...
%     'cbar_ticklabs',{{''}},...
%     'cbar_label','\Delta Power (dB)',...
%     'cbar_label_props',struct( ...
%         'LineStyle','none',...
%         'FontName','Arial', ...
%         'FontSize',8,...
%         'FontWeight','bold',...
%         'HorizontalAlignment','center', ...
%         'VerticalAlignment','middle',...
%         'Units','normalized',...
%         'Rotation',270), ...
%     'do_display_bandmarks',true,...
%     'bandmarks',{{'\theta','\alpha','\beta','\gamma'}},...
%     'bandmarks_shift',[0.075*psr,-0.05*(1/psr); ...
%                        0.075*psr,0.25*(1/psr); ...
%                        0.075*psr,0.53*(1/psr); ...
%                        0.075*psr,0.89*(1/psr)],...
%     'bandmarks_props',struct(...
%         'LineStyle','none',...
%         'FontName','Arial', ...
%         'FontSize',8,...
%         'FontWeight','bold',...
%         'HorizontalAlignment','left', ...
%         'VerticalAlignment','top', ...
%         'Units','normalized'));
% %## 
% FIG_STRUCT = struct();
% %% Define Parser
% p = inputParser;
% %## REQUIRED
% addRequired(p,'ax',@(x) isgraphics(x) || isempty(x));
% addRequired(p,'psd_dat_in',@isnumeric);
% %## PARAMETER
% addParameter(p,'PLOT_STRUCT',DEFAULT_PLOT_STRUCT,@(x) validate_struct(x,DEFAULT_PLOT_STRUCT,PRINT_CHKS));
% addParameter(p,'LINE_STRUCT',DEFAULT_PLOT_STRUCT,@(x) validate_struct(x,DEFAULT_PLOT_STRUCT,PRINT_CHKS));
% addParameter(p,'FIG_STRUCT',DEFAULT_PLOT_STRUCT,@(x) validate_struct(x,DEFAULT_PLOT_STRUCT,PRINT_CHKS));
% parse(p,ax,psd_dat_in,params,PLOT_STRUCT,allersp_mask,allersp_pcond,varargin{:});
% %## SET DEFAULTS
% PLOT_STRUCT = p.Results.PLOT_STRUCT;
% PLOT_STRUCT = set_defaults_struct(PLOT_STRUCT,DEFAULT_PLOT_STRUCT,PRINT_CHKS);
% 
% %## GET FREQUENCY & TIME LIMITS
% if isempty(PLOT_STRUCT.freq_lims)
%     PLOT_STRUCT.freq_lims = [allfreqs(1),allfreqs(end)];
% end
% if isempty(PLOT_STRUCT.time_lims)
%     PLOT_STRUCT.time_lims = [alltimes(1),alltimes(end)];
% end
% 
% %## CONVERT STRUCT ARGUMENTS TO CELL ARGS
% % ylabel_props = struct2args(PLOT_STRUCT.ylabel_props);
% % xlabel_props = struct2args(PLOT_STRUCT.xlabel_props);
% % title_props = struct2args(PLOT_STRUCT.title_props);
% bandmarks_props = struct2args(PLOT_STRUCT.bandmarks_props);
% % sgtitle_props = struct2args(PLOT_STRUCT.sgtitle_props);
%% EXTRACT PSD DATA ==================================================== %%
%--
XTICKS = [4,8,13,30,40];
XFREQ_LINES = [8,13,30];
STAT_ALPHA = 0.05;
% PSD_LINE_WIDTH = 2;
PSD_LINE_WIDTH = 1.5;

%## FUNC PARAMETER SETS
switch params.line_plot_opt
    case 'inter'
        csz = size(psd_dat_in,1);
        gsz = size(psd_dat_in,2);
    case 'group'
        csz = 1;
        gsz = size(psd_dat_in,2);
    case 'cond'
        csz = size(psd_dat_in,1);
        gsz = 1;
end

%--
paramso.stats_store = cell(4,1); 
paramso.ax_store = cell(4,1);
paramso.leg_store = cell(4,1);
cntl = 1;

%## STATS
pcond = {};
pgroup = {};
pinter = {};
%--
if ~isempty(params.stats_in)
    pcond = params.stats_in.pcond;
    pgroup = params.stats_in.pgroup;
    pinter = params.stats_in.pinter;
else
    [pcond, pgroup, pinter, ~, ~, ~] = ...
        std_stat(psd_dat_in, params.stats);
    paramso.stats_in = struct('pcond',{pcond}, ...
        'pgroup',{pgroup}, ...
        'pinter',{pinter});
end
%--
alphac = cell(size(pcond));
alphag = cell(size(pgroup));
alphai = cell(size(pinter));
for g_i = 1:length(pcond)
    for c_i = 1:length(pgroup)
        alphag{c_i} = pgroup{c_i} < STAT_ALPHA;        
    end
    alphai{g_i} = pinter{g_i} < STAT_ALPHA;
    alphac{g_i} = pcond{g_i} < STAT_ALPHA;
end

%## GROUP AVG PLOT
%-- data averaging
psd_in_c = zeros(size(psd_dat_in{1,1},1),csz,gsz);
err_bnd_in = zeros(size(psd_dat_in{1,1},1),csz,gsz,2);

for g_i = 1:gsz
    num_subj = size(psd_dat_in{c_i,g_i},2);
    for c_i = 1:csz
        switch params.line_plot_opt
            case 'inter'
                tmp = psd_dat_in{c_i,g_i};
                psd_in_c(:,c_i,g_i) = mean(tmp,2);
                %--
                % tmp = cat(2,psd_dat_in{:,g_i}); % all conds & subjs
                % tmp = mean(cat(3,psd_dat_in{:,g_i}),3); % mean cond for subjs
                % psd_in_g(:,c_i,g_i) = mean(tmp,2);
                %--
                N = sqrt(num_subj);
                err_bnd_in(:,c_i,g_i,1) = mean(tmp,2) + std(tmp,[],2)/(N);
                err_bnd_in(:,c_i,g_i,2) = mean(tmp,2) - std(tmp,[],2)/(N);
            case 'group'
                tmp = cat(2,psd_dat_in{:,g_i}); % all conds & subjs
                % tmp = mean(cat(3,psd_dat_in{:,g_i}),3); % mean cond for subjs
                psd_in_c(:,c_i,g_i) = mean(tmp,2);
                % N = sqrt(size(tmp,2));
                N = sqrt(num_subj);
                err_bnd_in(:,c_i,g_i,1) = mean(tmp,2) + std(tmp,[],2)/(N);
                err_bnd_in(:,c_i,g_i,2) = mean(tmp,2) - std(tmp,[],2)/(N);
            case 'cond'
                tmp = cat(2,psd_dat_in{c_i,:});
                psd_in_c(:,c_i,g_i) = mean(tmp,2);
                % N = sqrt(size(tmp,2));
                N = sqrt(num_subj);
                err_bnd_in(:,c_i,g_i,1) = mean(tmp,2) + std(tmp,[],2)/(N);
                err_bnd_in(:,c_i,g_i,2) = mean(tmp,2) - std(tmp,[],2)/(N);
        end
        
    end
end

%## PLOTTING PARAMETERS
% c_i = 1;
for g_i = 1:gsz
    for c_i = 1:csz
        tmp_line_struct = LINE_STRUCT;
        tmp_plot_struct = PLOT_STRUCT;
        if c_i == csz && g_i == gsz
            tmp_plot_struct.do_set_ax_props = true;
        else
            tmp_plot_struct.do_set_ax_props = false;
        end
        switch params.line_plot_opt
            case 'inter'
                psd_meanc = squeeze(psd_in_c(:,c_i,g_i));
                %-- line props
                tmp_line_struct.line_props = struct( ...
                    'LineWidth',PSD_LINE_WIDTH, ...
                    'LineStyle',params.line_styles{g_i}, ...
                    'DisplayName',params.xtick_label_c{c_i}, ...
                    'Color',[params.cmaps_scond(c_i,:),0.8] ...
                    );
                tmp_line_struct.err_bnd_props = struct( ...
                    'LineStyle',':', ...
                    'LineWidth',PSD_LINE_WIDTH, ...
                    'FaceAlpha',0.4, ...
                    'EdgeColor','none', ...
                    'FaceColor',params.cmaps_scond(c_i,:));
                tmp_line_struct.err_bnd_vec = [squeeze(err_bnd_in(:,c_i,g_i,1)), ...
                            squeeze(err_bnd_in(:,c_i,g_i,2))];
            case 'group'
                psd_meanc = squeeze(psd_in_c(:,c_i,g_i));
                %-- line props
                tmp_line_struct.line_props = struct( ...
                    'LineWidth',PSD_LINE_WIDTH, ...
                    'LineStyle','-', ...
                    'DisplayName',params.xtick_label_g{g_i}, ...
                    'Color',[params.cmaps_sgroup(g_i,:),0.8] ...
                    );
                tmp_line_struct.err_bnd_props = struct( ...
                    'LineStyle',':', ...
                    'LineWidth',3, ...
                    'FaceAlpha',0.4, ...
                    'EdgeColor','none', ...
                    'FaceColor',params.cmaps_sgroup(g_i,:));
                %--
                tmp_line_struct.err_bnd_vec = [squeeze(err_bnd_in(:,c_i,g_i,1)), ...
                            squeeze(err_bnd_in(:,c_i,g_i,2))];
                
            case 'cond'        
                psd_meanc = squeeze(psd_in_c(:,c_i,g_i));
                %-- line props
                tmp_line_struct.line_props = struct( ...
                    'LineWidth',PSD_LINE_WIDTH, ...
                    'LineStyle','-', ...
                    'DisplayName',params.xtick_label_c{c_i}, ...
                    'Color',[params.cmaps_scond(c_i,:),0.8] ...
                    );
                tmp_line_struct.err_bnd_props = struct( ...
                    'LineStyle',':', ...
                    'LineWidth',PSD_LINE_WIDTH, ...
                    'FaceAlpha',0.4, ...
                    'EdgeColor','none', ...
                    'FaceColor',params.cmaps_scond(c_i,:));
                tmp_line_struct.err_bnd_vec = [squeeze(err_bnd_in(:,c_i,g_i,1)), ...
                            squeeze(err_bnd_in(:,c_i,g_i,2))];
        end
        %## PLOT
        [~,~,Li] = plot_psd(ax,psd_meanc,params.freqs, ...
            'LINE_STRUCT',tmp_line_struct, ...
            'PLOT_STRUCT',tmp_plot_struct);
        paramso.leg_store{cntl} = Li;
        cntl = cntl + 1;
        %--
        hold on;
    end
end
% npos = get(ax,'Position',PLOT_STRUCT.ax_props.Position);

%## PLOT BANDMARKS
% params.bandmark_shifts = [0.7,0.7,0.7,0.7;
%     0.235,0.325,0.55,0.85];
% BM_HZ = 0.7;
% BM_VR = [0.235,0.325,0.55,0.85];
%--
% BM_HZ = 0.9;
% BM_VR = [0.155,0.25,0.49,0.75];
BMO = params.bandmark_shifts;
%--
tmpfg = gcf();
psr = tmpfg.Position(4)/tmpfg.Position(3);
bm_struct = struct(...
    'do_display_bandmarks',true,...
    'bandmarks',{{'\theta','\alpha','\beta','\gamma'}},...
    'bandmarks_shift',[BMO(2,1)*psr,BMO(1,1)*(1/psr); ...
                       BMO(2,2)*psr,BMO(1,2)*(1/psr); ...
                       BMO(2,3)*psr,BMO(1,3)*(1/psr); ...
                       BMO(2,4)*psr,BMO(1,4)*(1/psr)],...
    'bandmarks_props',struct(...
        'LineStyle','none',...
        'FontName','Arial', ...
        'FontSize',PLOT_STRUCT.ax_props.FontSize*1.25,...
        'FontWeight','bold',...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'Units','normalized'));
%--
bandmarks_props = struct2args(bm_struct.bandmarks_props);
npos = PLOT_STRUCT.ax_props.Position;
ANN_DIM = [0.1,0.1];
if bm_struct.do_display_bandmarks
    npos = get(ax,'Position');
    for ii = 1:length(bm_struct.bandmarks)
        xx = npos(1)+npos(3)*bm_struct.bandmarks_shift(ii,1)-ANN_DIM(1)/2;
        yy = npos(2)+npos(4)*bm_struct.bandmarks_shift(ii,2)-ANN_DIM(2)/2;
        a = annotation(gcf,'textbox',[xx,yy,0.1,0.1],...
            'String',bm_struct.bandmarks{ii}, ...
            bandmarks_props{:});
    end
end

%## SET Y-LABEL PROPS
set(ax,'YLim',params.y_lims, ...
    'YTick',params.yticks, ...
    'YTickLabel',params.ytick_labs);

%## SET X-LABEL PROPS
xtick_labs = cellstr(string(XTICKS));
set(ax,'XLim',[min(XTICKS)-1e-6,max(XTICKS)+1e-6], ...
    'XTick',XTICKS, ...
    'XTickLabel',xtick_labs);

%## BAND LINES
for xx = 1:length(XFREQ_LINES)
    xline(XFREQ_LINES(xx),'--');
end

%% PLOT INTER STATS ==================================================== %%
SHADE_BUFFER = 1;
SIGGAP_SHRINK_SZ = 6;
SHADE_HT_FUDGE = 0.002;
SHADE_HT_FACTOR = 0.05;
SHADE_FAC_ALPHA = 1;
SHADE_EDG_ALPHA = 1;
SHADE_EDG_COLOR = 'none';
STAT_HT_OFFSET_FACTOR = 0;
% stats_titles = {'cond','group','inter'};
stats_titles = {'group','speed','inter'};
%--
y_lims = get(ax,'YLim');
ht = y_lims(2)-y_lims(1);
stat_ht_o = ht*STAT_HT_OFFSET_FACTOR+ht*SHADE_HT_FUDGE;
stat_ht = ht*SHADE_HT_FACTOR+stat_ht_o;
soht_store = zeros(length(alphai),2);
cnt = 1;
%--
fprintf('Plotting group & condition comparison stats\n');
% isempty(alphai)
for i_i = 1:length(alphai)    
    SHADE_FAC_COLOR = params.cmaps_stats(i_i,:);
    reg = alphai{i_i};
    %--
    regt = [0;reg];
    regt = regt == 1;
    reg_st = find(diff(regt) > 0);
    reg_en = find(diff(regt) < 0);
    reg_en = reg_en - 1;
    if length(reg_en)+1 == length(reg_st)
        reg_en = [reg_en;length(reg);length(reg)+1];
    else
        reg_en = [reg_en;length(reg)];
    end    
    %-- shading buffing to ensure visualization
    reg_st = floor(reg_st - SHADE_BUFFER/2);
    reg_st(reg_st < 1) = 1;
    reg_en = ceil(reg_en + SHADE_BUFFER/2);
    reg_en(reg_en > length(reg)) = length(reg);   
    %-- apply stats if true
    fprintf('%i) %i frequencies significant\n',i_i,sum(reg));
    if ~isempty(reg_st)
        %-- ensure clean blocks
        regcc = sort(cat(1,reg_st,reg_en));
        rego = zeros(size(reg_st,1),2);
        cnti = 1;
        cntj = cnti + 1;
        cntr = 1;
        while cnti < length(regcc) && cntj < length(regcc)-1        
            dd2 = regcc(cntj+1)-regcc(cntj);
            while dd2 < SIGGAP_SHRINK_SZ && cntj < length(regcc)-1
                cntj = cntj + 2;
                dd2 = regcc(cntj+1)-regcc(cntj);            
            end
            rego(cntr,1) = regcc(cnti);
            rego(cntr,2) = regcc(cntj);        
            cnti = cntj + 1;
            cntj = cnti + 1;
            cntr = cntr + 1;
        end
        ind = all(rego ~= 0,2);
        rego = rego(ind,:);    
        for r_i = 1:size(rego,1)
            stx = params.freqs(rego(r_i,1));
            enx = params.freqs(rego(r_i,2));
            %-- plot patch
            Pa = patch(ax,[stx,enx,enx,stx], ...
                [y_lims(1)+stat_ht_o,y_lims(1)+stat_ht_o,y_lims(1)+stat_ht,y_lims(1)+stat_ht], ...
                SHADE_FAC_COLOR);
            hold on;
            set(Pa,'edgecolor',SHADE_EDG_COLOR, ...
                'facealpha',SHADE_FAC_ALPHA, ...
                'edgealpha',SHADE_EDG_ALPHA, ...
                'DisplayName',stats_titles{i_i});
            hold on;
            if r_i == 1
                paramso.stats_store{cnt} = Pa;
            end
            % if r_i == 1 
            %     paramso.leg_store{cntl} = Pa;
            % end
            % cntl = cntl + 1;
            cnt = cnt + 1;          
        end        
    end
    %--
    soht_store(i_i,1) = stat_ht;
    soht_store(i_i,2) = stat_ht_o;
    %--
    stat_ht_o = stat_ht;
    stat_ht = stat_ht_o + ht*SHADE_HT_FACTOR + ht*SHADE_HT_FUDGE;
end
%(04/02/2025) JS, so far the leat buggy stats algorithm.

%## OUTPUTS
paramso.stats_store = paramso.stats_store(~cellfun(@isempty,paramso.stats_store));
paramso.leg_store = paramso.leg_store(~cellfun(@isempty,paramso.leg_store));

%## LEGEND ============================================================== %%
if params.do_display_leg
    %-- legend
    leg_store = [paramso.leg_store{:};paramso.stats_store{:}];
    legend(gca,leg_store);
    [lg2,~,~,~]  = legend('boxoff');
    tmp = get(lg2,'String');
    %-- specs
    set(lg2,'String',tmp, ...
        'Position',[lg2.Position(1)+params.leg_position(1),...
            lg2.Position(2)+params.leg_position(2), ...
            lg2.Position(3), ...
            lg2.Position(4)], ...
        params.legends_specs{:});
    lg2.ItemTokenSize(1) = params.leg_token_size;
    hold off;
end

end
%% ===================================================================== %%
function [stat_out] = stat_psd_interaction(psd_dat)
%(04/09/2025) JS, the original call to std stat seems to make calls to:
% statcondfieldtrip -> ft_freqstatistics -> ft_statistics_montecarlo -> ...
% for condition wise: ft_statfun_depsamplesFunivariate
% for group wise: ft_statfun_indepsamplesF
% for interaction: ft_statfun_indepsamplesF?
% this fieldtrip guide seems to appraoch the interaction differently using 
% ft_statfun_indepsamplesT (maybe in my case this would be
% ft_statfun_depsamplesT?) 

cfg = [];
cfg.parameter = {'powspctrm','powspctrm_b'};
cfg.operation = 'subtract';

% sedation_respon_d = ft_math(cfg, base_sedation_respon, mode_sedation_respon);
% sedation_drowsy_d = ft_math(cfg, base_sedation_drowsy, mode_sedation_drowsy);

cfg = [];
cfg.channel          = 'all';
cfg.frequency        = foi_contrast;
cfg.avgovergfreq     = 'no';
cfg.parameter        = 'powspctrm_b';
cfg.method           = 'ft_statistics_montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.clusterthreshold = 'nonparametric_common';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = cfg.tail;
cfg.alpha            = 0.05;
cfg.correcttail      = 'alpha';
cfg.computeprob      = 'yes';
cfg.numrandomization = 1000;
cfg.neighbours       = cfg_neigh.neighbours;

gind = 2;
cind = 1;
sz = size(psd_dat);

design = zeros(1,size(respon_group,1) + size(drowsy_group,1));
design(1,1:size(respon_group,1)) = 1;
design(1,(size(respon_group,1)+1):(size(respon_group,1)+size(drowsy_group,1))) = 2;

cfg.design = design;
cfg.ivar   = 1;

%## design sample: this is what EEGLAB puts in for the
%ft_statfun_depsamplesFunivariate
%--
% cfg.mcorrect = 'fdr';
% cfg.numrandomization = 4000;
% cfg.method = montecarlo;
% cfg.correctm = 'fdr';
% cfg.feedback = 'no';
% cfg.alpha = 0.0500;
% cfg.correcttail = 'no';
% cfg.tail = 1;
% cfg.statistic = 'depsamplesFunivariate';
% cfg.ivar = 1;
% cfg.uvar = 2;
% design = zeros(2,size(psd_dat,1)*size(psd_dat{1},2)) % for 4 conditions
% and 23 subjects this would be length 92
% cnt = 1;
% for i = 1:size(psd_dat,1)
%     for j = 1:size(psd_dat{1},2)
%         design(1,cnt) = i; % condition number
%         design(2,cnt) = j; % subject number
%     end
% end
% cfg.design = design;

stat = ft_freqstatistics(cfg, sedation_respon_d, sedation_drowsy_d);

end

