function [paramso] = local_psd_plot_inter(ax,psd_dat_in,params,PLOT_STRUCT,LINE_STRUCT)
%LOCAL_PSD_PLOT Summary of this function goes here
%   Detailed explanation goes here
%% EXTRACT PSD DATA ==================================================== %%
%--
XTICKS = [4,8,13,30,40];
XFREQ_LINES = [8,13,30];
STAT_ALPHA = 0.05;
% YLIM_NTICKS = 5;
% YLIM_SIG_FIGS = 2;
% YLIM_FAC = 0.15;
AX_FONTWEIGHT = 'bold';
%--
% params.line_plot_opt = 'group'; % 'cond', 'inter'
% params.do_display_leg = true;

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
% psd_meanc = squeeze(mean(psd_in_c,2));

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
                % psd_meang = squeeze(psd_in_g(:,c_i,g_i));                                
                %## PLOT
                % if c_i == 1
                %     %--
                %     tmp_line_struct.line_color = params.cmaps_sgroup(g_i,:);
                %     tmp_line_struct.line_label = params.xtick_label_g{g_i};
                %     tmp_line_struct.err_color = params.cmaps_sgroup(g_i,:);
                %     %--
                %     tmp_line_struct.err_bnd_vec = [squeeze(err_bnd_in(:,c_i,g_i,1)), ...
                %             squeeze(err_bnd_in(:,c_i,g_i,2))];
                %     [~,~,~] = plot_psd(ax,psd_meang,params.freqs, ...
                %         'LINE_STRUCT',tmp_line_struct, ...
                %         'PLOT_STRUCT',PLOT_STRUCT);
                % end
                %## PLOT
                %--
                % tmp_line_struct.line_color = params.cmaps_scond(c_i,:);
                % tmp_line_struct.line_style = params.line_styles{g_i};
                % % tmp_line_struct.line_color = params.cmaps_sgroup(g_i,:);
                % tmp_line_struct.line_label = params.xtick_label_c{c_i};
                % tmp_line_struct.err_color = params.cmaps_scond(c_i,:);
                % tmp_line_struct.line_alpha = 0.65;
                % tmp_line_struct.err_bnd_vec = [squeeze(err_bnd_in(:,c_i,g_i,1)), ...
                %             squeeze(err_bnd_in(:,c_i,g_i,2))];

                %-- line props
                tmp_line_struct.line_props = { ...
                    'LineWidth',2, ...
                    'LineStyle',params.line_styles{g_i}, ...
                    'DisplayName',params.xtick_label_c{c_i}, ...
                    'Color',[params.cmaps_scond(c_i,:),0.6], ...
                    };
                % tmp_line_struct.line_color = [cmaps(c_i,:),0.65];
                tmp_line_struct.err_props = { ...
                    'LineStyle',':', ...
                    'LineWidth',3, ...
                    'FaceAlpha',0.6, ...
                    'EdgeColor','none', ...
                    'FaceColor',params.cmaps_scond(c_i,:)};
                tmp_line_struct.err_bnd_vec = [squeeze(err_bnd_in(:,c_i,g_i,1)), ...
                            squeeze(err_bnd_in(:,c_i,g_i,2))]; 
                % tmp_line_struct.do_err_shading = false;
                %--
                [~,~,Li] = plot_psd(ax,psd_meanc,params.freqs, ...
                    'LINE_STRUCT',tmp_line_struct, ...
                    'PLOT_STRUCT',tmp_plot_struct);
                % paramso.ax_store{c_i,g_i} = ax;
                paramso.leg_store{cntl} = Li;
            case 'group'
                psd_meanc = squeeze(psd_in_c(:,c_i,g_i));
                % tmp_line_struct.line_color = params.cmaps_sgroup(g_i,:);
                % tmp_line_struct.line_label = params.xtick_label_g{g_i};
                % tmp_line_struct.err_color = params.cmaps_sgroup(g_i,:);
                %-- line props
                tmp_line_struct.line_props = { ...
                    'LineWidth',2, ...
                    'LineStyle','-', ...
                    'DisplayName',params.xtick_label_g{g_i}, ...
                    'Color',[params.cmaps_sgroup(g_i,:),0.7], ...
                    };
                % tmp_line_struct.line_color = [cmaps(c_i,:),0.65];
                tmp_line_struct.err_props = { ...
                    'LineStyle',':', ...
                    'LineWidth',3, ...
                    'FaceAlpha',0.6, ...
                    'EdgeColor','none', ...
                    'FaceColor',params.cmaps_sgroup(g_i,:)};
                %--
                tmp_line_struct.err_bnd_vec = [squeeze(err_bnd_in(:,c_i,g_i,1)), ...
                            squeeze(err_bnd_in(:,c_i,g_i,2))];
                %-- plot props
                % tmp_plot_struct.ax_props = {...
                %     'box','off', ...
                %     'LineWidth',2, ...
                %     'FontName','Arial', ...
                %     'FontSize',9, ...
                %     'OuterPosition',[0 0 1 1], ...
                %     'Position',[x_shift,y_shift,ax_w,ax_h]};

                %## PLOT
                [~,~,Li] = plot_psd(ax,psd_meanc,params.freqs, ...
                    'LINE_STRUCT',tmp_line_struct, ...
                    'PLOT_STRUCT',tmp_plot_struct);
                % paramso.ax_store{c_i,g_i} = ax;
                paramso.leg_store{cntl} = Li;
            case 'cond'        
                psd_meanc = squeeze(psd_in_c(:,c_i,g_i));
                % tmp_line_struct.line_color = params.cmaps_scond(c_i,:);
                % tmp_line_struct.line_label = params.xtick_label_c{c_i};
                % tmp_line_struct.err_color = params.cmaps_scond(c_i,:);   
                % %--
                % tmp_line_struct.err_bnd_vec = [squeeze(err_bnd_in(:,c_i,g_i,1)), ...
                %             squeeze(err_bnd_in(:,c_i,g_i,2))];

                %-- line props
                tmp_line_struct.line_props = { ...
                    'LineWidth',2, ...
                    'LineStyle','-', ...
                    'DisplayName',params.xtick_label_c{c_i}, ...
                    'Color',[params.cmaps_scond(c_i,:),0.65], ...
                    };
                tmp_line_struct.err_props = { ...
                    'LineStyle',':', ...
                    'LineWidth',3, ...
                    'FaceAlpha',0.6, ...
                    'EdgeColor','none', ...
                    'FaceColor',params.cmaps_scond(c_i,:)};
                tmp_line_struct.err_bnd_vec = [squeeze(err_bnd_in(:,c_i,g_i,1)), ...
                            squeeze(err_bnd_in(:,c_i,g_i,2))];

                %## PLOT
                [~,~,Li] = plot_psd(ax,psd_meanc,params.freqs, ...
                    'LINE_STRUCT',tmp_line_struct, ...
                    'PLOT_STRUCT',tmp_plot_struct);
                % paramso.ax_store{c_i,g_i} = ax;
                paramso.leg_store{cntl} = Li;
        end
        
        cntl = cntl + 1;
        %--
        hold on;
    end
end
% set(ax,'FontWeight',AX_FONTWEIGHT)

%## Y-LIMIT SETTING
% YLIM_NTICKS = 5;
% YLIM_SIG_FIGS_I = 1;
% YLIM_SIG_FIGS = 2;
% YLIM_FAC = 0.15;
% %-- ylim
% y_lims = get(ax,'YLim');
% dy = YLIM_FAC*(y_lims(2)-y_lims(1));
% %-- change 0 point for ylim if exists
% [vy,iy] = min(y_lims);
% if vy == 0 && iy == 1
%     y_lims(iy) = -1e-8;
% elseif vy == 0 && iy == 2
%     y_lims(iy) = 1e-8;
% end
% u = round(y_lims(2),YLIM_SIG_FIGS_I,'significant');
% l = round(y_lims(1),YLIM_SIG_FIGS_I,'significant');
% % u = round((ax_ylim(2)+dy),YLIM_SIG_FIGS_I,'significant');
% % l = round((ax_ylim(1)-dy),YLIM_SIG_FIGS_I,'significant'); 
% bb = (u-l)/(YLIM_NTICKS-1);
% bbr = round(bb,YLIM_SIG_FIGS_I,'significant');
% dbb = abs(bbr-bb)*(YLIM_NTICKS-1);
% l = l-dbb/2;
% u = u+dbb/2;
% % l = dbb+1;
% yticks = unique(round(linspace(l,u,YLIM_NTICKS), ...
%     YLIM_SIG_FIGS,'significant'));
% y_lims = [min(yticks)-1e-8,max(yticks)+1e-8]; % sort([l-1e-8,u+1e-8]);
% yticks(yticks < 1e-8 & yticks > -1e-8) = 0;
% ytick_labs = cellstr(string(yticks));
% set(ax,'YLim',y_lims, ...
%     'YTick',yticks, ...
%     'YTickLabel',ytick_labs);
%--

set(ax,'YLim',params.y_lims, ...
    'YTick',params.yticks, ...
    'YTickLabel',params.ytick_labs);

%## X-LIMIT SETTING
xtick_labs = cellstr(string(XTICKS));
set(ax,'XLim',[min(XTICKS)-1e-6,max(XTICKS)+1e-6], ...
    'XTick',XTICKS, ...
    'XTickLabel',xtick_labs);

%## BAND LINES
% plot(PLOT_STRUCT.xlim,[0 0],'--','color','black');     
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
stat_ht_o = ht*STAT_HT_OFFSET_FACTOR+ht*SHADE_HT_FUDGE; %ax_ylim(1)*STAT_HT_OFFSET_FACTOR;
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
    % reg_en = unique([reg_en;length(reg)]);
    if length(reg_en)+1 == length(reg_st) %any(reg_en > length(reg))
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

%## PLOT GROUP STATS =====================================================
%{
% SHADE_BUFFER = 0;
% SIGGAP_SHRINK_SZ = 6;
% % SHADE_HT_FACTOR = 0.05;
% SHADE_FAC_ALPHA = 0.8;
% SHADE_EDG_ALPHA = 0.9;
% SHADE_EDG_COLOR = 'none';
% SHADE_FAC_COLOR = [0,0,0];
% STAT_HT_OFFSET_FACTOR = SHADE_HT_FACTOR; %SHADE_HT_FACTOR; %0.05;
%--
% ax_ylim = get(ax,'YLim');
% ht = ax_ylim(2)-ax_ylim(1);
% stat_ht_o = ht*STAT_HT_OFFSET_FACTOR; %ax_ylim(1)*STAT_HT_OFFSET_FACTOR;
% stat_ht = ht*SHADE_HT_FACTOR+stat_ht_o;

cnt = 1;
%--
fprintf('Plotting group comparison stats\n');
%-- plot over the ANOVA stats
stat_ht = soht_store(1,1);
stat_ht_o = soht_store(1,2);
for g_i = 1:length(alphag)
    %--
    % SHADE_FAC_COLOR = params.cmaps_scond(g_i,:);
    SHADE_FAC_COLOR = params.cmaps_stats(1,:);
    reg = alphag{g_i};
    %--
    regt = [0;reg];
    regt = regt == 1;
    reg_st = find(diff(regt) > 0);
    reg_en = find(diff(regt) < 0);
    reg_en = reg_en - 1;
    reg_en = unique([reg_en;length(reg)]);
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
                'DisplayName',sprintf('%s_s',params.xtick_label_c{i_i}));
            hold on;
            if r_i == 1
                paramso.stats_store{cnt} = Pa;
            end
            % if r_i == 1 
            %     paramso.leg_store{cntl} = Pa;
            % end
            % cntl = cntl + 1;
            cnt = cnt+1;            
        end
    end
    %--
    % stat_ht_o = stat_ht;
    % stat_ht = stat_ht_o + ht*SHADE_HT_FACTOR + ht*SHADE_HT_FUDGE;
end
%(04/16/2025) JS, this just seems to align with the ANOVA stats most of the
% time with little new information to expound on analysis
%}

%## PLOT COND STATS ======================================================
%{
% SHADE_BUFFER = 0;
% SIGGAP_SHRINK_SZ = 4;
% SHADE_HT_FACTOR = 0.05;
% SHADE_FAC_ALPHA = 0.8;
% SHADE_EDG_ALPHA = 0.9;
% SHADE_EDG_COLOR = 'none';
% SHADE_FAC_COLOR = [0,0,0];
% STAT_HT_OFFSET_FACTOR = SHADE_HT_FACTOR*2;
%--
% ax_ylim = get(ax,'YLim');
% ht = ax_ylim(2)-ax_ylim(1);
% stat_ht_o = ht*STAT_HT_OFFSET_FACTOR;
% stat_ht = ht*SHADE_HT_FACTOR+stat_ht_o;

%-- plot over the ANOVA stats
stat_ht = soht_store(2,1);
stat_ht_o = soht_store(2,2);
%--
cnt = 1;
fprintf('Plotting condition comparison stats\n');
for c_i = 1:length(alphac)
    % SHADE_FAC_COLOR = params.cmaps_sgroup(c_i,:);
    SHADE_FAC_COLOR = params.cmaps_stats(2,:);
    reg = alphac{c_i};
    %--
    regt = [0;reg];
    regt = regt == 1;
    reg_st = find(diff(regt) > 0);
    reg_en = find(diff(regt) < 0);
    reg_en = reg_en - 1;
    reg_en = unique([reg_en;length(reg)]);
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
                'DisplayName',sprintf('%s_s',params.xtick_label_g{i_i}));
            hold on;
            if r_i == 1
                paramso.stats_store{cnt} = Pa;
            end
            % if r_i == 1 
            %     paramso.leg_store{cntl} = Pa;
            % end
            % cntl = cntl + 1;
            cnt = cnt+1;            
        end
    end
    %--
    stat_ht_o = stat_ht;
    stat_ht = stat_ht_o + ht*SHADE_HT_FACTOR + ht*SHADE_HT_FUDGE;
end
%(04/16/2025) JS, this just seems to align with the ANOVA stats most of the
% time with little new information to expound on analysis
%}

%## OUTPUTS
paramso.stats_store = paramso.stats_store(~cellfun(@isempty,paramso.stats_store));
paramso.leg_store = paramso.leg_store(~cellfun(@isempty,paramso.leg_store));

%## LEGEND ============================================================== %%
if params.do_display_leg
    % params.legends_specs = {'FontName','Arial', ...
    %     'FontSize',8, ...
    %     'Orientation','Vertical', ...
    %     'Units','normalized'};
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

