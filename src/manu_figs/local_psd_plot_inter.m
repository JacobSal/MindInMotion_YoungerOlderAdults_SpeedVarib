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
csz = size(psd_dat_in,1);
gsz = size(psd_dat_in,2);
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
            case 'group'
                % tmp = cat(2,psd_dat_in{:,g_i}); % all conds & subjs
                tmp = mean(cat(3,psd_dat_in{:,g_i}),3); % mean cond for subjs
            case 'cond'
                tmp = cat(2,psd_dat_in{c_i,:});
        end
        psd_in_c(:,c_i,g_i) = mean(tmp,2);
        % N = sqrt(size(tmp,2));
        N = sqrt(num_subj);
        err_bnd_in(:,c_i,g_i,1) = mean(tmp,2) + std(tmp,[],2)/(N);
        err_bnd_in(:,c_i,g_i,2) = mean(tmp,2) - std(tmp,[],2)/(N);
    end
end
% psd_meanc = squeeze(mean(psd_in_c,2));

%## PLOTTING PARAMETERS
switch params.line_plot_opt
    case 'group'
        csz = 1;
    case 'cond'                
        gsz = 1;       
end
% c_i = 1;
for g_i = 1:gsz
    for c_i = 1:csz
        psd_meanc = squeeze(psd_in_c(:,c_i,g_i));
        switch params.line_plot_opt
            case 'inter'
            case 'group'
                LINE_STRUCT.line_color = params.cmaps(g_i,:);
                LINE_STRUCT.line_label = params.xtick_label{g_i};
                LINE_STRUCT.err_color = params.cmaps(g_i,:);
            case 'cond'                
                LINE_STRUCT.line_color = params.cmaps(c_i,:);
                LINE_STRUCT.line_label = params.xtick_label{c_i};
                LINE_STRUCT.err_color = params.cmaps(c_i,:);           
        end
        %--
        LINE_STRUCT.err_bnd_vec = [squeeze(err_bnd_in(:,c_i,g_i,1)), ...
                    squeeze(err_bnd_in(:,c_i,g_i,2))];
        %## PLOT
        [~,~,Li] = plot_psd(ax,psd_meanc,params.freqs, ...
            'LINE_STRUCT',LINE_STRUCT, ...
            'PLOT_STRUCT',PLOT_STRUCT);
        % paramso.ax_store{c_i,g_i} = ax;
        paramso.leg_store{cntl} = Li;
        cntl = cntl + 1;
        %--
        hold on;
    end
end
set(ax,'FontWeight',AX_FONTWEIGHT)

%## Y-LIMIT SETTING
YLIM_NTICKS = 5;
YLIM_SIG_FIGS_I = 1;
YLIM_SIG_FIGS = 2;
YLIM_FAC = 0.15;
%-- ylim
y_lims = get(ax,'YLim');
dy = YLIM_FAC*(y_lims(2)-y_lims(1));
%-- change 0 point for ylim if exists
[vy,iy] = min(y_lims);
if vy == 0 && iy == 1
    y_lims(iy) = -1e-8;
elseif vy == 0 && iy == 2
    y_lims(iy) = 1e-8;
end
u = round(y_lims(2),YLIM_SIG_FIGS_I,'significant');
l = round(y_lims(1),YLIM_SIG_FIGS_I,'significant');
% u = round((ax_ylim(2)+dy),YLIM_SIG_FIGS_I,'significant');
% l = round((ax_ylim(1)-dy),YLIM_SIG_FIGS_I,'significant'); 
bb = (u-l)/(YLIM_NTICKS-1);
bbr = round(bb,YLIM_SIG_FIGS_I,'significant');
dbb = abs(bbr-bb)*(YLIM_NTICKS-1);
l = l-dbb/2;
u = u+dbb/2;
% l = dbb+1;
yticks = unique(round(linspace(l,u,YLIM_NTICKS), ...
    YLIM_SIG_FIGS,'significant'));
y_lims = [min(yticks)-1e-8,max(yticks)+1e-8]; % sort([l-1e-8,u+1e-8]);
yticks(yticks < 1e-8 & yticks > -1e-8) = 0;
ytick_labs = cellstr(string(yticks));
set(ax,'YLim',y_lims, ...
    'YTick',yticks, ...
    'YTickLabel',ytick_labs);

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
stats_titles = {'cond','group','inter'};
%--
y_lims = get(ax,'YLim');
ht = y_lims(2)-y_lims(1);
stat_ht_o = ht*STAT_HT_OFFSET_FACTOR+ht*SHADE_HT_FUDGE; %ax_ylim(1)*STAT_HT_OFFSET_FACTOR;
stat_ht = ht*SHADE_HT_FACTOR+stat_ht_o;
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
            paramso.stats_store{cnt} = Pa;
            if r_i == 1 
                paramso.leg_store{cntl} = Pa;
            end
            cntl = cntl + 1;
            cnt = cnt + 1;          
        end        
    end
    %--
    stat_ht_o = stat_ht;
    stat_ht = stat_ht_o + ht*SHADE_HT_FACTOR + ht*SHADE_HT_FUDGE;
end
%(04/02/2025) JS, so far the leat buggy stats algorithm.

%## PLOT GROUP STATS
%{
SHADE_BUFFER = 0;
SIGGAP_SHRINK_SZ = 6;
SHADE_HT_FACTOR = 0.05;
SHADE_FAC_ALPHA = 0.8;
SHADE_EDG_ALPHA = 0.9;
SHADE_EDG_COLOR = 'none';
% SHADE_FAC_COLOR = [0,0,0];
STAT_HT_OFFSET_FACTOR = SHADE_HT_FACTOR; %SHADE_HT_FACTOR; %0.05;
%--
ax_ylim = get(ax,'YLim');
ht = ax_ylim(2)-ax_ylim(1);
stat_ht_o = ht*STAT_HT_OFFSET_FACTOR; %ax_ylim(1)*STAT_HT_OFFSET_FACTOR;
stat_ht = ht*SHADE_HT_FACTOR+stat_ht_o;
cnt = 1;
%--
fprintf('Plotting group comparison stats\n');
for g_i = 1:length(alphag)
    %--
    SHADE_FAC_COLOR = params.cmaps_scond(g_i,:);
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
                [ax_ylim(1)+stat_ht_o,ax_ylim(1)+stat_ht_o,ax_ylim(1)+stat_ht,ax_ylim(1)+stat_ht], ...
                SHADE_FAC_COLOR);
            hold on;
            set(Pa,'edgecolor',SHADE_EDG_COLOR, ...
                'facealpha',SHADE_FAC_ALPHA, ...
                'edgealpha',SHADE_EDG_ALPHA);
            hold on;
            paramso.stats_store{cnt} = Pa;
            cnt = cnt+1;            
        end
    end
end
%}

%## PLOT COND STATS
%{
SHADE_BUFFER = 0;
SIGGAP_SHRINK_SZ = 4;
SHADE_HT_FACTOR = 0.05;
SHADE_FAC_ALPHA = 0.8;
SHADE_EDG_ALPHA = 0.9;
SHADE_EDG_COLOR = 'none';
% SHADE_FAC_COLOR = [0,0,0];
STAT_HT_OFFSET_FACTOR = SHADE_HT_FACTOR*2;
%--
ax_ylim = get(ax,'YLim');
ht = ax_ylim(2)-ax_ylim(1);
stat_ht_o = ht*STAT_HT_OFFSET_FACTOR;
stat_ht = ht*SHADE_HT_FACTOR+stat_ht_o;
cnt = 1;
%--
fprintf('Plotting condition comparison stats\n');
for c_i = 1:length(alphac)
    SHADE_FAC_COLOR = params.cmaps_sgroup(c_i,:);
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
                [ax_ylim(1)+stat_ht_o,ax_ylim(1)+stat_ht_o,ax_ylim(1)+stat_ht,ax_ylim(1)+stat_ht], ...
                SHADE_FAC_COLOR);
            hold on;
            set(Pa,'edgecolor',SHADE_EDG_COLOR, ...
                'facealpha',SHADE_FAC_ALPHA, ...
                'edgealpha',SHADE_EDG_ALPHA);
            hold on;
            paramso.stats_store{cnt} = Pa;
            cnt = cnt+1;            
        end
    end
end
%}
%% LEGEND ============================================================== %%
if params.do_display_leg
    % params.legends_specs = {'FontName','Arial', ...
    %     'FontSize',8, ...
    %     'Orientation','Vertical', ...
    %     'Units','normalized'};
    %-- legend
    leg_store = [paramso.leg_store{:}];
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

