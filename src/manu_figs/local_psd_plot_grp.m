function [paramso] = local_psd_plot_grp(ax,psd_dat_in,params,PLOT_STRUCT,LINE_STRUCT)
%LOCAL_PSD_PLOT Summary of this function goes here
%   Detailed explanation goes here
%% EXTRACT PSD DATA ==================================================== %%
paramso.stats_store = cell(4,1); 
paramso.ax_store = cell(4,1);
paramso.leg_store = cell(4,1);

%## FUNC PARAMETER SETS
XFREQ_LINES = [3,8,13,30];
STAT_ALPHA = 0.05;
% AXES_DEFAULT_PROPS = {'box','off', ...
%     'xtick',[], ...
%     'ytick',[],...
%     'ztick',[], ...
%     'xcolor',[1,1,1], ...
%     'ycolor',[1,1,1]};
%--
csz = size(psd_dat_in,1);
gsz = size(psd_dat_in,2);
%## STATS
pcond = {};
pgroup = {};
pinter = {};
%--
[pcond, pgroup, pinter, ~, ~, ~] = ...
    std_stat(psd_dat_in, params.stats);
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
    for c_i = 1:csz
        psd_in_c(:,c_i,g_i) = mean(psd_dat_in{c_i,g_i},2);
        % err_bnd_in(:,c_i,g_i,1) = mean(psd_dat_in{c_i,g_i},2) + std(psd_dat_in{c_i,g_i},[],2);
        % err_bnd_in(:,c_i,g_i,2) = mean(psd_dat_in{c_i,g_i},2) - std(psd_dat_in{c_i,g_i},[],2);
        err_bnd_in(:,c_i,g_i,1) = mean(psd_dat_in{c_i,g_i},2) + std(psd_dat_in{c_i,g_i},[],2)/(sqrt(size(psd_dat_in{c_i,g_i},2)));
        err_bnd_in(:,c_i,g_i,2) = mean(psd_dat_in{c_i,g_i},2) - std(psd_dat_in{c_i,g_i},[],2)/(sqrt(size(psd_dat_in{c_i,g_i},2)));
    end
end
% psd_meanc = squeeze(mean(psd_in_c,2));

%## PLOTTING PARAMETERS
for g_i = 1:gsz
    for c_i = 1:csz
        psd_meanc = psd_in_c(:,c_i,g_i);
        % LINE_STRUCT.line_color = params.cmaps(g_i,:);
        % LINE_STRUCT.line_label = params.xtick_label{g_i};
        % LINE_STRUCT.err_color = params.cmaps(g_i,:);
        LINE_STRUCT.line_color = params.cmaps(c_i,:);
        LINE_STRUCT.line_label = params.xtick_label{c_i};
        LINE_STRUCT.err_color = params.cmaps(c_i,:);
        LINE_STRUCT.err_bnd_vec = [squeeze(err_bnd_in(:,c_i,g_i,1)), ...
                squeeze(err_bnd_in(:,c_i,g_i,2))];
        %## PLOT
        [~,~,Li] = plot_psd(ax,psd_meanc,params.freqs, ...
            'LINE_STRUCT',LINE_STRUCT, ...
            'PLOT_STRUCT',PLOT_STRUCT);
        % paramso.ax_store{c_i,g_i} = ax;
        paramso.leg_store{c_i,g_i} = Li;
        %--
        hold on;
    end
end
plot(PLOT_STRUCT.xlim,[0 0],'--','color','black');     
for xx = 1:length(XFREQ_LINES)
    xline(XFREQ_LINES(xx),'--');
end

% [axsignif,Pa] = plot_psd_stats(ax,params.freqs,alphac{1}, ...
%     'background','Frequency (Hz)');


%## PLOT GROUP STATS
% SHADE_BUFFER = 3;
% SHADE_HT_FACTOR = 0.05;
% SHADE_FAC_ALPHA = 0.8;
% SHADE_EDG_ALPHA = 0.9;
% SHADE_EDG_COLOR = 'none';
% SHADE_FAC_COLOR = [0,0,0];
% STAT_HT_OFFSET_FACTOR = -0.05;
% 
% %--
% ax_ylim = get(ax,'YLim');
% stat_ht_o = ax_ylim(2)*SHADE_HT_FACTOR; %ax_ylim(1)*STAT_HT_OFFSET_FACTOR;
% stat_ht = ax_ylim(2)*SHADE_HT_FACTOR;
% cnt = 1;
% %--
% for g_i = 1:length(alphag)
%     reg = alphac{g_i};
%     regt = [0;reg];
%     regt = regt == 1;
%     reg_st = find(diff(regt) > 0);
%     reg_en = find(diff(regt) < 0);
%     reg_en = reg_en - 1;
%     reg_en = unique([reg_en;length(reg)]);
%     %-- shading buffing to ensure visualization
%     reg_st = floor(reg_st - SHADE_BUFFER/2);
%     reg_st(reg_st < 1) = 1;
%     reg_en = ceil(reg_en + SHADE_BUFFER/2);
%     reg_en(reg_en > length(reg)) = length(reg);
%     %-- apply stats if true
%     if ~isempty(reg_st)
%         for r_i = 1:length(reg_st)            
%             stx = params.freqs(reg_st(r_i));
%             enx = params.freqs(reg_en(r_i));
%             %-- plot patch
%             Pa = patch(ax,[stx,enx,enx,stx], ...
%                 [ax_ylim(1)+stat_ht_o,ax_ylim(1)+stat_ht_o,ax_ylim(1)+stat_ht,ax_ylim(1)+stat_ht], ...
%                 SHADE_FAC_COLOR); 
%             hold on;
%             set(Pa,'edgecolor',SHADE_EDG_COLOR, ...
%                 'facealpha',SHADE_FAC_ALPHA, ...
%                 'edgealpha',SHADE_EDG_ALPHA);
%             hold on;
%             paramso.stats_store{cnt} = Pa;
%             cnt = cnt+1;
%         end
%     end
% end

%## PLOT COND STATS
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

%## LEGEND
%- lg2
% leg_store = [leg_store{:}];
leg_store = [paramso.leg_store{:,1}];
legend(gca,leg_store);
[lg2,~,~,~]  = legend('boxoff');
tmp = get(lg2,'String');
% set(lg2,'String',tmp, ...
%     'FontName',PLOT_STRUCT.font_name, ...
%     'FontSize',params.leg_txt_size, ...
%     'Orientation','horizontal', ...
%     'Units','normalized', ...
%     'Position',[params.leg_position(1),...
%         params.leg_position(2), ...
%         lg2.Position(3), ...
%         lg2.Position(4)]);
set(lg2,'String',tmp, ...
    'FontName',PLOT_STRUCT.font_name, ...
    'FontSize',params.leg_txt_size, ...
    'Orientation','vertical', ...
    'Units','normalized', ...
    'Position',[lg2.Position(1)+params.leg_position(1),...
        lg2.Position(2)+params.leg_position(2), ...
        lg2.Position(3), ...
        lg2.Position(4)]);
lg2.ItemTokenSize(1) = params.leg_token_size;
hold off;

end

