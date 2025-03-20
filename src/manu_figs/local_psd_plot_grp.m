function [paramso] = local_psd_plot(ax,psd_dat_in,params,PLOT_STRUCT,LINE_STRUCT)
%LOCAL_PSD_PLOT Summary of this function goes here
%   Detailed explanation goes here
%## EXTRACT PSD DATA =============================================== %%
%%
% params.im_resize = 0.8;
% params.ax_w = 0.3;
% params.ax_h = 0.25;
% params.ax_init_x = 0.09;
% params.ax_init_y = 0.7;    
% params.designs = {{'flat','low','med','high'},{'0p25','0p5','0p75','1p0'}};
% params.d_i = 1;
% params.des_i = 2;
% params.g_i = 1;
% params.cmaps = [];
% params.xtick_label = [];
% params.title;
% params.cl_i = 1;
% params.stats = stats;
% params.freqs = fooof_freqs;
%--
paramso.ax_store = cell(4,1);
paramso.leg_store = cell(4,1);
paramso.psd_char = [];
paramso.y_lim_store = zeros(4,2);
%--
% params.x_shift = params.ax_init_x;
% params.y_shift = params.ax_init_y;

%## FUNC PARAMETER SETS
csz = size(psd_dat_in,1);
gsz = size(psd_dat_in,2);
STAT_ALPHA = 0.05;
AXES_DEFAULT_PROPS = {'box','off', ...
    'xtick',[], ...
    'ytick',[],...
    'ztick',[], ...
    'xcolor',[1,1,1], ...
    'ycolor',[1,1,1]};

%## STATS
params.stats.groupstats = 'off';
[pcond, pgroup, pinter, ~, ~, ~] = ...
    std_stat(psd_dat_in, params.stats);
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
        err_bnd_in(:,c_i,g_i,1) = mean(psd_dat_in{c_i,g_i},2) + std(psd_dat_in{c_i,g_i},[],2);
        err_bnd_in(:,c_i,g_i,2) = mean(psd_dat_in{c_i,g_i},2) - std(psd_dat_in{c_i,g_i},[],2);
    end
end
psd_meanc = squeeze(mean(psd_in_c,2));

%## PLOTTING PARAMETERS
PLOT_STRUCT.title = params.title; %{sprintf('%s',cluster_titles{cluster_inds_plot(params.cl_i)})};
PLOT_STRUCT.ax_position = [params.x_shift, ...
    params.y_shift, ...
    params.ax_w*params.im_resize, ...
    params.ax_h*params.im_resize];
PLOT_STRUCT.xlim = [3,40];
%-- ylim calc
mu = mean(cat(2,psd_dat_in{:}),[2,1]);
sd = std(cat(2,psd_dat_in{:}),[],[2,1]);
paramso.y_lim_store= [mu-1.75*sd,mu+1.75*sd];
PLOT_STRUCT.ylim = paramso.y_lim_store;
%--
LINE_STRUCT.line_avg_fcn = @(x) mean(x,2);
LINE_STRUCT.line_color = params.cmaps(params.g_i,:);
LINE_STRUCT.line_alpha = 0.7;
LINE_STRUCT.line_label = params.xtick_label{params.g_i};
%--
LINE_STRUCT.do_err_shading = true;
LINE_STRUCT.err_color = params.cmaps(params.g_i,:);
LINE_STRUCT.err_alpha = 0.3;
LINE_STRUCT.err_bnd_vec = [mean(psd_in_c,2) + std(psd_in_c,[],2), ...
        mean(psd_in_c,2) - std(psd_in_c,[],2)];
%## PLOT
[ax,Pa,Li] = plot_psd(ax,psd_meanc,params.freqs, ...
    'LINE_STRUCT',LINE_STRUCT, ...
    'PLOT_STRUCT',PLOT_STRUCT);
paramso.ax_store{1} = ax;
paramso.leg_store{1} = Li;
%--
hold on;

% [axsignif,Pa] = plot_psd_stats(ax,params.freqs,alphac{1}, ...
%     'background','Frequency (Hz)');

%## PLOT GROUP STATS
%{
stat_color = [0,0,0];
freqs = params.freqs;
%--
axx = axes();
pos_in = ax.Position;
pos_in(4) = pos_in(4)*0.1;
set(axx,AXES_DEFAULT_PROPS{:}, ...
    'Position',pos_in);
ax_ylim = get(axx,'YLim');
stat_ht = ax_ylim(2)*0.1;
%--
for c_i = 1:csz
    reg = alphag{c_i};
    regt = [0;reg];
    regt = regt == 1;
    reg_st = find(diff(regt) > 0);
    reg_en = find(diff(regt) < 0);
    reg_en = reg_en - 1;
    reg_en = unique([reg_en;length(reg)]);
    if ~isempty(reg_st)
        for r_i = 1:length(reg_en)            
            stx = freqs(reg_st(r_i));
            enx = freqs(reg_en(r_i));
            Pa = patch(axx,[stx,enx,enx,stx], ...
                [ax_ylim(1),ax_ylim(1),ax_ylim(1)+stat_ht,ax_ylim(1)+stat_ht], ...
                stat_color); 
            hold on;
            set(Pa,'edgecolor','none', ...
                'facealpha',0.5, ...
                'edgealpha',0.2);
            hold on;
        end
    end
    %## PLOT    
end
%}

%## PLOT COND STATS
SHADE_BUFFER = 3;
SHADE_HT_FACTOR = 0.05;
SHADE_FAC_ALPHA = 0.8;
SHADE_EDG_ALPHA = 0.9;
SHADE_EDG_COLOR = 'none';
SHADE_FAC_COLOR = [0,0,0];
%--
ax_ylim = get(ax,'YLim');
stat_ht = ax_ylim(2)*SHADE_HT_FACTOR;
%--
for c_i = 1:length(alphac)
    reg = alphac{c_i};
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
    if ~isempty(reg_st)
        for r_i = 1:length(reg_st)            
            stx = params.freqs(reg_st(r_i));
            enx = params.freqs(reg_en(r_i));
            %-- plot patch
            Pa = patch(ax,[stx,enx,enx,stx], ...
                [ax_ylim(1),ax_ylim(1),ax_ylim(1)+stat_ht,ax_ylim(1)+stat_ht], ...
                SHADE_FAC_COLOR); 
            hold on;
            set(Pa,'edgecolor',SHADE_EDG_COLOR, ...
                'facealpha',SHADE_FAC_ALPHA, ...
                'edgealpha',SHADE_EDG_ALPHA);
            hold on;
        end
    end
end

end

