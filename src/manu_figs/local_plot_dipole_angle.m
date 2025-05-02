function [ax1] = local_plot_dipole_angle(fig,dip_fig_path, ...
    PLOT_STRUCT)
%PLOT_DIPOLE_SLICES Summary of this function goes here
%   Detailed explanation goes here
%## DIPOLE PLOTS =================================================== %%
hold on;
tps = PLOT_STRUCT;
%## OPEN DIPOLE FIGURE & EXTRACT AXES
AX1_ADJUST = [1,0,-0.145];
%--
tmp = openfig(dip_fig_path);
try
    %-- scatter edits
    h = findobj(tmp.Children(end).Children,'type','Scatter');    
    for ii = 1:length(h)
        h(ii).CData = tps.dip_border_col;
        % h(ii).CData = h(ii).MarkerFaceColor;
        h(ii).LineWidth = tps.dip_line_width;
        h(ii).MarkerFaceAlpha = tps.dip_alpha;
        % h(ii).MarkerFaceColor = [0,0,0];
        h(ii).SizeData = tps.dip_sz;
    end
    %-- line edits
    h = findobj(tmp.Children(end).Children,'type','Line');
    for ii = 1:length(h)
        h(ii).LineWidth = tps.lin_line_width;
        h(ii).LineStyle = tps.lin_line_style;
        h(ii).Color = tps.lin_line_col;
    end    
catch
    fprintf('not scatter 3d');
end
ax = get(tmp,'CurrentAxes');
daspect = get(ax,'DataAspectRatio');
ac = get(ax,'Children');

%## AXES 1
hold on;
ax1 = axes('Parent',fig, ...
    'DataAspectRatio',daspect, ...
    'Units','normalized');
copyobj(ac,ax1);
set(ax1,'box','off', ...
    'xtick',[], ...
    'ytick',[],...
    'ztick',[], ...
    'xcolor',[1,1,1], ...
    'ycolor',[1,1,1], ...
    'zcolor',[1,1,1]);
view(ax1,[45,30])
drawnow;
pp = get(ax1,'Position');
%-
AX1_ADJ = AX1_ADJUST(1);
AX1_Y_O = AX1_ADJUST(3);
AX1_X_O = AX1_ADJUST(2);
ax1_x = pp(3)*tps.im_resize*AX1_ADJ/tps.pg_size(3);
ax1_y = pp(4)*tps.im_resize*AX1_ADJ/tps.pg_size(4);
set(ax1,'OuterPosition',[0,0,1,1], ...
    'Position',[tps.position(1)+ax1_x*AX1_X_O,tps.position(2)+ax1_y*AX1_Y_O,ax1_x,ax1_y]);
%-- close copied fig
close(tmp);
end

