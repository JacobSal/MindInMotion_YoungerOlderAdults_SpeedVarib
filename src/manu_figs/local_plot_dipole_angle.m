function [ax1,ax2,ax3] = local_plot_dipole_angle(fig,dip_fig_path, ...
    im_resize,paper_size,ax_position)
%PLOT_DIPOLE_SLICES Summary of this function goes here
%   Detailed explanation goes here
%## DIPOLE PLOTS =================================================== %%
hold on;

%## OPEN DIPOLE FIGURE & EXTRACT AXES
AX1_ADJUST = [1,0,-0.145];
%--
DIP_SIZE = 8;
DIP_LINEWIDTH = 0.01;
DIP_BORDER_CDATA = [0,0,0];
DIP_MARKERFACEALPHA = 0.9;
LIN_LINEWIDTH = 0.5;
LIN_LINECOLOR = [0,0,0,0.3];
LIN_LINESTYLE = '--';
%--
tmp = openfig(dip_fig_path);
try
    %-- scatter edits
    h = findobj(tmp.Children(end).Children,'type','Scatter');    
    for ii = 1:length(h)
        h(ii).CData = DIP_BORDER_CDATA;
        % h(ii).CData = h(ii).MarkerFaceColor;
        h(ii).LineWidth = DIP_LINEWIDTH;
        h(ii).MarkerFaceAlpha = DIP_MARKERFACEALPHA;
        % h(ii).MarkerFaceColor = [0,0,0];
        h(ii).SizeData = DIP_SIZE;
    end
    %-- line edits
    h = findobj(tmp.Children(end).Children,'type','Line');
    for ii = 1:length(h)
        h(ii).LineWidth = LIN_LINEWIDTH;
        h(ii).LineStyle = LIN_LINESTYLE;
        h(ii).Color = LIN_LINECOLOR;
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
ax1_x = pp(3)*im_resize*AX1_ADJ/paper_size(3);
ax1_y = pp(4)*im_resize*AX1_ADJ/paper_size(3);
set(ax1,'OuterPosition',[0,0,1,1], ...
    'Position',[ax_position(1)+ax1_x*AX1_X_O,ax_position(2)+ax1_y*AX1_Y_O,ax1_x,ax1_y]);
%-- close copied fig
close(tmp);
end

