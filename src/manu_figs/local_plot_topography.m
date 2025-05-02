function [ax1] = local_plot_topography(fig,STUDY,cl_n,PLOT_STRUCT) %groups,g_chars,g_chars_topo, ...
    % ax_position,im_resize)
%PLOT_TOPOGRAPHY Summary of this function goes here
%   Detailed explanation goes here
%##
hold on;
%## PARAMS
AXES_DEFAULT_PROPS = {'box','off', ...
    'xtick',[], ...
    'ytick',[], ...
    'ztick',[], ...
    'xcolor',[1,1,1], ...
    'ycolor',[1,1,1]};
%--
tmp_plot_struct = PLOT_STRUCT;

%## TOPO PLOT
[~,h] = std_topoplot_CL(STUDY,cl_n,'together');
set(h,'color','w')
set(h,'Units','normalized');    
ax = get(h,'CurrentAxes');
colormap(h,linspecer);   
% pos = get(ax,'Position');
% opos = get(ax,'OuterPosition');
ac = get(ax,'Children');
ax1 = axes('Parent',fig,...
    'DataAspectRatio',get(ax,'DataAspectRatio'));
copyobj(ac,ax1);
set(ax1,AXES_DEFAULT_PROPS{:});
colormap(ax1,linspecer)
ax1.Title.String = tmp_plot_struct.title; %sprintf('N=%i',length(STUDY.cluster(cl_i).sets));
ax1.Title.Interpreter = 'none';
ax1.Title.FontSize = tmp_plot_struct.font_size;
ax1.FontWeight = tmp_plot_struct.font_weight;
%-- axes edits
ax1.FontSize = tmp_plot_struct.font_size; %PLOT_STRUCT.font_size;
ax1.FontName = tmp_plot_struct.font_name;
ax1.FontWeight = tmp_plot_struct.font_weight;
ax1.OuterPosition = [0,0,1,1];
ax1.Units = 'Normalized';
pp = get(ax1,'Position');    
ax1.Position = [tmp_plot_struct.position(1),tmp_plot_struct.position(2), ...
    tmp_plot_struct.position(3)+pp(3)*tmp_plot_struct.im_resize,tmp_plot_struct.position(4)+pp(4)*tmp_plot_struct.im_resize];  %[left bottom width height]
close(h)
end

