function [] = makeFigureNice( OPTS )
% simple function to make figure look nicer.
% 
% Inputs:
%     OPTS
%       type: struct
%       desc: OPTIONAL options
%             'grid_line_style': 
%                   any acceptable matlab grid line style,
%                   such as '--' for dashed or '-' for solid
%             'fig_size':
%                       Choose dimensions of figure [x size, y size]
%             'fig_pos':
%                       Choose position of figure [x pos, y pos]
%             'font_size':
%                   font size to use for all fonts
%             'line_width'
%             'axis_line_width'
%                   sets line width for ticks and grid
%             'fig_units':
%                   sets figure dimension units
%             'is_yyaxis':
%                   

% DEFAULTS
default_fig_units       = 'pixels';
default_fig_size        = [800, 600];
default_fig_pos         = [100, 100];
default_grid_line_style = '--';
default_font_size       = 16;
default_line_width      = 2.0;
default_axis_line_width = 1.5;

if nargin == 0
   OPTS.fig_units       = default_fig_units;
   OPTS.fig_size        = default_fig_size;
   OPTS.fig_pos         = default_fig_pos;
   OPTS.grid_line_style = default_grid_line_style;  % default to dashed lines
   OPTS.font_size       = default_font_size;  
   OPTS.line_width      = default_line_width; 
   OPTS.axis_line_width = default_axis_line_width;
else
    
   if ~isfield( OPTS, 'fig_units' )
       OPTS.fig_units = default_fig_units;
   end 
    
   if ~isfield( OPTS, 'fig_size' )
       OPTS.fig_size = default_fig_size;
   end
   
   if ~isfield( OPTS, 'fig_pos' )
       OPTS.fig_pos = default_fig_pos;
   end
   
   if ~isfield( OPTS, 'grid_line_style' )
       OPTS.grid_line_style = default_grid_line_style;
   end
   
   if ~isfield( OPTS, 'font_size' )
       OPTS.font_size = default_font_size;
   end
   
   if ~isfield( OPTS, 'line_width' )
       OPTS.line_width = default_line_width;
   end
   
   if ~isfield( OPTS, 'axis_line_width' )
       OPTS.axis_line_width = default_axis_line_width;
   end
    
end

% Size the figure
fig = gcf;
set( fig, 'Units', OPTS.fig_units );
set(fig,'Position',[OPTS.fig_pos(1), OPTS.fig_pos(2), OPTS.fig_size(1), OPTS.fig_size(2)]);

% Set the fonts
set(gca,'FontSize',OPTS.font_size);
set(get(gca,'YLabel'),'FontSize',OPTS.font_size);
set(get(gca,'XLabel'),'FontSize',OPTS.font_size);
set(get(gca,'Title'),'FontSize',OPTS.font_size);

% Make the lines thicker
hPlot = get(get(fig,'Children'),'Children');

try     % won't execute if there are no lines to thicken
    if iscell(hPlot)
        set(hPlot{1},'LineWidth',OPTS.line_width);
        set(hPlot{2},'LineWidth',OPTS.line_width);
    else
        set(hPlot,'LineWidth',OPTS.line_width);
    end
catch
    % do nothing
end

% turn the grid on
grid on;

% make tick marks thicker
ax = gca;
ax.TickLength       = [0.01, 0.01]; %[0.03, 0.075];
ax.LineWidth        = OPTS.axis_line_width;
ax.GridLineStyle    = OPTS.grid_line_style;

% make axes black
ax.XColor = [0,0,0];
if length(ax.YAxis) < 2 % don't re-color y axis if its 2 yaxes plot
    ax.YColor = [0,0,0];
end

% % Copy to clipboard
% print -dbitmap -r200;

end