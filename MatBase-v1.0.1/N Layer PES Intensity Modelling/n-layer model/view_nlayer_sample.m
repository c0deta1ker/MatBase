function fig = view_nlayer_sample(lyr_mat, lyr_thick, lyr_cols)
% fig = view_nlayer_sample(lyr_mat, lyr_thick, lyr_cols)
%   This function plots the solutions to the 'nlayer_pes_model()' function;
%   the n-layered sample stack and the corresponding photoelectron
%   contribution from each one of the layers. The user can input the colors
%   for each one of the layers using the 'lyr_cols' argument; this is
%   convenient when you want to color-match the schematic to your own
%   drawings.
%
%   IN:
%   -   lyr_mat:        Mx1 cell-vector of the material for each layer in the stack; e.g. "Si", "SiO2", "Al2O3"...
%   -   lyr_thick:      Mx1 cell-vector of the thickness of each layer in the stack (in nano-metres)
%   -   lyr_cols:       Mx1 cell-vector of the [R,G,B] color of each independent layer.
%
%   OUT:
%   -   fig:            MATLAB figure object with the ARPES data plotted.

%% Default parameters
def_cols = flipud(num2cell([0.25,0.25,0.25;lines(length(lyr_mat)-1)], 2));
% -- Defining the default parameters
if nargin < 3; lyr_cols = def_cols; end
if isempty(lyr_cols); lyr_cols = def_cols; end
% -- Extracting the total number of layers to be probed
Nlyrs       = length(lyr_mat);

%% 1 - Plotting the the material stack
fig= figure(); 
fig.Position(3) = 400;
fig.Position(4) = fig.Position(3);
hold on;
x_width	= 5;
bulk_thickness = 5;
if Nlyrs == 1
    y_cum   = cell2mat(lyr_thick);
else
    y_cum   = cumsum(cell2mat(lyr_thick));
    y_cum(isinf(y_cum)) = y_cum(end-1) + bulk_thickness;
end
% -- Plotting each layer from bottom-up
for i = Nlyrs:-1:1
    patch([-1, -1, 1, 1, -1].*x_width, [0, y_cum(i), y_cum(i), 0, 0].*-1,...
        lyr_cols{i}, 'edgecolor', [0 0 0]);
end
% -- Adding text for each material type
for i = 1:Nlyrs
    if i == 1;  y_loc = 0 - 0.5*(y_cum(i) - 0);
    else;       y_loc = -y_cum(i-1) - 0.5*(y_cum(i) - y_cum(i-1));
    end
    if i == Nlyrs
        text(0, y_loc, "(Bulk) "+string(lyr_mat{i}),...
        'color', 'k', 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'FontWeight','bold', 'FontSize',10);
    else
        text(0, y_loc, sprintf("(Layer %i) ",i)+string(lyr_mat{i}),...
        'color', 'k', 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'FontWeight','bold', 'FontSize',10); 
    end
end
% -- Box Styling and Axis properties
ylabel('Depth From Surface [nm]', 'fontweight', 'bold');
title('Sample model stack');
axis([-x_width, x_width, -1*max(y_cum(:)), 0]);
% - Formatting the axes
ax = gca;
% Font properties
ax.FontName         = 'Segoe UI'; 
ax.FontWeight       = 'normal'; 
ax.FontSize         = 10;
% Tick properties
ax.XMinorTick       = 'off'; 
ax.YMinorTick       = 'off';
ax.TickDir          = 'both';
ax.XColor           = [0 0 0]; 
ax.YColor           = [0 0 0];
% Ruler properties
ax.XAxisLocation    = 'bottom';             % 'bottom' | 'top' | 'origin'
ax.YAxisLocation    = 'left';               % 'left' | 'right' | 'origin'
% Box Styling properties
ax.Color            = [1 1 1];
ax.Box              = 'off';                % 'on' | 'off'
ax.LineWidth        = 0.75;     
ax.Layer            = 'Top';

end