function fig = view_nlayer_pes_model_vs_theta(pes_model, lyr_cols)
% fig = view_nlayer_pes_model_vs_theta(pes_model, lyr_cols)
%   This function plots the solutions to the 'nlayer_pes_model()' function;
%   the n-layered sample stack and the corresponding photoelectron
%   contribution from each one of the layers. The user can input the colors
%   for each one of the layers using the 'lyr_cols' argument; this is
%   convenient when you want to color-match the schematic to your own
%   drawings.
%
%   IN:
%   -   pes_model:      data structure that contains all the pes model parameters and variables (from 'nlayer_pes_model()').
%   -   lyr_cols:       Mx1 cell-vector of the [R,G,B] color of each independent layer.
%
%   OUT:
%   -   fig:            MATLAB figure object with the ARPES data plotted.

%% Default parameters
def_cols = flipud(num2cell([0.25,0.25,0.25;lines(length(pes_model.lyr_mat)-1)], 2));
% -- Defining the default parameters
if nargin < 2; lyr_cols = def_cols; end
if isempty(lyr_cols); lyr_cols = def_cols; end

%% 1 - Plotting the the model solutions
fig = figure();
fig.Position(3) = 900; 
fig.Position(4) = 450;
%% 1.1 - Plotting the intensity profiles for each layer
subplot(1,3,1); hold on;
x_width	= 5;
bulk_thickness = 5;
if pes_model.Nlyrs == 1
    y_cum   = cell2mat(pes_model.lyr_thick);
else
    y_cum   = cumsum(cell2mat(pes_model.lyr_thick));
    y_cum(isinf(y_cum)) = y_cum(end-1) + bulk_thickness;
end
% -- Plotting each layer from bottom-up
for i = pes_model.Nlyrs:-1:1
    patch([-1, -1, 1, 1, -1].*x_width, [0, y_cum(i), y_cum(i), 0, 0].*-1,...
        lyr_cols{i}, 'edgecolor', [0 0 0]);
end
% -- Adding text for each material type
for i = 1:pes_model.Nlyrs
    if i == 1;  y_loc = 0 - 0.5*(y_cum(i) - 0);
    else;       y_loc = -y_cum(i-1) - 0.5*(y_cum(i) - y_cum(i-1));
    end
    if i == pes_model.Nlyrs
        text(0, y_loc, "(Bulk) "+string(pes_model.lyr_mat{i}),...
        'color', 'k', 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'FontWeight','bold', 'FontSize',10);
    else
        text(0, y_loc, sprintf("(Layer %i) ",i)+string(pes_model.lyr_mat{i}),...
        'color', 'k', 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'FontWeight','bold', 'FontSize',10);
    end
end
% - Box Styling and Axis properties
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
%% 4.2 - Plotting the intensity profiles for each layer
subplot(1,3,[2,3]); hold on;
% - Filing through each layer
for i = 1:pes_model.Nlyrs
    % -- If only one photon energy is given
    if length(pes_model.hv) == 1
        plot(pes_model.theta, pes_model.lyr_ints0{i}, 'k.-', 'color', lyr_cols{i}, 'markerfacecolor', lyr_cols{i}, 'linewidth', 2.5);
    % -- If a range of angles are to be plotted
    else
        num_cols    = length(pes_model.hv);
        start_col   = lyr_cols{i};
        end_col     = lyr_cols{i} + 0.50*(1-lyr_cols{i});
        colors_p = [linspace(start_col(1),end_col(1),num_cols)', linspace(start_col(2),end_col(2),num_cols)', linspace(start_col(3),end_col(3),num_cols)'];
        for j = 1:length(pes_model.hv)
            plot(pes_model.theta, pes_model.lyr_ints0{i}(:,j), 'k.-', 'color', colors_p(j,:), 'markerfacecolor', colors_p(j,:), 'linewidth', 2.5);
        end
    end
    leg_labels{i}   = string(pes_model.lyr_mat{i})+"("+string(pes_model.lyr_cls{i})+")";
    leg_cols{i}     = plot(nan, 'k.-', 'color', lyr_cols{i}, 'markerfacecolor', lyr_cols{i}, 'linewidth', 2.5);
end
% -- Formatting the axes
legend([leg_cols{:}], leg_labels, 'location', 'best');
title('Photoelectron intensities vs theta');
xlabel('Emission Angle [deg.]', 'Interpreter', 'none', 'FontWeight', 'bold');
ylabel('Relative Contribution', 'Interpreter', 'none', 'FontWeight', 'bold');
xlim([0, 90]);
ylim([0, 1.00]);
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
ax.YAxisLocation    = 'right';               % 'left' | 'right' | 'origin'
% Box Styling properties
ax.Color            = [1 1 1];
ax.Box              = 'off';                % 'on' | 'off'
ax.LineWidth        = 0.75;     
ax.Layer            = 'Top';