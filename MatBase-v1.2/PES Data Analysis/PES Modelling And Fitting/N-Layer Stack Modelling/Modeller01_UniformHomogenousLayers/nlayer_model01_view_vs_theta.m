function fig = nlayer_model01_view_vs_theta(pes_model)
% fig = nlayer_model01_view_vs_theta(pes_model, lyr_cols)
%   This function plots the solutions to the 'nlayer_pes_model()' function;
%   the n-layered sample stack and the corresponding photoelectron
%   contribution from each one of the layers. Specifically used for
%   plotting the emission angle dependence.
%
%   IN:
%   -   pes_model:      data structure that contains all the pes model parameters and variables (from 'nlayer_pes_model()').
%
%   OUT:
%   -   fig:            output MATLAB figure.

%% Defining constants
lyr_cols    = flipud(num2cell([0.25,0.25,0.25;lines(pes_model.Nlyrs-1)], 2));
%% Creating a figure object
fig = figure(); 
fig.Position(1) = 100; fig.Position(2) = 100;
fig.Position(3) = 900; 
fig.Position(4) = 450;
% -- Creating a tiled axis
t = tiledlayout(1,3);
t.TileSpacing = 'compact';
t.Padding = 'compact';
%% 1.1    :    Plotting the the material stack
nexttile([1]); hold on;
% -- Extracting the extent of each layer
x_width	= 5;
bulk_thickness = 5;
if pes_model.Nlyrs == 1
    y_cum   = cell2mat(pes_model.lyr_thick);
else
    y_cum   = cumsum(cell2mat(pes_model.lyr_thick));
    y_cum(isinf(y_cum)) = y_cum(end-1) + bulk_thickness;
end
% -- Plotting each layer from bottom-to-up
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
        'color', 'k', 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'FontWeight','bold', 'FontSize',9);
    else
        text(0, y_loc, sprintf("(L%i) %s [%.2f nm]",i, pes_model.lyr_mat{i}, pes_model.lyr_thick{i}),...
        'color', 'k', 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'FontWeight','bold', 'FontSize',9); 
    end
end
%% 1.2    :    Box Styling and Axis properties
ylabel('Depth From Surface [nm]', 'fontweight', 'bold');
axis([-x_width, x_width, -1*max(y_cum(:)), 0]);
% -- Formatting the axes
ax = gca;
% -- Font properties
ax.FontName         = 'Segoe UI'; 
ax.FontWeight       = 'normal'; 
ax.FontSize         = 10;
% -- Tick properties
ax.XMinorTick       = 'off'; 
ax.YMinorTick       = 'off';
ax.TickDir          = 'both';
ax.XColor           = [0 0 0]; 
ax.YColor           = [0 0 0];
% -- Ruler properties
ax.XAxisLocation    = 'bottom';             % 'bottom' | 'top' | 'origin'
ax.YAxisLocation    = 'left';               % 'left' | 'right' | 'origin'
% -- Box Styling properties
ax.Color            = [1 1 1];
ax.Box              = 'off';                % 'on' | 'off'
ax.LineWidth        = 0.75;     
ax.Layer            = 'Top';
%% 2.1    :    Plotting the intensity profiles for each layer
nexttile([1,2]); hold on; grid on; grid minor;
% - Filing through each layer
for i = 1:pes_model.Nlyrs
    % -- If only one photon energy is given
    if length(pes_model.hv) == 1
        plot(pes_model.theta, pes_model.lyr_Inorm{i}, 'k.-', 'color', lyr_cols{i}, 'markerfacecolor', lyr_cols{i}, 'linewidth', 2.5);
    % -- If a range of angles are to be plotted
    else
        num_cols    = length(pes_model.hv);
        start_col   = lyr_cols{i};
        end_col     = lyr_cols{i} + 0.50*(1-lyr_cols{i});
        colors_p = [linspace(start_col(1),end_col(1),num_cols)', linspace(start_col(2),end_col(2),num_cols)', linspace(start_col(3),end_col(3),num_cols)'];
        for j = 1:length(pes_model.hv)
            plot(pes_model.theta, pes_model.lyr_Inorm{i}(j,:), 'k.-', 'color', colors_p(j,:), 'markerfacecolor', colors_p(j,:), 'linewidth', 2.5);
        end
    end
    leg_labels{i}   = string(pes_model.lyr_mat{i})+"("+string(pes_model.lyr_cls{i})+")";
    leg_cols{i}     = plot(nan, 'k.-', 'color', lyr_cols{i}, 'markerfacecolor', lyr_cols{i}, 'linewidth', 2.5);
end
%% 2.2    :    Box Styling and Axis properties
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