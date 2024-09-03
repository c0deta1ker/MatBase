function fig = nlayer_model02_view_vs_hv(pes_model)
% fig = nlayer_model01_view_vs_hv(pes_model)
%   This function plots the solutions to the 'nlayer_pes_model()' function;
%   the n-layered sample stack and the corresponding photoelectron
%   contribution from each one of the layers. Specifically used for
%   plotting the photon energy dependence.
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
bulk_thickness = 5;
if pes_model.Nlyrs == 1
    y_cum   = cell2mat(pes_model.lyr_thick);
else
    y_cum   = cumsum(cell2mat(pes_model.lyr_thick));
    y_cum(isinf(y_cum)) = y_cum(end-1) + bulk_thickness;
end
% -- Plotting each layer from bottom-to-up
max_conc = []; labels = {};
for i = 1:pes_model.Nlyrs
    for j = 1:length(pes_model.lyr_conc{i})
        if j == 1
            if i == 1
                patch([0, 0, 1, 1, 0].*pes_model.lyr_conc{i}{j}, [0, y_cum(i), y_cum(i), 0, 0].*-1,...
                    lyr_cols{i}, 'edgecolor', [0 0 0]);
            else
                patch([0, 0, 1, 1, 0].*pes_model.lyr_conc{i}{j}, [y_cum(i-1), y_cum(i), y_cum(i), y_cum(i-1), y_cum(i-1)].*-1,...
                    lyr_cols{i}, 'edgecolor', [0 0 0]);
            end
        else
            if i == 1
                patch(pes_model.lyr_conc{i}{j-1}+[0, 0, pes_model.lyr_conc{i}{j}, pes_model.lyr_conc{i}{j}, 0], [0, y_cum(i), y_cum(i), 0, 0].*-1,...
                    lyr_cols{i}./j, 'edgecolor', [0 0 0]);
            else
                patch(pes_model.lyr_conc{i}{j-1}+[0, 0, pes_model.lyr_conc{i}{j}, pes_model.lyr_conc{i}{j}, 0], [y_cum(i-1), y_cum(i), y_cum(i), y_cum(i-1), y_cum(i-1)].*-1,...
                    lyr_cols{i}./j, 'edgecolor', [0 0 0]);
            end
        end
        labels{end+1} = pes_model.lyr_mat{i}{j};
    end
    max_conc = vertcat(max_conc, sum(cell2mat(pes_model.lyr_conc{i})));
end
% -- Finding the maximum concentration
max_conc = max(max_conc(:));
axis([0, max_conc, -1*max(y_cum(:)), 0]);
legend(labels, 'location', 'eastoutside');
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
xlabel('Concentration (at. %)', 'fontweight', 'bold');
ylabel('Depth [nm]', 'fontweight', 'bold');

%% 2.1    :    Plotting the intensity profiles for each layer
nexttile([1,2]); hold on; grid on; grid minor;
% - Filing through each layer
leg_labels = {}; leg_cols = {};
for i = 1:pes_model.Nlyrs
    for j = 1:length(pes_model.lyr_conc{i})
        plot(pes_model.hv, pes_model.lyr_Inorm{i}{j}, 'k.-', 'color', lyr_cols{i}./j, 'markerfacecolor', lyr_cols{i}./j, 'linewidth', 2.5);
        leg_labels{end+1}   = string(pes_model.lyr_mat{i}{j})+"("+string(pes_model.lyr_cls{i}{j})+")";
        leg_cols{end+1}     = plot(nan, 'k.-', 'color', lyr_cols{i}./j, 'markerfacecolor', lyr_cols{i}./j, 'linewidth', 2.5);
    end
end
%% 2.2    :    Box Styling and Axis properties
legend([leg_cols{:}], leg_labels, 'location', 'best');
title('Photoelectron intensities vs photon energy');
xlabel('Photon energy [eV]', 'Interpreter', 'none', 'FontWeight', 'bold');
ylabel('Relative Contribution', 'Interpreter', 'none', 'FontWeight', 'bold');
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