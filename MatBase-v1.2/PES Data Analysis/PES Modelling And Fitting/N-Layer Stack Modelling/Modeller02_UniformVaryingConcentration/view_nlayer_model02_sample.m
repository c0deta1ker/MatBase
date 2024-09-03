function fig = view_nlayer_model02_sample(lyr_mat, lyr_thick, lyr_conc)
% fig = view_nlayer_model02_sample(lyr_mat, lyr_thick, lyr_conc)
%   This function generates a plot of a user-defined sample stack. The stack 
%   consists of N layers, each with a specified material and thickness. The 
%   layers are defined using a top-down approach, starting from the surface 
%   and ending with the bulk, where 'lyr_mat' is used to define the materials 
%   for each layer and 'lyr_thick' is used to specify the thickness of each
%   layer. In order to specify the bulk layer, the final 'lyr_thick' value
%   should be 'Inf'.
%
%   IN:
%   -   lyr_mat:        Mx1 cell-vector of the material for each layer in the stack; e.g. "Si", "SiO2", "Al2O3"...
%   -   lyr_thick:      Mx1 cell-vector of the thickness of each layer in the stack (nm)
%   -   lyr_conc:       Mx1 cell-vector of the Gaussian FWHM intermixing of each layer in the stack (nm)
%
%   OUT:
%   -   fig:            output MATLAB figure.

%% Validity check on inputs
% -- Default Inputs
if nargin < 3; lyr_conc = num2cell(100*ones(size(lyr_thick))); end
% -- Ensuring the user-defined material input is a cell-array
if ~iscell(lyr_mat);    lyr_mat = cellstr(lyr_mat); end
if ~iscell(lyr_thick);  lyr_thick = num2cell(lyr_thick); end
if ~iscell(lyr_conc);   lyr_conc = num2cell(lyr_conc); end
%% Defining constants
Nlyrs       = length(lyr_mat);
lyr_cols    = flipud(num2cell([0.25,0.25,0.25;lines(Nlyrs-1)], 2));
%% 1    :   Extracting the atomic concentration vs depth
xpad = 2.5;
for i = 1:Nlyrs
    if i == 1;          lyr_x0{i} = 0.5 * sum(cell2mat(lyr_thick(i)));
    elseif i == Nlyrs;  lyr_x0{i} = sum(cell2mat(lyr_thick(1:i-1))) + xpad;
    else;               lyr_x0{i} = sum(cell2mat(lyr_thick(1:i-1))) + 0.5 * sum(cell2mat(lyr_thick(i)));
    end
end
%% 2    :    Plotting the the material stack
% -- Initialising the figure
fig = figure(); hold on;
fig.Position(3) = 400;
fig.Position(4) = 400;
%% 2.1    :   Plotting the the material stack
nexttile([1,1]); hold on; 
bulk_thickness = 5;
if Nlyrs == 1
    y_cum   = cell2mat(lyr_thick);
else
    y_cum   = cumsum(cell2mat(lyr_thick));
    y_cum(isinf(y_cum)) = y_cum(end-1) + bulk_thickness;
end
% -- Plotting each layer from bottom-to-up
max_conc = []; labels = {};
for i = 1:Nlyrs
    for j = 1:length(lyr_conc{i})
        if j == 1
            if i == 1
                patch([0, 0, 1, 1, 0].*lyr_conc{i}{j}, [0, y_cum(i), y_cum(i), 0, 0].*-1,...
                    lyr_cols{i}, 'edgecolor', [0 0 0]);
            else
                patch([0, 0, 1, 1, 0].*lyr_conc{i}{j}, [y_cum(i-1), y_cum(i), y_cum(i), y_cum(i-1), y_cum(i-1)].*-1,...
                    lyr_cols{i}, 'edgecolor', [0 0 0]);
            end
        else
            if i == 1
                patch(lyr_conc{i}{j-1}+[0, 0, lyr_conc{i}{j}, lyr_conc{i}{j}, 0], [0, y_cum(i), y_cum(i), 0, 0].*-1,...
                    lyr_cols{i}./j, 'edgecolor', [0 0 0]);
            else
                patch(lyr_conc{i}{j-1}+[0, 0, lyr_conc{i}{j}, lyr_conc{i}{j}, 0], [y_cum(i-1), y_cum(i), y_cum(i), y_cum(i-1), y_cum(i-1)].*-1,...
                    lyr_cols{i}./j, 'edgecolor', [0 0 0]);
            end
        end
        labels{end+1} = lyr_mat{i}{j};
    end
    max_conc = vertcat(max_conc, sum(cell2mat(lyr_conc{i})));
end
% -- Adding text for each material type
for i = 1:Nlyrs
    if i == 1;  y_loc = 0 - 0.5*(y_cum(i) - 0);
    else;       y_loc = -y_cum(i-1) - 0.5*(y_cum(i) - y_cum(i-1));
    end
    for j = 1:length(lyr_conc{i})
        if j == 1
            if i == Nlyrs
                text(0.5*lyr_conc{i}{j}, y_loc, "(Bulk) "+string(lyr_mat{i}{j}),...
                'color', 'k', 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'FontWeight','bold', 'FontSize',9);
            else
                text(0.5*lyr_conc{i}{j}, y_loc, sprintf("(L%i) %s",i, lyr_mat{i}{j}),...
                'color', 'k', 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'FontWeight','bold', 'FontSize',9); 
            end
        else
            if i == Nlyrs
                text(lyr_conc{i}{j-1}+0.5*lyr_conc{i}{j}, y_loc, "(Bulk) "+string(lyr_mat{i}{j}),...
                'color', 'k', 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'FontWeight','bold', 'FontSize',9);
            else
                text(lyr_conc{i}{j-1}+0.5*lyr_conc{i}{j}, y_loc, sprintf("(L%i) %s",i, lyr_mat{i}{j}),...
                'color', 'k', 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'FontWeight','bold', 'FontSize',9); 
            end
        end
    end
end
% -- Finding the maximum concentration
max_conc = max(max_conc(:));
axis([0, max_conc, -1*max(y_cum(:)), 0]);
% legend(labels, 'location', 'eastoutside');
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

end