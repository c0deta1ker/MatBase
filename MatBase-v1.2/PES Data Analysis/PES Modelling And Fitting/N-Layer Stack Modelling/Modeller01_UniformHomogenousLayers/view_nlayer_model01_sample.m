function fig = view_nlayer_model01_sample(lyr_mat, lyr_thick, lyr_fwhm)
% fig = view_nlayer_model01_sample(lyr_mat, lyr_thick, lyr_fwhm)
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
%   -   lyr_fwhm:       Mx1 cell-vector of the Gaussian FWHM intermixing of each layer in the stack (nm)
%
%   OUT:
%   -   fig:            output MATLAB figure.

%% Validity check on inputs
% -- Default Inputs
if nargin < 3; lyr_fwhm = num2cell(zeros(size(lyr_thick))); end
% -- Ensuring the user-defined material input is a cell-array
if ~iscell(lyr_mat);    lyr_mat = cellstr(lyr_mat); end
if ~iscell(lyr_thick);  lyr_thick = num2cell(lyr_thick); end
if ~iscell(lyr_fwhm);   lyr_fwhm = num2cell(lyr_fwhm); end
%% Defining constants
Nlyrs       = length(lyr_mat);
lyr_cols    = flipud(num2cell([0.25,0.25,0.25;lines(Nlyrs-1)], 2));
%% 1    :   Extracting the atomic concentration vs depth
xpad = 2.5;
xdat = 0:0.01:sum(cell2mat(lyr_thick(1:end-1)))+xpad;
for i = 1:Nlyrs
    if i == 1;          lyr_x0{i} = 0.5 * sum(cell2mat(lyr_thick(i)));
    elseif i == Nlyrs;  lyr_x0{i} = sum(cell2mat(lyr_thick(1:i-1))) + xpad;
    else;               lyr_x0{i} = sum(cell2mat(lyr_thick(1:i-1))) + 0.5 * sum(cell2mat(lyr_thick(i)));
    end
    if i == Nlyrs;      lyr_atmconc{i} = TopHatModel_Gauss(xdat, lyr_x0{i}, 100, 2*xpad, lyr_fwhm{i});
    else;               lyr_atmconc{i} = TopHatModel_Gauss(xdat, lyr_x0{i}, 100, lyr_thick{i}, lyr_fwhm{i});
    end
end
%% 2    :    Plotting the the material stack
% -- Initialising the figure
fig = figure(); hold on;
fig.Position(3) = 400;
fig.Position(4) = 400;
%% 2.1    :   Plotting the the material stack
nexttile([1,1]); hold on; 
x_width	= 5;
bulk_thickness = 5;
if Nlyrs == 1
    y_cum   = cell2mat(lyr_thick);
else
    y_cum   = cumsum(cell2mat(lyr_thick));
    y_cum(isinf(y_cum)) = y_cum(end-1) + bulk_thickness;
end
% -- Plotting each layer from bottom-to-up
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
        'color', 'k', 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'FontWeight','bold', 'FontSize',9);
    else
        text(0, y_loc, sprintf("(L%i) %s [%.2f nm]",i, lyr_mat{i}, lyr_thick{i}),...
        'color', 'k', 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'FontWeight','bold', 'FontSize',9); 
    end
end
axis([-x_width, x_width, -1*max(y_cum(:)), 0]);
axis square;
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
ylabel('Depth [nm]', 'fontweight', 'bold');

end