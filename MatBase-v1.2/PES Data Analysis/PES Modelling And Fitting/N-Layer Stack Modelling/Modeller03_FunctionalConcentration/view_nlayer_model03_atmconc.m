function fig = view_nlayer_model03_atmconc(lyr_mat, lyr_type, lyr_conc, lyr_thick, lyr_cdl)
% fig = view_nlayer_model03_atmconc(lyr_mat, lyr_type, lyr_conc, lyr_thick, lyr_cdl)
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
% -- Ensuring the user-defined material input is a cell-array
if ~iscell(lyr_mat);    lyr_mat = cellstr(lyr_mat); end
if ~iscell(lyr_type);   lyr_type = cellstr(lyr_type); end
if ~iscell(lyr_conc);   lyr_conc = num2cell(lyr_conc); end
if ~iscell(lyr_thick);  lyr_thick = num2cell(lyr_thick); end
if ~iscell(lyr_cdl);    lyr_cdl = num2cell(lyr_cdl); end
%% Extracting constants
Nlyrs       = length(lyr_mat);
lyr_cols    = flipud(num2cell([0.25,0.25,0.25;lines(Nlyrs-1)], 2));
%% 1    :   Mid-point depth of all the layers
for i = 1:Nlyrs
    if i == 1;          lyr_z0{i} = 0.5 * sum(cell2mat(lyr_thick(i)));
    elseif i == Nlyrs;  lyr_z0{i} = sum(cell2mat(lyr_thick(1:i-1)));
    else;               lyr_z0{i} = sum(cell2mat(lyr_thick(1:i-1))) + 0.5 * sum(cell2mat(lyr_thick(i)));
    end
end
%% 2    :   Atomic Concentrations vs depth
zpad = 6.5;
lyr_z = 0:0.001:sum(cell2mat(lyr_thick(1:end-1)))+zpad; lyr_z = lyr_z';
for i = 1:Nlyrs
    % -- For the bulk layer
    if i == Nlyrs
        if strcmpi(lyr_type{i}, "StLHS");               lyr_atmconc{i} = StepModel_LHS(lyr_z, lyr_z0{i}, lyr_conc{i});
        elseif strcmpi(lyr_type{i}, "StExLHS");         lyr_atmconc{i} = StepModel_Exp_LHS(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "StGaLHS");         lyr_atmconc{i} = StepModel_Gauss_LHS(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "ToHaExTrLHS");     lyr_atmconc{i} = StepModel_Exp_LHS_trunc(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i}, lyr_thick{i-1});
        elseif strcmpi(lyr_type{i}, "ToHaGaTrLHS");     lyr_atmconc{i} = StepModel_Gauss_LHS_trunc(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i}, lyr_thick{i-1});
        end
    % -- For each layer
    else
        if strcmpi(lyr_type{i}, "StLHS");           lyr_atmconc{i} = StepModel_LHS(lyr_z, lyr_z0{i}, lyr_conc{i});
        elseif strcmpi(lyr_type{i}, "StRHS");       lyr_atmconc{i} = StepModel_RHS(lyr_z, lyr_z0{i}, lyr_conc{i});
        elseif strcmpi(lyr_type{i}, "StExLHS");     lyr_atmconc{i} = StepModel_Exp_LHS(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "StExRHS");     lyr_atmconc{i} = StepModel_Exp_RHS(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "StGaLHS");     lyr_atmconc{i} = StepModel_Gauss_LHS(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "StGaRHS");     lyr_atmconc{i} = StepModel_Gauss_RHS(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "StExTrLHS");   lyr_atmconc{i} = StepModel_Exp_LHS_trunc(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i}, lyr_thick{i-1});
        elseif strcmpi(lyr_type{i}, "StExTrRHS");   lyr_atmconc{i} = StepModel_Exp_RHS_trunc(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i}, lyr_thick{i+1});
        elseif strcmpi(lyr_type{i}, "StGaTrLHS");   lyr_atmconc{i} = StepModel_Gauss_LHS_trunc(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i}, lyr_thick{i-1});
        elseif strcmpi(lyr_type{i}, "StGaTrRHS");   lyr_atmconc{i} = StepModel_Gauss_RHS_trunc(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i}, lyr_thick{i+1});
        elseif strcmpi(lyr_type{i}, "ToHa");        lyr_atmconc{i} = TopHatModel(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i});
        elseif strcmpi(lyr_type{i}, "ToHaEx");      lyr_atmconc{i} = TopHatModel_Exp(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "ToHaExLHS");   lyr_atmconc{i} = TopHatModel_Exp_LHS(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "ToHaExRHS");   lyr_atmconc{i} = TopHatModel_Exp_RHS(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "ToHaGa");      lyr_atmconc{i} = TopHatModel_Gauss(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "ToHaGaLHS");   lyr_atmconc{i} = TopHatModel_Gauss_LHS(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "ToHaGaRHS");   lyr_atmconc{i} = TopHatModel_Gauss_LHS(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "ToHaExTrLHS"); lyr_atmconc{i} = TopHatModel_Exp_LHS_trunc(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i}, lyr_cdl{i}, lyr_thick{i-1});
        elseif strcmpi(lyr_type{i}, "ToHaExTrRHS"); lyr_atmconc{i} = TopHatModel_Exp_RHS_trunc(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i}, lyr_cdl{i}, lyr_thick{i+1});
        elseif strcmpi(lyr_type{i}, "ToHaGaTrLHS"); lyr_atmconc{i} = TopHatModel_Gauss_LHS_trunc(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i}, lyr_cdl{i}, lyr_thick{i-1});
        elseif strcmpi(lyr_type{i}, "ToHaGaTrRHS"); lyr_atmconc{i} = TopHatModel_Gauss_RHS_trunc(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i}, lyr_cdl{i}, lyr_thick{i+1});
        end
    end
end
%% 3    :   Initialising the figure
% -- Initialising the figure
fig = figure(); hold on;
fig.Position(3) = 900;
fig.Position(4) = 400;
% -- Creating a tiled axis
t = tiledlayout(1,3);
t.TileSpacing = 'compact';
t.Padding = 'compact';
%% 3.1    :   Plotting the the material stack
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
title(sprintf("Sample Stack"),'FontWeight','bold','interpreter', 'none', 'fontsize', 12);

%% 3.2    :   Plotting the atomic concentration with depth
nexttile([1,2]); hold on; 
for i = 1:Nlyrs
    area(lyr_z, lyr_atmconc{i}, 0, 'facecolor', lyr_cols{i}, 'edgecolor', 'none', 'facealpha', 0.5);
    plot(lyr_z, lyr_atmconc{i}, 'k-', 'color', lyr_cols{i}, 'linewidth', 2);
end
% -- Adding text for each material type
for i = 1:Nlyrs
    if i == Nlyrs
        text(lyr_z0{i}+0.5*zpad, 120, sprintf("(Bulk) \n %s", lyr_mat{i}),...
        'color', 'k', 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'FontWeight','bold', 'FontSize',9);
    else
        text(lyr_z0{i}, 120, sprintf("(L%i)\n %s \n %.1f nm",i, lyr_mat{i}, lyr_thick{i}),...
        'color', 'k', 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'FontWeight','bold', 'FontSize',9); 
    end
end
% -- Formatting the axes
xlabel('Depth (nm)', 'fontweight', 'bold');
ylabel('Concentration (at. %)', 'fontweight', 'bold');
axis([0, max(lyr_z(:)), 0, 140]);
% -- Formatting the axes
ax = gca;
% -- Font properties
ax.FontName         = 'Segoe UI'; 
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
ax.Box              = 'on';                % 'on' | 'off'
title(sprintf("Atomic Concentration versus Depth"),'FontWeight','bold','interpreter', 'none', 'fontsize', 12);

end