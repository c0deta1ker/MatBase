function fig = view_be(elements, energy_lims)
% fig = view_be(elements)
%   This is a function that plots the binding energies of all the core-levels 
%   for a particular element defined by the user. This can be used to 
%   quickly view all the available data for a particular element.
%
%   IN:
%   -   element:	    M×1 cell vector of strings of the element names; e.g. "H", "He", "Si", "In"...
% 	-   energy_lims:    [1×2] row vector of energy-axis limits
%
%   OUT:
%   -   fig:	    figure output

%% - Loading the MATLAB data structure
BE_DB_Constantinou2023	= load('BE_DB_Constantinou2023.mat'); BE_DB_Constantinou2023 = BE_DB_Constantinou2023.BE_DB_Constantinou2023;
ATOM_SYMB   = BE_DB_Constantinou2023.ATOM_SYMB;
ATOM_CL     = BE_DB_Constantinou2023.ATOM_CL;
ATOM_BE     = BE_DB_Constantinou2023.BE;
%% Default parameters
if nargin < 1; elements = ATOM_SYMB;  end
if nargin < 2; energy_lims = [0.1, 130000];  end
if isempty(elements); elements = ATOM_SYMB; end
if isempty(energy_lims); energy_lims = [0.1, 130000]; end
% -- Ensuring the input is a string
elements = convertCharsToStrings(elements);
%% - 1 - Filing through all the elements and extracting binding energies
be = {}; cls = {};
for i = 1:length(elements);  [be{i}, cls{i}] = calc_be(elements{i}, [], [], 0); end
%% - 2 - Plotting the binding energy spectrum
fig = figure(); hold on;
fig.Position(3) = 800; 
fig.Position(4) = 400;
for i = 1:length(be)
    cols    = lines(length(be));
    for j = 1:length(be{i})
        if j == 1; stem(be{i}(j), i/length(be), '-', 'linewidth', 1.0, 'marker', 'none', 'color', cols(i,:));
        else;  stem(be{i}(j), i/length(be), '-', 'linewidth', 1.0, 'marker', 'none', 'color', cols(i,:), 'HandleVisibility','off');
        end
        text(be{i}(j), i/length(be), sprintf('%s%s(%.2f)', elements{i}, cls{i}{j}, be{i}(j)), 'Rotation',90, 'FontWeight','normal', 'FontSize',8);
    end
end
legend(elements, 'location', 'southeastoutside', 'FontSize', 6);
title('Binding Energy Spectrum');
xlabel('$$ \bf Binding\ Energy\ [eV] $$', 'interpreter', 'latex');
ylabel('$$ \bf arb. $$', 'interpreter', 'latex');
ax = gca; ax.YScale = 'linear'; ax.XScale = 'log';
axis([energy_lims(1), energy_lims(2), 0, 1.40]);
end