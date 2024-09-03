function overlay_be(elements, parity, energy_lims)
% overlay_be(elements, parity, energy_lims)
%   This is a function that plots the binding energies of all the core-levels 
%   for a particular element defined by the user. This can be used to 
%   quickly view all the available data for a particular element.
%
%   IN:
%   -   element:	    M×1 cell vector of strings of the element names; e.g. "H", "He", "Si", "In"...
%   -   parity:	        either 1 or -1; 1 plots the binding energies as positive. -1 plots the binding energis as negative.
% 	-   energy_lims:    [1×2] row vector of binding energy limits of core-levels to be plotted
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
if nargin < 2; parity = 1;  end
if nargin < 3; energy_lims = [1, 130000];  end
if isempty(elements); elements = ATOM_SYMB; end
if isempty(parity); parity = 1; end
if isempty(energy_lims); energy_lims = [1, 130000]; end
% -- Ensuring the input is a string
elements = convertCharsToStrings(elements);
%% - 1 - Filing through all the elements and extracting binding energies
be = {}; cls = {};
for i = 1:length(elements);  [be{i}, cls{i}] = calc_be(elements{i}, [], [], 0); end
%% - 2 - Ensuring consistency in parity and energy limits
if parity == -1
    energy_lims = sort(-1 .* abs(energy_lims));
    for i = 1:length(be); be{i} = -1 .* be{i}; end
elseif parity == +1
    energy_lims = sort(+1 .* abs(energy_lims));
    for i = 1:length(be); be{i} = +1 .* be{i}; end
end
%% - 2 - Plotting the binding energy lines
hold on;
ax_lims = axis;
yval    = mean(ax_lims(3:4));
for i = 1:length(be)
    cols    = lines(length(be));
    for j = 1:length(be{i})
        be_ij = be{i}(j);
        if be_ij > energy_lims(1) && be{i}(j) < energy_lims(2)
            xline(be_ij, ':', 'linewidth', 1, 'color', cols(i,:), 'HandleVisibility','off');
            str = sprintf('%s%s(%.2f)', elements{i}, cls{i}{j}, be_ij);
            text(be_ij+0.1, yval, str, 'Rotation',90, 'FontWeight','normal', 'FontSize', 8, 'color', cols(i,:));
        end
    end
end
end