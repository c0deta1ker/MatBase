function [be, cls] = be_trzh2018(element, corelevel, plot_results)
% [be, cls] = be_trzh2018(element, corelevel, plot_results)
%   This is a function that extracts the electron binding energies from elements
%   with Z from 1 to 98 of the individual subshells. This is from the original 
%   work of M.B. Trzhaskovskaya and V.G. Yarzhemsky [1], which has now been 
%   digitised here for use in MATLAB.
%   [1] M.B. Trzhaskovskaya, V.G. Yarzhemsky. Atomic Data and Nuclear Data Tables 119 (2018) 99–174. Web. doi:10.1016/j.adt.2017.04.003
%
%   IN:
%   -   element:    	string of the element; e.g. "H", "He", "Si", "In"...
%   -   corelevel:      M×1 cell of strings of the core-levels to be probed; e.g. "5s1", "5p1", "5p3", "5d3", "5d5", "5f5', "5f7"...
%   -   plot_results:   if 1, will plot figure summary, otherwise it wont.
%
%   OUT:
%   -   be:             M×1 vector of the binding energies of the chosen core-levels [eV]. Returns NaN for undefined core-level energies.
%   -   cls:            M×1 vector of the core-level labels.

%% Default parameters
if nargin < 2; corelevel = [];  end
if nargin < 3; plot_results = 0;  end
if isempty(corelevel); corelevel = []; end
if isempty(plot_results); plot_results = 0; end
%% Validity checks on the input parameters
corelevel   = string(corelevel);
% -- Ensuring the inputs are in the form of 1D column vectors
if size(corelevel, 2) > 1;  corelevel = corelevel'; end
%% - 1 - Loading the MATLAB data structure
BE_DB_Trzh2018	= load('BE_DB_Trzh2018.mat'); BE_DB_Trzh2018 = BE_DB_Trzh2018.BE_DB_Trzh2018;
ATOM_SYMB   = BE_DB_Trzh2018.ATOM_SYMB;
ATOM_CL     = BE_DB_Trzh2018.ATOM_CL;
ATOM_BE     = BE_DB_Trzh2018.BE;
%% - 2 - Find the database index of the element
% -- Find the element
if ~isempty(element)
    % -- Find the index for the element of choice (case-insensitive)
    element  	= char(element);
    ele_indx 	= find(strcmpi(ATOM_SYMB, element));
    if isempty(ele_indx); msg = 'Element could not be identified. Only use atomic-symbols for elements 1 - 98; H, He, Li, Be..., Bk, Cf'; error(msg); end
% -- If no element is defined, return an error
else; msg = 'Element could not be identified. Only use atomic-symbols for elements 1 - 98; H, He, Li, Be..., Bk, Cf'; error(msg);
end
%% - 3 - Find the database index of the core-level
% -- Find the index for the element of choice (case-insensitive)
if ~isempty(corelevel)
    cl_indx = [];
    % --- If 1 core-level is entered
    if size(corelevel, 1) == 1
        % -- Try to match the user input directly with database
        corelevel   = char(corelevel);
        cl_indx 	= find(strcmpi(ATOM_CL, corelevel));
        if isempty(cl_indx)
            msg = sprintf("Core-level %s not found. Only the following exist for %s : %s . NaN values returned.", corelevel, element, join(string(ATOM_CL), ', ')); 
            warning(msg); cl_indx = 0;
        end
    % --- If >1 core-level is entered
    else
        % ---- Filing through all the core-levels defined
        for i = 1:size(corelevel, 1)
            % -- Try to match the user input directly with database
            ith_corelevel   = char(corelevel(i));
            cl_indx(i)  = find(strcmpi(ATOM_CL, ith_corelevel));
            if isempty(ith_cl_indx)
                msg = sprintf("Core-level %s not found. Only the following exist for %s : %s . NaN values returned.", corelevel(i), element, join(string(ATOM_CL), ', ')); 
                warning(msg); cl_indx(i) = 0;
            else; cl_indx(i) = ith_cl_indx; 
            end
        end
    end
% -- If no core-level is defined, use them all
else
    corelevel = string(ATOM_CL); cl_indx = 1:length(ATOM_CL); 
    if size(corelevel, 2) > 1;  corelevel = corelevel'; end
end
%% - 4 - Extracting the relevant binding energies
be = []; cls = "";
ATOM_BE_ARRAY = table2array(ATOM_BE);
ATOM_CL_ARRAY = string(ATOM_CL);
for i = 1:length(cl_indx)
    if cl_indx(i) == 0
        cls(i)  = NaN(1);
        be(i)   = NaN(1);
    else
        cls(i)  = ATOM_CL_ARRAY(cl_indx(i));
        be(i)   = ATOM_BE_ARRAY(ele_indx,cl_indx(i));
    end
end
%% Validity check o nthe outputs
if size(cls, 2) > 1;  cls = cls'; end
if size(be, 2) > 1;  be = be'; end
%% -- Plot for debugging
if plot_results == 1
    nCL = length(be);
    figX = figure(); figX.Position(3) = 750; figX.Position(4) = 300;
    hold on; grid on;
    if nCL == 1
        stem(be, 1, '-', 'linewidth', 1.5, 'marker', 'none');
        text(be, 1, sprintf('%s(%.2f)', cls(1), be), 'Rotation',45, 'FontWeight','bold', 'FontSize',8);
    else
        for i = 1:nCL
            stem(be(i), i/nCL, '-', 'linewidth', 1.5, 'marker', 'none');
            text(be(i), i/nCL, sprintf('%s(%.2f)', cls(i), be(i)), 'Rotation',45, 'FontWeight','bold', 'FontSize',8);
        end
    end
    text(0.03, 0.95, sprintf("%s (Z=%i)", element, ele_indx),...
        'FontSize', 12, 'color', 'k', 'Units','normalized', 'FontWeight', 'bold', 'HorizontalAlignment','left');
    title('Trzhaskovskaya 2018 - Binding Energy (eV)');
    xlabel('$$ \bf Binding\ Energy\ [eV] $$', 'interpreter', 'latex');
    ylabel('$$ \bf Intensity\ [arb.] $$', 'interpreter', 'latex');
    ax = gca; ax.YScale = 'linear'; ax.XScale = 'log';
    axis([1, 130000, 0, 1.40]);
end
end