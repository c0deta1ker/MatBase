function xsect_delta = xsect_delta_trzh2018(element, corelevel, hv, plot_results)
% xsect_delta = xsect_delta_trzh2018(element, corelevel, hv, plot_results)
%   This is a function that calculates the cross sections for photoionization
%   for atomic subshells with binding energies lower than 1.5 keV of all elements with 1 ≤ Z ≤
%   100 in the photon energy range 1.5–10 keV. The calculations were performed in an effort to provide handy
%   theoretical data for experimental studies by hard X-ray photoelectron spectroscopy (HAXPES).
%   The relativistic treatment of atomic photoeffect and the Dirac–Fock method with proper consideration of
%   the electron exchange for computing the electron wave functions. The photoionization cross sections were
%   determined including all multipoles of the radiative field while the photoelectron angular distribution
%   parameters were obtained within the quadrupole approximation. The effect of the hole resulting in the
%   atomic subshell after photoionization was taken into account by the use of the frozen orbital model.
%   This is from the original work of M.B. Trzhaskovskaya and V.G. Yarzhemsky [1], 
%   which has now been digitised here for use in MATLAB.
%   [1] M.B. Trzhaskovskaya, V.G. Yarzhemsky. Atomic Data and Nuclear Data Tables 119 (2018) 99–174. Web. doi:10.1016/j.adt.2017.04.003
%
%   IN:
%   -   element:    	string of the element; e.g. "H", "He", "Si", "In"...
%   -   corelevel:      string or M×1 string vector of the core-levels to be probed; e.g. "1s1", "2p1", "2p3", "3d3", "3d5", "5f5', "5f7"...
%   -   hv:             scalar or N×1 vector of the incident photon energies [eV]
%   -   plot_results:   if 1, will plot figure summary, otherwise it wont.
%
%   OUT:
%   -   xsect_delta:     N×M vector of the photoelectron angular distribution asymmetry (non-dipole) parameter delta.

%% Default parameters
if nargin < 2; corelevel = [];  end
if nargin < 3; hv = 1e3:10:1e4;  end
if nargin < 4; plot_results = 0;  end
if isempty(corelevel); corelevel = []; end
if isempty(hv); hv = 1e3:10:1e4; end
if isempty(plot_results); plot_results = 0; end
%% Validity checks on the input parameters
% -- Forcing the photon energy vector to have no duplicates and in ascending order
hv          = sort(unique(hv));
corelevel   = string(corelevel);
% -- Ensuring the inputs are in the form of 1D column vectors
if size(hv, 2) > 1;         hv = hv'; end
if size(corelevel, 2) > 1;  corelevel = corelevel'; end
%% - 1 - Loading the MATLAB data structure that contains the 1973 Schofield Cross-Sections
XS_DB_Trzh2018	= load('XS_DB_Trzh2018.mat'); XS_DB_Trzh2018 = XS_DB_Trzh2018.XS_DB_Trzh2018;
ATOM_SYMB   = XS_DB_Trzh2018.ATOM_SYMB;
HV          = XS_DB_Trzh2018.HV;
XSECT       = XS_DB_Trzh2018.XSECT_DELTA;
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
% -- Extracting the relevant table
table_hv    = HV;
table_xs    = XSECT{ele_indx};
%% - 3 - Find the database index of the core-level
% -- Finding the list of possible core-levels
ATOM_CL = table_xs.Properties.VariableNames;
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
            % warning(msg); 
            cl_indx = 0;
        end
    % --- If >1 core-level is entered
    else
        % ---- Filing through all the core-levels defined
        for i = 1:size(corelevel, 1)
            % -- Try to match the user input directly with database
            ith_corelevel   = char(corelevel(i));
            ith_cl_indx     = find(strcmpi(ATOM_CL, ith_corelevel));
            if isempty(ith_cl_indx)
                msg = sprintf("Core-level %s not found. Only the following exist for %s : %s . NaN values returned.", corelevel(i), element, join(string(ATOM_CL), ', ')); 
                % warning(msg); 
                cl_indx(i) = 0;
            else; cl_indx(i) = ith_cl_indx; 
            end
        end
    end
% -- If no core-level is defined, use them all
else
    corelevel = string(ATOM_CL); cl_indx = 1:length(ATOM_CL); 
    if size(corelevel, 2) > 1;  corelevel = corelevel'; end
end
%% - 4 - Extracting the relevant photoionisation parameter
hv_data     = table_hv.(1);
xsect_data  = [];
for i = 1:length(cl_indx)
    if cl_indx(i) == 0;     xsect_data(:,i)   = NaN(size(hv_data));
    else;                   xsect_data(:,i)   = table_xs.(cl_indx(i));
    end
end
%% - 5 - Interpolating the cross section values from the database
% - Ensuring the data is in the form of 1D column vectors
if size(hv_data, 2) > 1;    hv_data = hv_data'; end
if size(hv, 2) > 1;         hv = hv'; end
xsect_delta = interp1(hv_data, xsect_data, hv, 'linear');
%% -- Plot for debugging
if plot_results == 1
    nCL     = size(xsect_delta, 2);
    cols    = lines(nCL);
    figure(); hold on; grid on;
    for i = 1:nCL
        plot(hv, xsect_delta(:,i), 'o-', 'markersize', 2, 'markeredgecolor', cols(i,:), 'markerfacecolor', cols(i,:), 'color', cols(i,:)); 
        scatter(hv_data, xsect_data(:,i), 'SizeData', 40, 'markeredgecolor', cols(i,:), 'markerfacecolor', cols(i,:), 'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.6, 'HandleVisibility','off');
    end
    legend(corelevel, 'location', 'eastoutside', 'FontSize', 6);
    text(0.45, 0.07, sprintf("%s (Z=%i)", element, ele_indx),...
        'FontSize', 12, 'color', 'k', 'Units','normalized', 'FontWeight', 'bold', 'HorizontalAlignment','center');
    title('Trzh. & Yarz. 2018 - Photoionisation Parameter (delta)');
    xlabel('$$ \bf Photon\ Energy\ [eV] $$', 'interpreter', 'latex');
    ylabel('$$ \bf \delta $$', 'interpreter', 'latex');
    ax = gca; ax.YScale = 'linear'; ax.XScale = 'linear';
    yline(0, 'k-', 'LineWidth',1, 'HandleVisibility','off');
    axis([900, 11000, 0, 0.55]);
end
end