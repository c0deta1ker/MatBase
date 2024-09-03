function xsect_sigma = xsect_sigma_cant2022(element, corelevel, hv, plot_results)
% xsect_sigma = xsect_sigma_cant2022(element, corelevel, hv, plot_results)
%   This is a function that calculates the cross sections for photoionization 
%   from elements with Z from 3 to 98 and for photon energies from 1.5 to 10 keV. 
%   This is an empirical functions used to describe discrete theoretically 
%   calculated values for photoemission cross sections and asymmetry parameters
%   from the individual subshells are determined. This is from the original
%   work of David J. H. Cant [1], which has now been digitised here for use in
%   MATLAB.
%   [1] David J. H. Cant, Ben F. Spencer, Wendy R. Flavell, Alexander G. Shard. Surf Interface Anal. 2022; 54(4): 442-454. doi:10.1002/sia.7059
%
%   IN:
%   -   element:    	string of the element; e.g. "H", "He", "Si", "In"...
%   -   corelevel:      string or M×1 string vector of the core-levels to be probed; e.g. "5s1/2", "5p1/2", "5p3/2", "5d3/2", "5d5/2", "5f5/2', "5f7/2"...
%   -   hv:             scalar or N×1 vector of the incident photon energies [eV]
%   -   plot_results:   if 1, will plot figure summary, otherwise it wont.
%
%   OUT:
%   -   xsect_sigma:    N×M vector of the photoionisation cross-sections [barn/atom]

%% Default parameters
if nargin < 2; corelevel = [];  end
if nargin < 3; hv = 1e2:10:1.5e4;  end
if nargin < 4; plot_results = 0;  end
if isempty(corelevel); corelevel = []; end
if isempty(hv); hv = 1e2:10:1.5e4; end
if isempty(plot_results); plot_results = 0; end
%% Validity checks on the input parameters
% -- Forcing the photon energy vector to have no duplicates and in ascending order
hv          = sort(unique(hv));
corelevel   = string(corelevel);
% -- Ensuring the inputs are in the form of 1D column vectors
if size(hv, 2) > 1;         hv = hv'; end
if size(corelevel, 2) > 1;  corelevel = corelevel'; end
%% - 1 - Defining the table of coefficients
T1_XSect_Sigma_Coeff = table();
T1_XSect_Sigma_Coeff.Shell = {'1s1';'2s1';'3s1';'4s1';'5s1'; '2p1';'2p3';'3p1';'3p3';'4p1';'4p3';'5p1';'5p3'; '3d3';'3d5';'4d3';'4d5';'4f5';'4f7'};
T1_XSect_Sigma_Coeff.a_sig = [-1.63E-10;-2.66E-11;2.76E-12;-1.02E-13;9.46E-14;-1.27E-08;-2.87E-08;9.83E-10;5.22E-09;-1.02E-10;-1.05E-10;7.97E-14;-6.01E-13;-2.10E-06;-3.90E-06;1.29E-06;1.78E-06;-3.48E-05;-4.39E-05];
T1_XSect_Sigma_Coeff.m_sig = [1.07E-10;4.30E-12;-1.61E-14;4.61E-15;-9.44E-16;1.70E-09;3.76E-09;7.18E-11;6.28E-11;3.69E-12;2.89E-12;3.79E-15;2.91E-14;8.48E-08;1.56E-07;-8.41E-09;-1.08E-08;1.20E-06;1.52E-06];
T1_XSect_Sigma_Coeff.b_sig = [2.42E+02;9.08E+01;6.42E+02;1.26E+03;8.87E+02;2.74E+02;1.91E+02;8.73E+02;3.51E+02;3.08E+03;-1.64E+02;1.29E+03;2.30E+03;9.90E+01;5.37E+01;2.70E+02;3.36E+02;1.01E+02;1.01E+02];
T1_XSect_Sigma_Coeff.n_sig = [-1.14E+01;1.07E+01;-1.48E+01;-3.90E+01;-5.51E+00;-3.46E+00;4.87E+00;-3.31E+01;5.98E-01;-1.03E+02;1.20E+01;-3.33E+01;-3.96E+01;5.18E+00;8.05E+00;5.20E-02;5.21E-02;3.79E+00;3.74E+00];
T1_XSect_Sigma_Coeff.p_sig = [1.17E+00;8.79E-01;6.05E-01;6.61E-01;-4.97E-02;6.46E-01;4.61E-01;9.51E-01;5.94E-01;1.10E+00;3.08E-01;6.30E-01;5.13E-01;3.36E-01;2.70E-01;4.84E-01;4.52E-01;2.23E-01;2.16E-01];
T1_XSect_Sigma_Coeff.q_sig = [-2.37E+00;-2.28E+00;-1.82E+00;-1.97E+00;-1.49E+00;-2.42E+00;-2.64E+00;-2.70E+00;-2.52E+00;2.10E-01;7.94E-01;-1.97E+00;-2.01E+00;-2.96E+00;-3.05E+00;-1.84E+00;-1.68E+00;-4.09E+00;-4.09E+00];
T1_XSect_Sigma_Coeff.r_sig = [-2.63E+00;-2.14E+00;-1.76E+00;-1.43E+00;-7.56E-01;-5.82E+00;-5.47E+00;-2.00E+00;-6.45E+00;-9.01E+00;-7.06E+00;-2.12E+00;-1.67E+00;-7.90E+00;-8.49E+00;-1.66E+01;-1.70E+01;-2.17E+04;-2.16E+04];
T1_XSect_Sigma_Coeff.s_sig = [4.88E+00;7.42E+00;3.14E+01;1.64E+01;9.62E+01;3.80E+00;3.37E+00;2.92E+01;4.42E+00;3.05E+01;2.10E+03;4.21E+01;5.45E+01;4.55E+00;3.44E+00;2.12E+00;1.69E+00;2.88E-02;2.87E-02];
T1_XSect_Sigma_Coeff.t_sig = [3.62E-01;3.78E-01;6.52E-01;6.14E-01;4.52E+00;2.61E-01;2.85E-01;7.78E-01;3.48E-01;1.18E-01;1.02E-01;1.01E+00;9.84E-01;3.00E-01;2.93E-01;2.20E-01;1.99E-01;3.03E-01;3.03E-01];
%% - 2 - Defining all of the atomic variables
% -- Defining the elements whose photoionisation cross-sections we have
ATOM_SYMB = {...
    'H','He',...
    'Li','Be','B','C','N','O','F','Ne',...
    'Na','Mg','Al','Si','P','S','Cl','Ar',...
    'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',...
    'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',...
    'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',...
    'Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf'};
% -- Defining the Z number for each one of these elements
ATOM_ZNUM = 1:length(ATOM_SYMB);
% -- Defining all of the core-levels that are probed
% --- Defining in terms of orbital notation
ATOM_CL = {...
    '1s1', '2s1', '3s1', '4s1', '5s1',...
    '2p1', '2p3','3p1', '3p3','4p1', '4p3','5p1', '5p3'...
    '3d3', '3d5','4d3', '4d5',...
    '4f5', '4f7',...
    };
% -- Extracting material properties
mpd = get_mpd_props(element); 
Z   = mpd.ATOM_ZNUM;
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
            ith_cl_indx     = find(strcmpi(ATOM_CL, ith_corelevel));
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
%% - 4 - Calculating the photoionization cross-sections from the equation
cm22barn = 1e+24;
% -- Filing through all core-levels
xsect_sigma = [];
for i = 1:length(cl_indx)
    if cl_indx(i) == 0; xsect_sigma(:,i)   = NaN(size(hv));
    else             
        % --- Extract the coefficients
        a = T1_XSect_Sigma_Coeff.a_sig(cl_indx(i));
        m = T1_XSect_Sigma_Coeff.m_sig(cl_indx(i));
        b = T1_XSect_Sigma_Coeff.b_sig(cl_indx(i));
        n = T1_XSect_Sigma_Coeff.n_sig(cl_indx(i));
        p = T1_XSect_Sigma_Coeff.p_sig(cl_indx(i));
        q = T1_XSect_Sigma_Coeff.q_sig(cl_indx(i));
        r = T1_XSect_Sigma_Coeff.r_sig(cl_indx(i));
        s = T1_XSect_Sigma_Coeff.s_sig(cl_indx(i));
        t = T1_XSect_Sigma_Coeff.t_sig(cl_indx(i));
        % --- Calculate sigma in barns
        xsect_sigma(:,i) = (a + m*Z).*(hv + b + n*Z + (p*Z).^2).^(q + r.*exp(-(Z/s).^t)) .* cm22barn;
    end
end
%% -- Plot for debugging
if plot_results == 1
    nCL     = length(cl_indx);
    cols    = lines(nCL);
    figure(); hold on; grid on;
    for i = 1:nCL
        plot(hv, xsect_sigma(:,i), 'o-', 'markersize', 2, 'markeredgecolor', cols(i,:), 'markerfacecolor', cols(i,:), 'color', cols(i,:)); 
    end
    legend(corelevel, 'location', 'eastoutside', 'FontSize', 6);
    text(0.45, 0.07, sprintf("%s (Z=%i)", element,Z),...
        'FontSize', 12, 'color', 'k', 'Units','normalized', 'FontWeight', 'bold', 'HorizontalAlignment','center');
    title('Cant 2022 - Photoionisation Cross Section (sigma)');
    xlabel('$$ \bf Photon\ Energy\ [eV] $$', 'interpreter', 'latex');
    ylabel('$$ \bf Cross\ Section\ [barn/atom] $$', 'interpreter', 'latex');
    ax = gca; ax.YScale = 'log'; ax.XScale = 'log';
    axis([100, 15000, 1e-8, 1e8]);
end
end