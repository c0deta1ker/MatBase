function xsect_gamma = xsect_gamma_cant2022(element, corelevel, hv, plot_results)
% xsect_gamma = xsect_gamma_cant2022(element, corelevel, hv, plot_results)
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
%   -   xsect_delta:     N×M vector of the photoelectron angular distribution asymmetry (non-dipole) parameter delta.

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
%% - 1 - Extracting the table of coefficients
% -- Table of S core-level coefficients
T3_XSect_GammaS_Coeff = table();
T3_XSect_GammaS_Coeff.Shell = {'1s1';'2s1';'3s1';'4s1';'5s1'};
T3_XSect_GammaS_Coeff.a_gamma = [1.06E+00;2.31E+00;2.58E+00;1.22E+00;1.12E+00];
T3_XSect_GammaS_Coeff.b_gamma = [4.68E-01;4.23E-01;4.46E-01;3.39E-01;3.20E-01];
T3_XSect_GammaS_Coeff.m_gamma = [5.09E-01;2.21E-01;2.36E-01;1.18E-02;1.07E-02];
T3_XSect_GammaS_Coeff.n_gamma = [2.41E+00;3.40E+00;1.78E+00;2.84E+00;3.74E+00];
T3_XSect_GammaS_Coeff.q_gamma = [1.83E+00;3.06E+00;3.04E+00;-4.69E-01;-1.96E+00];
T3_XSect_GammaS_Coeff.r_gamma = [-4.71E-01;3.90E-03;3.90E-03;1.92E-03;9.29E-03];
T3_XSect_GammaS_Coeff.p_gamma = [2.54E+00;2.18E+00;1.16E+00;5.53E-01;4.33E-01];
T3_XSect_GammaS_Coeff.s_gamma = [1.62E-01;5.98E-01;3.27E-01;4.47E-03;-2.38E-04];
% -- Table of all other core-level coefficients
T3_XSect_Gamma_Coeff = table();
T3_XSect_Gamma_Coeff.Shell = {'2p1';'2p3';'3p1';'3p3';'4p1';'4p3';'5p1';'5p3'; '3d3';'3d5';'4d3';'4d5';'5d3';'5d5';'4f5';'4f7'};
T3_XSect_Gamma_Coeff.da_gamma = [1.42E-03;1.42E-03;2.79E-03;5.94E-03;1.71E-03;1.81E-03;1.94E-01;2.76E-01;4.45E-01;3.63E-01;8.61E+05;8.61E+05;7.98E-04;8.98E-04;8.38E+05;8.62E+05];
T3_XSect_Gamma_Coeff.ga_gamma = [-1.28E+00;-1.28E+00;5.93E-01;5.09E-01;6.35E-01;6.34E-01;1.66E-01;1.61E-01;-5.08E-01;-8.37E-01;-2.10E+00;-2.13E+00;8.64E-01;7.79E-01;-2.18E+00;-2.10E+00];
T3_XSect_Gamma_Coeff.dm_gamma = [6.23E-02;5.35E-02;1.24E-02;7.02E-03;3.64E-03;3.19E-03;3.52E-01;1.42E-01;1.35E-02;1.63E-02;2.40E-03;2.56E-03;1.23E-03;1.38E-03;3.19E-03;3.91E-03];
T3_XSect_Gamma_Coeff.gm_gamma = [3.60E-01;3.82E-01;4.53E-01;5.21E-01;5.83E-01;6.10E-01;9.49E-02;1.91E-01;5.14E-01;4.91E-01;6.67E-01;6.59E-01;8.14E-01;7.36E-01;6.31E-01;6.09E-01];
T3_XSect_Gamma_Coeff.du_gamma = [3.41E+01;3.12E+01;1.76E+02;2.07E+02;1.62E+02;1.22E+02;2.26E+01;2.81E+01;9.48E+01;9.61E+01;1.18E+02;1.12E+02;6.32E+01;8.99E+01;8.40E+01;8.54E+01];
T3_XSect_Gamma_Coeff.gu_gamma = [-2.34E-01;-2.27E-01;-3.91E-01;-4.13E-01;-3.96E-01;-3.82E-01;-1.92E-01;-2.24E-01;-3.70E-01;-3.83E-01;-3.90E-01;-3.97E-01;-4.01E-01;-4.33E-01;-3.91E-01;-3.95E-01];
T3_XSect_Gamma_Coeff.cw_gamma = [4.29E-01;4.41E-01;-2.29E-01;-1.33E-01;-6.99E-02;4.91E-01;1.36E+00;1.47E+00;6.04E+00;6.22E+00;5.81E+00;6.11E+00;1.91E+00;1.42E+00;5.89E+00;5.88E+00];
T3_XSect_Gamma_Coeff.dw_gamma = [-2.82E-05;-3.20E-05;1.11E-06;-1.46E-05;9.50E-06;-1.60E-05;-8.81E-05;-9.26E-05;-7.39E-06;-6.97E-06;-2.18E-06;-1.26E-05;2.52E-02;2.51E-02;-1.17E-05;-9.40E-06];
T3_XSect_Gamma_Coeff.n_gamma = [-8.11E-01;-8.89E-01;2.40E+00;9.20E-01;5.45E+00;4.38E+00;9.22E-01;6.96E-01;8.41E-01;5.72E-01;2.97E+00;2.80E+00;2.20E-01;2.32E-01;2.55E+00;3.10E+00];
T3_XSect_Gamma_Coeff.cr_gamma = [8.98E+02;6.79E+02;6.75E+01;8.18E+01;-1.12E+02;1.52E+01;1.51E+03;1.33E+03;2.48E+02;2.49E+02;-5.22E+02;-4.92E+02;1.23E+03;1.22E+03;-4.94E+02;-4.66E+02];
T3_XSect_Gamma_Coeff.dr_gamma = [3.15E-01;3.90E-01;1.45E-01;1.70E-01;2.71E-01;2.71E-01;-1.51E-01;-1.33E-01;4.21E-02;7.04E-03;6.04E-01;5.85E-01;-1.19E-01;-1.22E-01;7.87E-01;7.21E-01];
T3_XSect_Gamma_Coeff.p_gamma = [-2.06E-01;-2.08E-01;-2.90E+00;-1.63E+00;-7.62E+00;-5.61E+00;-2.96E+01;-3.85E+01;-2.50E-01;-8.86E-03;-7.01E+00;-6.82E+00;-1.07E-01;-1.05E-01;-7.43E+00;-5.67E+00];
T3_XSect_Gamma_Coeff.cs_gamma = [-3.60E+03;-3.55E+03;1.18E+02;7.67E+01;-3.01E+01;7.64E+01;1.30E+02;2.82E+01;3.45E+02;1.01E+02;-1.34E+02;-1.43E+02;-3.36E+01;-3.36E+01;-1.22E+02;-1.22E+02];
T3_XSect_Gamma_Coeff.ds_gamma = [2.40E+00;2.37E+00;3.73E-02;3.73E-02;1.74E-01;1.91E-01;1.09E-01;2.08E-01;-3.00E-02;-1.01E-02;2.42E-01;2.31E-01;4.08E-01;3.98E-01;3.14E-01;3.60E-01];
%% - 2 - Defining all of the variables
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
    '3d3', '3d5','4d3', '4d5','5d3', '5d5',...
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
%% - 4 - Calculating the non-dipole parameter, gamma
% -- Filing through all core-levels
xsect_gamma = [];
for i = 1:length(cl_indx)
    if cl_indx(i) == 0; xsect_gamma(:,i)   = NaN(size(hv));
    else
        gamma = [];
        % --- For all S core levels, user the BetaS table of coefficients
        if cl_indx(i) < 6
            % --- Extract the coefficients
            a = T3_XSect_GammaS_Coeff.a_gamma(cl_indx(i));
            b = T3_XSect_GammaS_Coeff.b_gamma(cl_indx(i));
            m = T3_XSect_GammaS_Coeff.m_gamma(cl_indx(i));
            n = T3_XSect_GammaS_Coeff.n_gamma(cl_indx(i));
            q = T3_XSect_GammaS_Coeff.q_gamma(cl_indx(i));
            r = T3_XSect_GammaS_Coeff.r_gamma(cl_indx(i));
            p = T3_XSect_GammaS_Coeff.p_gamma(cl_indx(i));
            s = T3_XSect_GammaS_Coeff.s_gamma(cl_indx(i));
            % --- Calculate the equation constants
            Zadj = Z ./ (a.*hv.^b);
            % --- Calculate GAMMA
            gamma = (m .* sin((n.*Zadj.^p + q)) + s .* Zadj.^p + r).*sqrt(hv);
        else
            % --- Extract the coefficients
            da = T3_XSect_Gamma_Coeff.da_gamma(cl_indx(i)-5);
            ga = T3_XSect_Gamma_Coeff.ga_gamma(cl_indx(i)-5);
            dm = T3_XSect_Gamma_Coeff.dm_gamma(cl_indx(i)-5);
            gm = T3_XSect_Gamma_Coeff.gm_gamma(cl_indx(i)-5);
            du = T3_XSect_Gamma_Coeff.du_gamma(cl_indx(i)-5);
            gu = T3_XSect_Gamma_Coeff.gu_gamma(cl_indx(i)-5);
            cw = T3_XSect_Gamma_Coeff.cw_gamma(cl_indx(i)-5);
            dw  = T3_XSect_Gamma_Coeff.dw_gamma(cl_indx(i)-5);
            n  = T3_XSect_Gamma_Coeff.n_gamma(cl_indx(i)-5);
            cr = T3_XSect_Gamma_Coeff.cr_gamma(cl_indx(i)-5);
            dr = T3_XSect_Gamma_Coeff.dr_gamma(cl_indx(i)-5);
            p  = T3_XSect_Gamma_Coeff.p_gamma(cl_indx(i)-5);
            cs = T3_XSect_Gamma_Coeff.cs_gamma(cl_indx(i)-5);
            ds = T3_XSect_Gamma_Coeff.ds_gamma(cl_indx(i)-5);
            % --- Calculate the equation constants
            a_gamma = da.*(hv).^ga;
            m_gamma = dm.*(hv).^gm;
            u_gamma = du.*(hv).^gu;
            w_gamma = cw + dw.*(hv);
            n_gamma = n;
            Ebe = calc_be(element, ATOM_CL{cl_indx(i)}); 
            if isempty(Ebe); Ebe = 0; end
            if isnan(Ebe); Ebe = 0; end
            Ek = hv - Ebe;
            r_gamma = Ek ./ (cr + dr.*hv);
            p_gamma = p;
            s_gamma = Ek ./ (cs + ds.*hv);
            % --- Calculate GAMMA
            gamma = a_gamma + m_gamma.*sin((u_gamma .* (Z/100) + w_gamma)) + n_gamma .*exp(-1.*r_gamma) + p_gamma.*exp(-1.*(s_gamma));
            if size(gamma, 2) > 1;  gamma = gamma'; end
            % gamma(isnan(gamma)) = 0;
        end
        % --- Calculate beta
        xsect_gamma(:,end+1) = gamma;
    end
end

%% -- Plot for debugging
if plot_results == 1
    nCL     = length(cl_indx);
    cols    = lines(nCL);
    figure(); hold on; grid on;
    for i = 1:nCL
        plot(hv, xsect_gamma(:,i), 'o-', 'markersize', 2, 'markeredgecolor', cols(i,:), 'markerfacecolor', cols(i,:), 'color', cols(i,:)); 
    end
    legend(corelevel, 'location', 'eastoutside', 'FontSize', 6);
    text(0.45, 0.07, sprintf("%s (Z=%i)", element,Z),...
        'FontSize', 12, 'color', 'k', 'Units','normalized', 'FontWeight', 'bold', 'HorizontalAlignment','center');
    title('Cant 2022 - Photoionisation Parameter (gamma)');
    xlabel('$$ \bf Photon\ Energy\ [eV] $$', 'interpreter', 'latex');
    ylabel('$$ \bf \gamma $$', 'interpreter', 'latex');
    ax = gca; ax.YScale = 'linear'; ax.XScale = 'linear';
    yline(0, 'k-', 'LineWidth',1, 'HandleVisibility','off');
    axis([100, 15000, -1, 2.75]);
end
end