function xsect_beta = xsect_beta_cant2022(element, corelevel, hv, plot_results)
% xsect_beta = xsect_beta_cant2022(element, corelevel, hv, plot_results)
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
%   -   xsect_beta:     N×M vector of the photoelectron angular distribution asymmetry (dipole) parameter beta.

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
T2_XSect_BetaS_Coeff = table();
T2_XSect_BetaS_Coeff.Shell = {'1s1';'2s1';'3s1';'4s1';'5s1'};
T2_XSect_BetaS_Coeff.ca_beta = [2.00E+00;2.01E+00;2.00E+00;1.44E+00;1.34E+00];
T2_XSect_BetaS_Coeff.da_beta = [-1.03E-02;-1.04E-02;-2.73E-02;1.89E-02;1.61E-02];
T2_XSect_BetaS_Coeff.dm_beta = [6.17E-01;5.42E-01;4.72E-02;5.35E-01;6.27E-01];
T2_XSect_BetaS_Coeff.gm_beta = [1.58E-01;1.18E-01;9.01E-01;-5.01E-02;-1.92E-02];
T2_XSect_BetaS_Coeff.dq_beta = [1.78E+00;1.59E+00;2.39E-01;1.20E-03;1.19E-03];
T2_XSect_BetaS_Coeff.gq_beta = [2.74E-02;3.70E-02;2.32E-01;2.43E+00;2.30E+00];
T2_XSect_BetaS_Coeff.dn_beta = [0.00E+00;6.40E-03;1.72E-03;5.87E-05;1.41E-04];
T2_XSect_BetaS_Coeff.gn_beta = [0.00E+00;3.40E-02;-1.07E-01;1.59E+00;1.22E+00];
T2_XSect_BetaS_Coeff.ds_beta = [0.00E+00;9.29E+00;6.70E+00;1.04E+01;8.77E+00];
T2_XSect_BetaS_Coeff.gs_beta = [0.00E+00;-2.04E-01;-4.52E-02;-2.94E-01;-2.26E-01];
% -- Table of all other core-level coefficients
T2_XSect_Beta_Coeff = table();
T2_XSect_Beta_Coeff.Shell = {'2p1';'2p3';'3p1';'3p3';'4p1';'4p3';'5p1';'5p3'; '3d3';'3d5';'4d3';'4d5';'5d3';'5d5';'4f5';'4f7'};
T2_XSect_Beta_Coeff.ca_beta = [1.48E+00;1.55E+00;1.64E+00;2.21E+00;1.74E+00;1.75E+00;1.72E+00;1.77E+00;1.32E+00;1.24E+00;1.53E+00;1.41E+00;1.44E+00;1.46E+00;1.06E+00;1.05E+00];
T2_XSect_Beta_Coeff.da_beta =[-6.67E-03;0.00E+00;-4.40E-03;0.00E+00;-1.07E-03;0.00E+00;-2.35E-03;0.00E+00;0.00E+00;0.00E+00;0.00E+00;0.00E+00;-2.50E-05;-2.49E-05;0.00E+00;-4.43E-06];
T2_XSect_Beta_Coeff.dm_beta = [3.98E+01;3.49E+01;9.59E+00;5.04E+00;3.57E+00;5.75E+00;3.52E+00;3.25E+00;6.86E+00;8.76E+00;5.81E+00;4.88E+00;4.55E+00;2.51E+00;4.41E+00;4.10E+00];
T2_XSect_Beta_Coeff.gm_beta = [-8.14E-01;-7.85E-01;-5.60E-01;-2.87E-01;-3.56E-01;-3.87E-01;-3.54E-01;-1.86E-01;-2.70E-01;-5.74E-01;-3.17E-01;-6.01E-01;-1.82E-01;-4.04E-01;-6.23E-01;-6.43E-01];
T2_XSect_Beta_Coeff.db_beta = [2.24E-01;2.37E-01;3.56E-01;4.91E-01;4.76E-01;4.56E-01;5.20E-01;4.57E-01;4.67E-01;4.22E-01;5.28E-01;5.36E-01;6.06E-01;6.60E-01;6.20E-01;6.31E-01];
T2_XSect_Beta_Coeff.gb_beta = [4.57E-01;4.57E-01;3.89E-01;2.88E-01;3.48E-01;3.16E-01;3.12E-01;3.28E-01;3.40E-01;4.17E-01;3.38E-01;4.09E-01;2.73E-01;3.50E-01;3.95E-01;3.99E-01];
T2_XSect_Beta_Coeff.dq_beta = [1.11E-02;1.07E-02;0.00E+00;0.00E+00;0.00E+00;3.03E-04;3.06E-04;3.05E-04;5.21E-02;4.06E-02;5.01E-02;1.92E-02;3.96E-02;1.60E-10;6.37E-02;7.21E-02];
T2_XSect_Beta_Coeff.n_beta = [1.10E+00;1.31E+00;1.40E+00;2.90E+00;2.01E+00;8.94E-01;7.87E-01;7.83E-01;9.45E-01;7.36E-01;2.34E+00;2.93E+00;3.25E+00;2.92E+00;5.52E+00;5.60E+00];
T2_XSect_Beta_Coeff.cs_beta = [5.01E-01;5.25E-01;6.96E-01;3.72E-01;9.81E-01;1.41E+00;8.11E+00;8.20E+00;1.02E+00;1.03E+00;2.90E-01;3.96E-01;5.67E-01;4.81E-01;3.54E-01;3.26E-01];
T2_XSect_Beta_Coeff.ds_beta = [2.31E-01;3.31E-03;1.35E-01;9.16E-02;0.00E+00;0.00E+00;0.00E+00;0.00E+00;5.28E-01;4.45E-01;6.14E-01;3.12E-01;1.36E-02;0.00E+00;1.46E+00;1.34E+00];
T2_XSect_Beta_Coeff.cr_beta = [5.76E+01;6.09E+01;2.23E+01;7.58E+01;0.00E+00;0.00E+00;0.00E+00;0.00E+00;2.06E+02;2.12E+02;3.44E+02;1.39E+02;1.76E+02;1.30E+02;2.44E+02;2.06E+02];
T2_XSect_Beta_Coeff.dr_beta = [1.75E+01;2.51E+01;7.56E+01;1.37E+02;2.23E+02;2.27E+02;7.37E+02;7.12E+02;2.13E+02;1.34E+02;0.00E+00;0.00E+00;1.35E-03;4.44E-09;0.00E+00;0.00E+00];
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
%% - 4 - Calculating the dipole parameter, beta
cm22barn = 1e+24;
% -- Filing through all core-levels
xsect_beta = [];
for i = 1:length(cl_indx)
    if cl_indx(i) == 0; xsect_beta(:,i)   = NaN(size(hv));
    else
        beta = [];
        % --- For all S core levels, user the BetaS table of coefficients
        if cl_indx(i) < 6
            % --- Extract the coefficients
            ca = T2_XSect_BetaS_Coeff.ca_beta(cl_indx(i));
            da = T2_XSect_BetaS_Coeff.da_beta(cl_indx(i));
            dm = T2_XSect_BetaS_Coeff.dm_beta(cl_indx(i));
            gm = T2_XSect_BetaS_Coeff.gm_beta(cl_indx(i));
            dq = T2_XSect_BetaS_Coeff.dq_beta(cl_indx(i));
            gq = T2_XSect_BetaS_Coeff.gq_beta(cl_indx(i));
            dn = T2_XSect_BetaS_Coeff.dn_beta(cl_indx(i));
            gn = T2_XSect_BetaS_Coeff.gn_beta(cl_indx(i));
            ds = T2_XSect_BetaS_Coeff.ds_beta(cl_indx(i));
            gs = T2_XSect_BetaS_Coeff.gs_beta(cl_indx(i));
            % --- Calculate the equation constants
            a_beta = ca + da.*(hv./1000);
            m_beta = dm.*(hv./1000).^gm;
            q_beta = dq.*(hv./1000).^gq;
            n_beta = dn.*(hv./1000).^gn;
            s_beta = ds.*(hv./1000).^gs;
            % --- Calculate BETA
            beta = a_beta + m_beta.*(Z/100).^q_beta - n_beta.*exp(s_beta.*(Z/100));
        else
            % --- Extract the coefficients
            ca = T2_XSect_Beta_Coeff.ca_beta(cl_indx(i)-5);
            da = T2_XSect_Beta_Coeff.da_beta(cl_indx(i)-5);
            dm = T2_XSect_Beta_Coeff.dm_beta(cl_indx(i)-5);
            gm = T2_XSect_Beta_Coeff.gm_beta(cl_indx(i)-5);
            db = T2_XSect_Beta_Coeff.db_beta(cl_indx(i)-5);
            gb = T2_XSect_Beta_Coeff.gb_beta(cl_indx(i)-5);
            dq = T2_XSect_Beta_Coeff.dq_beta(cl_indx(i)-5);
            n  = T2_XSect_Beta_Coeff.n_beta(cl_indx(i)-5);
            cs = T2_XSect_Beta_Coeff.cs_beta(cl_indx(i)-5);
            ds = T2_XSect_Beta_Coeff.ds_beta(cl_indx(i)-5);
            cr = T2_XSect_Beta_Coeff.cr_beta(cl_indx(i)-5);
            dr = T2_XSect_Beta_Coeff.dr_beta(cl_indx(i)-5);
            % --- Calculate the equation constants
            a_beta = ca + da.*(hv./1000);
            m_beta = dm.*(hv./1000).^gm;
            b_beta = db.*(hv./1000).^gb;
            q_beta = dq.*(hv./1000);
            n_beta = n;
            s_beta = cs + ds ./(hv./1000);
            r_beta = cr + dr.*(hv./1000);
            Ebe = calc_be(element, ATOM_CL{cl_indx(i)}); 
            if isempty(Ebe); Ebe = 0; end
            if isnan(Ebe); Ebe = 0; end
            Ek = hv - Ebe;
            t_beta = (Ek./r_beta).^s_beta;
            % --- Calculate BETA
            beta = a_beta - m_beta.*((Z/100 - b_beta).^2).*(Z/100).^q_beta - n_beta.*exp(-1*t_beta);
            if size(beta, 2) > 1;  beta = beta'; end
            % beta(isnan(beta)) = 0;
        end   
        % --- Calculate beta
        xsect_beta(:,i) = beta;
    end
end
%% -- Plot for debugging
if plot_results == 1
    nCL     = length(cl_indx);
    cols    = lines(nCL);
    figure(); hold on; grid on;
    for i = 1:nCL
        plot(hv, xsect_beta(:,i), 'o-', 'markersize', 2, 'markeredgecolor', cols(i,:), 'markerfacecolor', cols(i,:), 'color', cols(i,:)); 
    end
    legend(corelevel, 'location', 'eastoutside', 'FontSize', 6);
    text(0.45, 0.07, sprintf("%s (Z=%i)", element,Z),...
        'FontSize', 12, 'color', 'k', 'Units','normalized', 'FontWeight', 'bold', 'HorizontalAlignment','center');
    title('Cant 2022 - Photoionisation Parameter (beta)');
    xlabel('$$ \bf Photon\ Energy\ [eV] $$', 'interpreter', 'latex');
    ylabel('$$ \bf \beta $$', 'interpreter', 'latex');
    ax = gca; ax.YScale = 'linear'; ax.XScale = 'linear';
    yline(0, 'k-', 'LineWidth',1, 'HandleVisibility','off');
    axis([100, 15000, -1, 2.75]);
end
end