function xsect_gamma = calc_xsect_gamma(element, corelevel, hv, formalism, plot_results)
% sect_gamma = calc_xsect_gamma(element, corelevel, hv, formalism, plot_results)
%   This is a general function that calculates the cross sections for photoionization
%   from different sources in the literature. The Scofield (1973), Yeh &
%   Lindau (1985), Trzhaskovskaya & Yarzhemsky (2018) and Cant (2022) formalisms 
%   are available. The user can define a scalar or vector of incident photon energies and a
%   interpolation will be made between the dataset to determine the best
%   estimate of the cross sections for photoionization.
%
%   IN:
%   -   element:    	string of the element; e.g. "H", "He", "Si", "In"...
%   -   corelevel:      string or M×1 string vector of the core-levels to be probed; e.g. "5s1", "5p1", "5p3", "5d3", "5d5", "5f5', "5f7"...
%   -   hv:             scalar or N×1 vector of the incident photon energies [eV]
%   -   formalism:      string of the photoionization cross-section formalism to use. Default:"Cant2022" ["Scofield1973","YehLindau1985","Trzhaskovskaya2018","Cant2022"]
%   -   plot_results:   if 1, will plot figure summary, otherwise it wont.
%
%   OUT:
%   -   xsect_gamma:     N×M vector of the photoelectron angular distribution asymmetry (non-dipole) parameter gamma.
%
%   SEE REFERENCES:
%       [1] Scofield, J H. Theoretical photoionization cross sections from 1 to 1500 keV.. United States: N. p., 1973. Web. doi:10.2172/4545040.
%       [2] J.J. Yeh, I. Lindau. ATOMIC DATA AND NUCLEAR DATA TABLES 32, l-l 55 (1985). Web. doi:10.1016/0092-640X(85)90016-6
%       [3] M.B. Trzhaskovskaya, V.G. Yarzhemsky. Atomic Data and Nuclear Data Tables 119 (2018) 99–174. Web. doi:10.1016/j.adt.2017.04.003
%       [4] David J. H. Cant, Ben F. Spencer, Wendy R. Flavell, Alexander G. Shard. Surf Interface Anal. 2022; 54(4): 442-454. doi:10.1002/sia.7059

%% Default parameters
if nargin < 2; corelevel = [];  end
if nargin < 3; hv = 1e2:10:1.5e4;  end
if nargin < 4; formalism = "Cant2022";  end
if nargin < 5; plot_results = 0;  end
if isempty(corelevel); corelevel = []; end
if isempty(hv); hv = 1e2:10:1.5e4; end
if isempty(formalism); formalism = "Cant2022"; end
if isempty(plot_results); plot_results = 0; end
%% Validity checks on the input parameters
% -- Forcing the photon energy vector to have no duplicates and in ascending order
hv          = sort(unique(hv));
corelevel   = string(corelevel);
% -- Ensuring the inputs are in the form of 1D column vectors
if size(hv, 2) > 1;         hv = hv'; end
if size(corelevel, 2) > 1;  corelevel = corelevel'; end
%% - 1 - Determination of the photoionization cross section sigma
% -- Scofield1973 formalism
if strcmpi(formalism, "Scofield1973") || strcmpi(formalism, "Sco") || strcmpi(formalism, "S") || strcmpi(formalism, "1973")  || strcmpi(formalism, "S1973")
    xsect_gamma = []; msg = 'The non-dipole asymmetry parameter gamma is not available in the Scofield1973 database. Gamma is empty.'; 
    % disp(msg);
% -- YehLindau1985 formalism
elseif strcmpi(formalism, "YehLindau1985") || strcmpi(formalism, "Yeh") || strcmpi(formalism, "Lindau") || strcmpi(formalism, "YL") || strcmpi(formalism, "1985")  || strcmpi(formalism, "YL1985")
    xsect_gamma = []; msg = 'The non-dipole asymmetry parameter gamma is not available in the YehLindau1985 database. Gamma is empty.'; 
    % disp(msg);
% -- Trzhaskovskaya2018 formalism
elseif strcmpi(formalism, "Trzhaskovskaya2018") || strcmpi(formalism, "Trz") || strcmpi(formalism, "T") || strcmpi(formalism, "2018")  || strcmpi(formalism, "T2018")
    xsect_gamma = xsect_gamma_trzh2018(element, corelevel, hv, plot_results);
% -- Cant2020 formalism
elseif strcmpi(formalism, "Cant2022") || strcmpi(formalism, "Cant") || strcmpi(formalism, "C") || strcmpi(formalism, "2022")  || strcmpi(formalism, "C2022")
    xsect_gamma = xsect_gamma_cant2022(element, corelevel, hv, plot_results);
else; msg = 'Formalism not found. One of the following must be used: "Scofield1973", "YehLindau1985", "Trzhaskovskaya2018" or "Cant2022".'; error(msg);
end
end