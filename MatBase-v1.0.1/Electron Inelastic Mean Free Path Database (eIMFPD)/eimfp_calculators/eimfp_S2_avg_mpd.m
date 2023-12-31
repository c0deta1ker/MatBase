function imfp = eimfp_S2_avg_mpd(ke_dat, material)
% imfp = eimfp_S2_avg_mpd(ke_dat, material)
%   Function that determines the electron inelastic mean free path (IMFP) 
%   based on the S2 equation described by M. P. Seah [1]. The IMFP here
%   only depends on the value of Z. In this function, you can define the 
%   material as a string and it will look it up the relevant parameters in 
%   the Material Properties Database (MPD) ('MPD_PCC.mat') and determine 
%   the imfp using the S2 equation.
%
%   REFERENCE:
%   [1] M. P. Seah, “An accurate and simple universal curve for the 
%       energy-dependent electron inelastic mean free path,” 
%       Surf. Interface Anal., vol. 44, no. 4, pp. 497–503, 2012, 
%       doi: 10.1002/sia.4816.
%
%   IN:
%   -   ke_dat:     N×1 vector of the input electron kinetic energy (for PES; KE = BE - PHI) [eV]
%   -   material:  	string of the material whose imfp is to be determined; e.g. "Si", "SiO2", "InAs", "Al2O3"...
%
%   OUT:
%   -   imfp:       N×1 column vector of the electron IMFP values [Angstroms]

%% Default parameters (Parameters for Silicon)
if nargin < 2; material = "Si";  end
if isempty(material); material = "Si"; end
%% - 1 - Extracting the material parameters from the materials database
% - Extracting the material properties
material_props = get_mpd_props(material);
% - Extracting the material properties required for the S2 formalism
Z       = material_props.ATOM_ZNUM;
%% - 2 - Determination of the IMFP TPP-2M
% If the kinetic energy is negative, assume it is zero
ke_dat(ke_dat<0) = 0;
% Evaluating the imfp
a           = 2.5;  % Average atomic spacing in Angstroms
% Relativistic corrected equation
imfp = a .* ((1.52 + 0.167 .* Z.^(0.5) + 0.0394 .* ke_dat.^(0.872)) ./ Z.^(0.3)); % IMFP in Angstroms
% If isnan, return zero
imfp(isnan(imfp)) = 0;
%% Ensuring the imfp is a column vector
if size(imfp, 2) >1; imfp = imfp'; end
end