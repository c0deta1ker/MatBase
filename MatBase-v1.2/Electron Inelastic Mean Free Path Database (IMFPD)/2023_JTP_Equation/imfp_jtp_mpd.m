function imfp = imfp_jtp_mpd(ke_dat, material)
% imfp = imfp_jtp_mpd(ke_dat, material)
%   Function that determines the electron inelastic mean free path (IMFP) 
%   based on the JTP method, which is a more accurate method than the 
%   standard universal curve and TPP-2M determination. In this function, 
%   you can define the material as a string and it will look it up the relevant 
%   parameters in the Material Properties Database (MPD) ('MPD_PCC.mat'). See
%   references [1] for more information.
%   [1] Jablonski A, Tanuma S, Powell CJ. Surf Interface Anal. 2023; 55(8): 609-637. doi:10.1002/sia.7217
%
%   IN:
%   -   ke_dat:     N×1 vector of the input electron kinetic energy (for PES; KE = HV - BE - PHI) [eV]
%   -   material:  	string of the material whose imfp is to be determined; e.g. "Si", "SiO2", "Al2O3"...
%
%   OUT:
%   -   imfp:       N×1 column vector of the electron IMFP values [Angstroms]

%% Default parameters (Parameters for Silicon)
if nargin < 2; material = "Si";  end
if isempty(material); material = "Si"; end
%% - 1 - Extracting the material parameters from the materials database
material = string(material);                % Ensure the input is a string
material_props = get_mpd_props(material);   % Extracting the material properties
% - Extracting the material properties required for the TPP-2M formalism
rho     = material_props.DENSITY;
Nv      = material_props.ELECT_VALENCY;
M       = material_props.ATOM_MASS;
Egap    = material_props.ELE_BGAP;
%% - 2 - Determination of the IMFP via TPP-2M formalism
imfp = imfp_jtp(ke_dat, rho, Nv, M, Egap);
end