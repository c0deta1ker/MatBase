function FP = calc_angular_anisotropy(formalism, beta, gamma, delta, theta, phi, P, material)
% FP = calc_angular_anisotropy(formalism, beta, gamma, delta, theta, phi, P, material)
%   This function computes the angular anisotropy factor, FL, for 
%   linearly polarized light that is incident on the sample. FL depends on the 
%   experimental geometry, which is specified by two angles: theta and phi. 
%   Theta is the angle between the emitted photoelectrons and the electric field 
%   vector, and phi is the angle between the direction of the incident photons 
%   and the direction of the emitted photoelectrons. 
%   This formalism is common to use in the Hard X-Ray regime, where the
%   David Cant Cross-Sections are also used. This is from the original 
%   work of David J. H. Cant [1], see below.
%   [1] David J. H. Cant, Ben FP. Spencer, Wendy R. Flavell, Alexander G. Shard. Surf Interface Anal. 2022; 54(4): 442-454. doi:10.1002/sia.7059
%
%   IN:
%   -   formalism:  string of the angular asymmetry formalism to use. Default:"Cant2022" ["Scofield1973","YehLindau1985","Trzhaskovskaya2018","Cant2022"]
%   -   beta:  	    N×1 column vector of the dipole asymmetry factor.
%   -   gamma:      N×1 column vector of the first non-dipole asymmetry factor.
%   -   delta:      N×1 column vector of the second non-dipole asymmetry factor.
%   -   theta:      scalar or N×1 column vector of the polar angle between the photoelectron vector relative to electric field vector (i.e. at normal emission: LV (p-pol, E//MP) = 0, LH (s-pol, E⊥MP) = 90) [degree]
%   -   phi:        scalar of the azimuthal angle between the photon momentum vector relative to the projection of the photoelectron vector on to the plane perpendicular to the electric field vector (i.e. normal emission = 0) [degree]
%   -   P:          scalar of degree of polarization, where 1 or 0 is equivalent to full linear polarization, and 0.5 is equivalent to unpolarized light.
%   -   material:   char/string of the material beng probed; e.g. "Si", "InAs", "Al2O3"... This is used to determine the average Z number for the Scofield1973 &YehLindau1985 formalisms.
%
%   OUT:
%   -   FP:      	N×1 column vector of the angular anisotropy factor

%% Default parameters
if nargin < 5; theta = 0; end
if nargin < 6; phi = 0; end
if nargin < 7; P = 0.5; end
if nargin < 8; material = []; end
if isempty(theta); theta = 0; end
if isempty(phi); phi = 0; end
if isempty(P); P = 0.5; end
if isempty(material); material = []; end
%% -- Validity check on inputs
if size(theta, 1) >1;   theta = theta'; end
%% - 1 - Determination of the angular anisotropy factor
% -- Scofield1973 formalism
if strcmpi(formalism, "Scofield1973") || strcmpi(formalism, "Sco") || strcmpi(formalism, "S") || strcmpi(formalism, "1973")  || strcmpi(formalism, "S1973")
    material = string(material); 
    material_props = get_mpd_props(material);
    avg_Z = material_props.ATOM_ZNUM;
    FP = calc_angular_anisotropy_Fadj(beta, theta, avg_Z);
% -- YehLindau1985 formalism
elseif strcmpi(formalism, "YehLindau1985") || strcmpi(formalism, "Yeh") || strcmpi(formalism, "Lindau") || strcmpi(formalism, "YL") || strcmpi(formalism, "1985")  || strcmpi(formalism, "YL1985")
    material = string(material); 
    material_props = get_mpd_props(material);
    avg_Z = material_props.ATOM_ZNUM;
    FP = calc_angular_anisotropy_Fadj(beta, theta, avg_Z);
% -- Trzhaskovskaya2018 formalism
elseif strcmpi(formalism, "Trzhaskovskaya2018") || strcmpi(formalism, "Trz") || strcmpi(formalism, "T") || strcmpi(formalism, "2018")  || strcmpi(formalism, "T2018")
    FP = calc_angular_anisotropy_FP(beta, gamma, delta, theta, phi, P);
% -- Cant2020 formalism
elseif strcmpi(formalism, "Cant2022") || strcmpi(formalism, "Cant") || strcmpi(formalism, "C") || strcmpi(formalism, "2022")  || strcmpi(formalism, "C2022")
    FP = calc_angular_anisotropy_FP(beta, gamma, delta, theta, phi, P);
else; msg = 'Formalism not found. One of the following must be used: "Scofield1973", "YehLindau1985", "Trzhaskovskaya2018" or "Cant2022".'; error(msg);
end
%% Validity checks on the output parameters
FP = real(FP);
end