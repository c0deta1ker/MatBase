function Fadj = calc_angular_anisotropy_Fadj(beta, omega, avg_Z)
% Fadj = calc_angular_anisotropy_Fadj(beta, omega, avg_Z)
%   This function computes the angular anisotropy factor, Fadj, for 
%   unpolarized light that is incident on the sample. Fadj depends on the 
%   experimental geometry, which is defined by spherical coordinates. 
%   The only angle needed for this function is the angle between the 
%   direction of the incoming photons and the direction of the emitted 
%   photoelectrons. First, the adjusted asymmetry factor for solid
%   materials is determined via the equation:
%       β* = β (a - bZ + cZ^2), where a = 0.781; b = 0.00514; c = 0.000031;
%   This is then used to determine the angular anisotropy via the equation:
%       L = 1 + 1/2 β* (3/2 sin^2(ω) - 1)
%   This formalism is common to use in the Soft X-Ray regime, where the
%   Yeh-Lindau Cross-Sections are also used. See reference below for more
%   information.
%   [1] Angle Resoled XPS, Thermo Scientific, Application Note: 31014
%
%   IN:
%   -   beta:  	    N×1 column vector of the dipole asymmetry factor.
%   -   omega:      scalar of the angle between the photon momentum vector relative to the photoelectron vector [degree]
%   -   avg_Z:      scalar of the aveage mass number of the compound.
%
%   OUT:
%   -   Fadj:      	angular anisotropy factor for unpolarised light.

%% -- Validity check on inputs
if size(beta, 2) > 1;    beta = beta'; end
if isnan(beta);  beta = 2; end
%% 1 : Determine the adjusted asymmetry factor for solid materials
a = 0.781; b = 0.00514; c = 0.000031;
beta_adj = beta .* (a - b.*avg_Z + c.*avg_Z.^2);    
%% 2 : Determine the angular correction factor using the asymmetry parameter
Fadj = 1 + 0.5.*beta_adj .* (1.5.*(cos(deg2rad(omega))).^2 - 1);
%% -- Validity check on outputs
% if size(Fadj, 2) >1; Fadj = Fadj'; end
end