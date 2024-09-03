function FU = calc_angular_anisotropy_FU(beta, gamma, delta, omega)
% FU = calc_angular_anisotropy_FU(beta, gamma, delta, omega)
%   This function calculates the angular anisotropy factor, FU, for 
%   unpolarized light that is incident on the sample. FU depends on the 
%   experimental geometry, which is defined by spherical coordinates. 
%   The only angle needed for this function is the angle between the 
%   direction of the incoming photons and the direction of the emitted 
%   photoelectrons.
%   This formalism is common to use in the Hard X-Ray regime, where the
%   David Cant Cross-Sections are also used. This is from the original 
%   work of David J. H. Cant [1], see below.
%   [1] David J. H. Cant, Ben F. Spencer, Wendy R. Flavell, Alexander G. Shard. Surf Interface Anal. 2022; 54(4): 442-454. doi:10.1002/sia.7059
%
%   IN:
%   -   beta:  	    N×1 column vector of the dipole asymmetry factor.
%   -   gamma:      N×1 column vector of the first non-dipole asymmetry factor.
%   -   delta:      N×1 column vector of the second non-dipole asymmetry factor.
%   -   omega:      scalar or N×1 column vector of the angle between the photon momentum vector relative to the photoelectron vector [degree]
%
%   OUT:
%   -   FU:      	angular anisotropy factor for unpolarized light.

%% -- Validity check on inputs
if nargin < 4; omega = 90;  end
if isempty(omega);      omega = 90; end
if size(beta, 2) >1;    beta = beta'; end
if size(delta, 2) >1;   delta = delta'; end
if size(gamma, 2) >1;   gamma = gamma'; end
if size(omega, 2) >1; omega = omega'; end
%% 1 : Angular anisotropy factor for unpolarized light
Omega = deg2rad(omega);
FU = 1 - 0.25 .*beta .* (3 .* cos(Omega).^2 - 1) + (delta + 0.5 .* gamma .* sin(Omega).^2) .* cos(Omega);
%% -- Validity check on outputs
if size(FU, 2) >1; FU = FU'; end
end