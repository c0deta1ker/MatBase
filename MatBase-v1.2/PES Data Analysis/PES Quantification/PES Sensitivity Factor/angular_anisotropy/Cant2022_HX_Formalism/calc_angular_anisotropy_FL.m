function FL = calc_angular_anisotropy_FL(beta, gamma, delta, theta, phi)
% FL = calc_angular_anisotropy_FL(beta, gamma, delta, theta, phi)
%   This function computes the angular anisotropy factor, FL, for 
%   linearly polarized light that is incident on the sample. FL depends on the 
%   experimental geometry, which is specified by two angles: theta and phi. 
%   Theta is the angle between the emitted photoelectrons and the electric field 
%   vector, and phi is the angle between the direction of the incident photons 
%   and the direction of the emitted photoelectrons. 
%   This formalism is common to use in the Hard X-Ray regime, where the
%   David Cant Cross-Sections are also used. This is from the original 
%   work of David J. H. Cant [1], see below.
%   [1] David J. H. Cant, Ben F. Spencer, Wendy R. Flavell, Alexander G. Shard. Surf Interface Anal. 2022; 54(4): 442-454. doi:10.1002/sia.7059
%
%   IN:
%   -   beta:  	    N×1 column vector of the dipole asymmetry factor.
%   -   gamma:      N×1 column vector of the first non-dipole asymmetry factor.
%   -   delta:      N×1 column vector of the second non-dipole asymmetry factor.
%   -   theta:      scalar or N×1 column vector of the polar angle between the photoelectron vector relative to electric field vector (i.e. at normal emission: LV (p-pol, E//MP) = 0, LH (s-pol, E⊥MP) = 90) [degree]
%   -   phi:        scalar of the azimuthal angle between the photon momentum vector relative to the projection of the photoelectron vector on to the plane perpendicular to the electric field vector (i.e. normal emission = 0) [degree]
%
%   OUT:
%   -   FL:      	angular anisotropy factor for linearly polarized light.

%% -- Validity check on inputs
if nargin < 4; theta = 0;  end
if nargin < 5; phi = 0;  end
if isempty(theta);      theta = 0; end
if isempty(phi);        phi = 0; end
if size(beta, 2) >1;    beta = beta'; end
if size(gamma, 2) >1;   gamma = gamma'; end
if size(delta, 2) >1;   delta = delta'; end
if size(theta, 2) >1; theta = theta'; end
%% 1 : Angular anisotropy factor for linearly polarized light
FL = 1 + 0.5 .* beta .* (3 .* cos(deg2rad(theta)).^2 - 1) + (delta + gamma .* cos(deg2rad(theta)).^2) .* sin(deg2rad(theta)) .* cos(deg2rad(phi));
%% -- Validity check on outputs
if size(FL, 2) >1; FL = FL'; end
end