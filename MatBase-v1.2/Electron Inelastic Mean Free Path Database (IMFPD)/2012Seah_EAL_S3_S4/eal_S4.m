function eal = eal_S4(ke_dat, Z)
% eal = eal_S4(ke_dat, Z)
%   Function that determines the electron effective attenuation length (EAL) 
%   based on S4 formalism described by M. P. Seah [1]. The EAL here
%   only depends on the value of Z. See reference [1] for more information.
%   [1] M. P. Seah, Simple universal curve for the energy‐dependent electron attenuation length (2012)
%
%   IN:
%   -   ke_dat:     N×1 vector of the input electron kinetic energy (for PES; KE = BE - PHI) [eV]
%   -   Z:          scalar of the atomic mass number (Z) (or average for compounds)
%
%   OUT:
%   -   eal:        N×1 column vector of the electron eal values [Angstroms]

%% Default parameters (Parameters for Silicon)
if nargin < 2;  Z = 14; end
if isempty(Z); 	Z = 14; end
%% - 1 - Determination of the Effective Attenuation Length (EAL) based on S4 formalism
eal = (0.65 + 0.007.*ke_dat.^(0.93)) ./ Z.^(0.38);      % EAL in nm
eal = 10 .* eal;                                        % Convert to Angstroms
%% -- Validity check on outputs
eal(isnan(eal)) = 0; if size(eal, 2) >1; eal = eal'; end % Ensuring the eal is a column vector
end