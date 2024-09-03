function eal = eal_S3_organic(ke_dat)
% eal = eal_S3_organic(ke_dat)
%   Function that determines the electron effective attenuation length (EAL) 
%   based on S3 formalism described by M. P. Seah [1] for organic materials.
%   This approximation is good for all organic matter.
%   See reference [1] for more information.
%   [1] M. P. Seah, Simple universal curve for the energy‐dependent electron attenuation length (2012)
%
%   IN:
%   -   ke_dat:     N×1 vector of the input electron kinetic energy (for PES; KE = BE - PHI) [eV]
%
%   OUT:
%   -   eal:        N×1 column vector of the electron eal values [Angstroms]

%% - 1 - Determination of the Effective Attenuation Length (ELA) based on S3 formalism
eal = 0.00837 .* ke_dat.^(0.842);   % IMFP in nm
eal = 10 .* eal;                    % Convert IMFP to Angstroms
%% -- Validity check on outputs
eal(isnan(eal)) = 0; if size(eal, 2) >1; eal = eal'; end % Ensuring the eal is a column vector
end