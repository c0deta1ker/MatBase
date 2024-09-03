function eal = eal_S3(ke_dat, rho, M, Z, stoic)
% eal = eal_S3(ke_dat, rho, M, Z, stoic)
%   Function that determines the electron effective attenuation length (EAL) 
%   based on S3 formalism described by M. P. Seah [1]. The input parameters 
%   of the materials are now required and this formalism is compatible
%   for elemental or binary materials. See reference [1] for more information.
%   [1] M. P. Seah, Simple universal curve for the energy‐dependent electron attenuation length (2012)
%
%   IN:
%   -   ke_dat:     N×1 vector of the input electron kinetic energy (for PES; KE = BE - PHI) [eV]
%   -   rho:        scalar of the density of the material (g/cc)
%   -   M:          scalar of the atomic or molecular weight (in amu == g/mol)
%   -   Z:          scalar of the atomic mass number of the element (or average for compound) (Z) 
%   -   stoic:      scalar of the stoiciometry of the material (e.g. for elements, stoic = 1, for molecules of the form G_gH_h, the stoiciometry is g + h)
%
%   OUT:
%   -   eal:        N×1 column vector of the electron eal values [Angstroms]

%% Default parameters (Parameters for Silicon)
if nargin < 5; stoic = 1; end
if nargin < 4; Z = 14; end
if nargin < 3; M = 28.085000;   end
if nargin < 2; rho = 2.3300000; end
if isempty(Z);      Z = 14; end
if isempty(M);      M = 28.085000; end      % amu == g/mol
if isempty(rho); 	rho = 2.3300000; end    % g/cc
if isempty(stoic); 	stoic = 1; end
% -- Defining the constants
pc = physics_constants(); NA  = pc.NA;    % Avogadro constant
%% - 1 - Determination of the Effective Attenuation Length (EAL) based on S4 formalism
% Extracting parameters for EAL
a = ((1e21 .* M) ./ (rho .* NA .* stoic)).^(1./3); % nm
W = 0;  % for elements
% Calculating the EAL
eal = (5.8 + 0.0041 .*Z.^(1.7) + 0.088 .* ke_dat.^(0.93)) .* a.^(1.82) ./ (Z.^(0.38) .* (1-W));
eal = 10 .* eal;                                      % Convert eal from nm to Angstroms
%% -- Validity check on outputs
if rho == 0; eal(:) = Inf; end
eal(isnan(eal)) = 0; if size(eal, 2) >1; eal = eal'; end % Ensuring the eal is a column vector
end