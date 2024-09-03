function imfp = imfp_S1(ke_dat, rho, M, Egap, Z, stoic)
% imfp = imfp_S1(ke_dat, rho, M, Egap, Z, stoic)
%   Function that determines the electron inelastic mean free path (IMFP) 
%   based on the S1 equation described by M. P. Seah [1]. The input parameters 
%   of the materials are now required and this formalism is compatible
%   for elemental or binary materials. See reference [1] for more information.
%   [1] M. P. Seah, An accurate and simple universal curve for the energy-dependent electron inelastic mean free path (2011)
%   
%   IN:
%   -   ke_dat:     N×1 vector of the input electron kinetic energy (for PES; KE = BE - PHI) [eV]
%   -   rho:        scalar of the density of the material (g/cc)
%   -   M:          scalar of the atomic or molecular weight (in amu == g/mol)
%   -   Egap:       scalar of the band gap energy (eV)
%   -   Z:          scalar of the atomic mass number of the element (or average for compound) (Z) 
%   -   stoic:      scalar of the stoiciometry of the material (e.g. for elements, stoic = 1, for molecules of the form G_gH_h, the stoiciometry is g + h)
%
%   OUT:
%   -   imfp:       N×1 column vector of the electron IMFP values [Angstroms]

%% Default parameters (Parameters for Silicon)
if nargin < 6; stoic = 1; end
if nargin < 5; Z = 14; end
if nargin < 4; Egap = 1.12000; end
if nargin < 3; M = 28.085000;   end
if nargin < 2; rho = 2.3300000; end
if isempty(Z);      Z = 14; end
if isempty(Egap);   Egap = 1.12000; end     % eV
if isempty(M);      M = 28.085000; end      % amu == g/mol
if isempty(rho); 	rho = 2.3300000; end    % g/cc
if isempty(stoic); 	stoic = 1; end
% -- Defining the constants
pc = physics_constants(); NA  = pc.NA;    % Avogadro constant
%% - 1 - Determination of the IMFP using S1 formalism (M P Seah)
% Extracting parameters for IMFP
a           = ((1e21 .* M) ./ (rho .* NA .* stoic)).^(1./3); % nm
% Calculating the IMFP
imfp = (4 + 0.44 .* Z.^(0.5) + 0.104 .* ke_dat.^(0.872)) .* a.^(1.7) ./ (Z.^(0.3) .* (1 - 0.02 .*Egap)); % IMFP in nanometres
imfp = 10 .* imfp;                                      % Convert IMFP from nm to Angstroms
%% -- Validity check on outputs
if rho == 0; imfp(:) = Inf; end
imfp(isnan(imfp)) = 0; if size(imfp, 2) >1; imfp = imfp'; end % Ensuring the imfp is a column vector
end