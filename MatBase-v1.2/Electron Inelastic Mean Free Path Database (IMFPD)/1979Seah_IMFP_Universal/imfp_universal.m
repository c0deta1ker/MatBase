function imfp = imfp_universal(ke_dat)
% imfp = imfp_universal(ke_dat)
%   Function that determines the universal electron inelastic mean free
%   path (IMFP) in elements based on the KE^0.5 equation described by M. P. 
%   Seah [1]. This gives a good, first order approximation of the IMFP values, 
%   with the only required input being the electron kinetic energy. 
%   See reference [1] for more information.
%   [1] M. P. Seah, Quantitative electron spectroscopy of surfaces A Standard Data Base for Electron Inelastic Mean Free Paths in Solids (1979)
%
%   IN:
%   -   ke_dat:  	N×1 column vector of the input electron kinetic energy (for PES; KE = BE - PHI) [eV]
%
%   OUT:
%   -   imfp:       N×1 column vector of the electron IMFP values [Angstroms]

%% - 1 - Determination of the IMFP using universal formula
imfp = (143 ./ ke_dat.^2) + 0.054 .* sqrt(ke_dat);      % IMFP in nano-metres
imfp = 10 .* imfp;                                      % Convert IMFP from nm to Angstroms
%% -- Validity check on outputs
imfp(isnan(imfp)) = 0; if size(imfp, 2) >1; imfp = imfp'; end % Ensuring the imfp is a column vector
end