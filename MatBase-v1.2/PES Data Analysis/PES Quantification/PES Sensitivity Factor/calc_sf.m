function [sf, params] = calc_sf(element, corelevel, material, hv, theta, phi, P, formalism_xsect, formalism_imfp)
% [sf, params] = calc_sf(element, corelevel, material, hv, theta, phi, P, formalism_xsect, formalism_imfp)
%   Function that determines the Sensitivity Factor (SF) in X-ray photoelectron
%   spectroscopy (XPS) using only fundamental parameters. The SFs in XPS are 
%   empirically derived factors (from compounds of known composition) 
%   by which the peak intensity is normalized (divided by) to provide an atomic 
%   concentration. The peak area measured for the "main" peak, I, is dependent 
%   on fundamental parameters from the equation below [1], but also on 
%   several experimental parameters specific to the instrument being used: 
%       I = (n)*(SF) = (n)*(JσλTS[1+sF]exp(-c/λ))
%   where:
%       n: Atomic concentration of the element.
%       σ: Cross section for a specific subshell.
%       λ: Inelastic mean free path (IMFP) for electrons with a given kinetic energy.
%       F: Parameter for angular anisotropy of electrons emitted from a subshell.
%       T = 1: Transmission function, describes intensity reduction of emitted electrons due to instrumental losses.
%       S = 1: Factor (0 - 1) for intensity reduction of emitted electrons due to elastic scattering, depends on element and emission angle.
%       s = 1: Factor (0.5 - 1) for change in angular distribution of emitted electrons due to elastic scattering, depends on element and emission angle.
%       c = 0: Contamination factor to correct for surface contamination attenuation.
%       J = 1: Photon flux incident on the sample.
%   Here, only the relevant σ, λ, F variables are used to calculate the SF.
%   SEE REFERENCES:
%       [1] David J. H. Cant, Ben F. Spencer, Wendy R. Flavell, Alexander G. Shard. Surf Interface Anal. 2022; 54(4): 442-454. doi:10.1002/sia.7059
%
%   IN:
%   -   element:    	    string of the element; e.g. "H", "He", "Si", "In"...
%   -   corelevel:          string of the core-level to be probed; e.g. "5s1", "5p1", "5p3", "5d3", "5d5", "5f5', "5f7"...
%   -   material:           char/string of the material; e.g. "Si", "InAs", "Al2O3"...
%   -   hv:                 scalar or N×1 vector of the incident photon energies [eV]
%   -   theta:              scalar or 1×M of the polar angle between the photoelectron vector relative to electric field vector (i.e. at normal emission: LV (p-pol, E//MP) = 0, LH (s-pol, E⊥MP) = 90) [degree]
%   -   phi:                scalar of the azimuthal angle between the photon momentum vector relative to the projection of the photoelectron vector on to the plane perpendicular to the electric field vector (i.e. normal emission = 0) [degree]
%   -   P:                  scalar of degree of polarization, where 1 or 0 is equivalent to full linear polarization, and 0.5 is equivalent to unpolarized light.
%   -   formalism_xsect:    string of the photoionization cross-section formalism_xsect to use. Default:"Cant2022" ["Scofield1973","YehLindau1985","Trzhaskovskaya2018","Cant2022"]
%   -   formalism_imfp:     string for imfp calculator formalism_xsect. Default:"S2" ["Universal","TPP2M","TPP2M-avg","Optical","S1","S2","S3","S3-organic","S4"]
%
%   OUT:
%   -   sf:      	        scalar or N×1 vector of the empirical XPS sensitivity factor 
%   -   params:      	    data structure containing all the parameters used in calculating the XPS sensitivity factor 

%% -- Validity check on inputs
if nargin < 9; formalism_imfp = "S2"; end
if nargin < 8; formalism_xsect = "Cant2022"; end
if nargin < 7; P = 0.5;  end
if nargin < 6; phi = 0;  end
if nargin < 5; theta = 0;  end
if isempty(formalism_imfp); formalism_imfp = "S2"; end
if isempty(formalism_xsect); formalism_xsect = "Cant2022"; end
if isempty(P); P = 0.5; end
if isempty(phi); phi = 0; end
if isempty(material); material = element; end
if isempty(theta); theta = 0; end
if size(hv, 2) >1; hv = hv'; end
if size(theta, 1) >1;   theta = theta'; end
%% 1 : Determination of the XPS Sensitivity Factor (SF)
% -- Extracting the XSECT
[sigma, beta, gamma, delta] = calc_xsect(element, corelevel, hv, formalism_xsect, 0);
if isempty(beta); beta = 0 .* ones(size(sigma)); end
if isempty(gamma); gamma = 0 .* ones(size(sigma)); end
if isempty(delta); delta = 0 .* ones(size(sigma)); end
% -- Extracting the IMFP
be          = calc_be(element,corelevel); 
if isempty(be); fprintf('Eb not found in any database, defaulting to 0 eV.\n'); be = 0; end
ke          = hv - be;
imfp        = calc_imfp(ke,formalism_imfp,material);
imfp        = imfp ./ 10;  % convert to nm
% -- Extracting the ANGULAR ANISOTROPY
F           = calc_angular_anisotropy(formalism_xsect, beta, gamma, delta, theta, phi, P, material);
% -- Calculating the SENSITIVITY FACTOR
sf          = sigma .* imfp .* (1 + F);
% sf          = sigma  .* (1 + F);
%% -- Validity check on outputs
sf = real(sf);
% if size(sf, 2) >1; sf = sf'; end
%% -- Saving all the params to a structure for debugging
params                      = struct();
params.hv                   = hv;
params.theta                = theta;
params.phi                  = phi;
params.P                    = P;
params.formalism_xsect      = formalism_xsect;
params.formalism_imfp       = formalism_imfp;
params.element              = element;
params.corelevel            = corelevel;
params.material             = material;
params.sigma                = sigma;
params.beta                 = beta;
params.gamma                = gamma;
params.delta                = delta;
params.be                   = be;
params.ke                   = ke;
params.imfp                 = imfp;
params.FP                   = F;
end