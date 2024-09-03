function pes_model = nlayer_model01_run(lyr_mat, lyr_thick, lyr_ele, lyr_cls, hv, theta, phi, P, formalism_xsect, formalism_imfp, plot_result)
% pes_model = nlayer_model01_run(lyr_mat, lyr_thick, lyr_ele, lyr_cls, hv, theta, phi, P, formalism_xsect, formalism_imfp, plot_result)
%   This function calculates the total photoelectron intensity originating
%   from N independent layers in a sample. The sample is composed of various
%   materials, each with a user-defined thickness. The stack is made up of N layers,
%   each with a specified material and thickness. The layers are defined using a 
%   top-down approach, starting from the surface and ending with the bulk. 
%   'lyr_mat' is used to define the materials for each layer, and 'lyr_thick' 
%   is used to specify the thickness of each layer. To define the bulk layer, 
%   the final 'lyr_thick' value should be set to 'Inf'. The function also 
%   requires the definition of the core-level being probed and the experimental 
%   geometry. The formalism of the cross-sections and inelastic mean free
%   path (IMFP) can also be defined.
%   --------------------------------------------------------------------
%   The model does the following:
%
%       (1)  Using the input photon energies, the electron kinetic energies
%       are determined so that the attenuation length / inelastic-mean free
%       path (IMFP) of electrons in each layer can be determined. 
%       Using:
%           Ek = Ehv - Ebe ,
%       Ek is the electron kinetic energy; Ehv is the input photon energy;
%       Ebe is the binding energy of the electron, which is determined by
%       in the Photoionisation Energy and Fluorescence Database (PIEFD)
%       whose value depends on the input core-level defined; note, the sample
%       work function is omitted as it leads to only a very small change to 
%       the electron energy and is not well defined for all materials.
%
%       (2) The total intensity of photoelectrons from each layer can be 
%       determined via:
%           I = (n)*(SF) = (n)*(σλ[1+F])
%       where n is the atomic concentration of the element concerned; σ is the 
%       cross section for a specific subshell; λ is the inelastic mean free
%       path (IMFP) for electrons with a given kinetic energy; F is a parameter 
%       describing the angular anisotropy of electrons emitted from a given
%       subshell.
%       (2.1) The atomic concentration (n) can be determined by:
%               n = N_AVOGADRO * DENSITY / MOLECULAR WEIGHT.
%       The density and molecular weight are both contained within the MPD
%       for each material defined.
%       (2.2) The Sensitivity Factor (SF) is determined using the
%       'calc_sf()' function, which takes into account the photoionisation
%       cross-sections & asymmetry parameters, for a given photon energy
%       and experimental geometry. The transmission and surface containation 
%       correction factors are neglected in this model.
%       Thus, the product of n and SF yields the TOTAL intensity I0 of each 
%       layer.
%
%       (3) The Beer-Lambert law is used to model the attenuation of the
%       photoelectron intensity as a function of depth. This model states
%       that the photoelectron intensity decreases exponentially due to 
%       inelastic scattering, with a decay constant equal to the IMFP.
%       Thus, for a bulk, single-layer system, the Beer-Lambert law states:
%           I(z) = I0 * exp(-z / IMFP),
%       which gives the photoelectron intensity, I, that emerges from a
%       depth, z, within the sample. This is applied to each layer in the
%       n-layered stack defined by the user, where I0 is determined from
%       part (2), above. The total integral intensity for each
%       layer in the sample stack then gives the total emitted
%       photoelectron intensity.
%           
%       (4) Once the integral intensites from each layer are determined,
%       the final intensities are expressed in terms of a 'relative
%       contribution'. Here, each layer intensity is normalised by the sum
%       of all the intensities at each photon energy.
%   --------------------------------------------------------------------
%
%   IN:
%   -   lyr_mat:  	        Mx1 cell-vector of the material for each layer in the stack; e.g. "Si", "SiO2", "Al2O3"...
%   -   lyr_thick:	        Mx1 cell-vector of the thickness of each layer in the stack (nm)
%   -   lyr_ele:            Mx1 cell-vector of strings of the element to be probed; e.g. "H", "He", "Si", "In"...
%   -   lyr_cls:  	        Mx1 cell-vector of strings of the core-level to be probed; e.g. "5s1", "5p1", "5p3", "5d3", "5d5", "5f5', "5f7"...
%   -   hv:                 scalar or N×1 vector of the incident photon energies [eV]
%   -   theta:              scalar or 1×L vector of the polar angle between the photoelectron vector relative to electric field vector (i.e. at normal emission: LV (p-pol, E//MP) = 0, LH (s-pol, E⊥MP) = 90) [degree]
%   -   phi:                scalar of the azimuthal angle between the photon momentum vector relative to the projection of the photoelectron vector on to the plane perpendicular to the electric field vector (i.e. normal emission = 0) [degree]
%   -   P:                  scalar of degree of polarization, where 1 or 0 is equivalent to full linear polarization, and 0.5 is equivalent to unpolarized light.
%   -   formalism_xsect:    string of the photoionization cross-section formalism to use. Default:"Cant2022" ["Scofield1973","YehLindau1985","Trzhaskovskaya2018","Cant2022"]
%   -   formalism_imfp:     string for imfp calculator formalism. Default:"S2" ["Universal","Optical","TPP2M","TPP2M-avg","S1","S2","S3","S4"]
%   -   plot_results:       if 1, will plot figure summary, otherwise it wont.
%
%   OUT:
%   -   pes_model:          data structure that contains all the pes model parameters and variables.

%% Default parameters
% -- Default Inputs
if nargin < 11; plot_result = 0; end
if nargin < 10; formalism_imfp = "JTP"; end
if nargin < 9; formalism_xsect = "Cant2022"; end
if nargin < 8; P = 0.5;  end
if nargin < 7; phi = 0;  end
if nargin < 6; theta = 0;  end
if isempty(plot_result); plot_result = 0; end
if isempty(formalism_imfp); formalism_imfp = "JTP"; end
if isempty(formalism_xsect); formalism_xsect = "Cant2022"; end
if isempty(P); P = 0.5; end
if isempty(phi); phi = 0; end
if isempty(theta); theta = 0; end
%% Validity check on inputs
% -- Ensuring the user-defined material input is a cell-array
if ~iscell(lyr_mat);    lyr_mat = cellstr(lyr_mat); end
if ~iscell(lyr_thick);  lyr_thick = num2cell(lyr_thick); end
if ~iscell(lyr_ele);    lyr_ele = cellstr(lyr_ele); end
if ~iscell(lyr_cls);    lyr_cls = cellstr(lyr_cls); end
% -- Ensuring the photon energy and emission angle are the correct size
if size(hv, 2) > 1; hv = hv'; end
if size(theta, 1) > 1; theta = theta'; end
% -- Ensuring the formalisms are strings
formalism_xsect = string(formalism_xsect);
formalism_imfp  = string(formalism_imfp);
%% Defining constants
Nlyrs       = length(lyr_mat);
% -- Verify that the input layer variables are consistent in size
lyr_thick   = lyr_thick(1:Nlyrs);
lyr_ele     = lyr_ele(1:Nlyrs);
lyr_cls   	= lyr_cls(1:Nlyrs);
%% 1    :   Determine the maximum intensity I0 for each layer material
for i = 1:Nlyrs
    % - Extracting all relevant material parameters
    lyr_mpd{i}	    = get_mpd_props(lyr_mat{i});
    lyr_density{i}	= lyr_mpd{i}.DENSITY;
    lyr_atmass{i}   = lyr_mpd{i}.ATOM_MASS;
    lyr_atznum{i}   = lyr_mpd{i}.ATOM_ZNUM;
    % - Determination of the atomic concentration (n)
    pc = physics_constants();
    lyr_atstoic{i}  = get_ele_ratio_from_mat(lyr_mat{i}, lyr_ele{i});              % number of atoms per molecule (scales the intensity according the ratio of elements in the compound)
    lyr_n{i}        = lyr_atstoic{i} .* pc.NA.*(lyr_density{i} ./ lyr_atmass{i});   % number of atoms per unit volume
    % - Determination of the Sensitivity Factor (SF)
    [lyr_sf{i}, lyr_sf_params{i}] = calc_sf(lyr_ele{i}, lyr_cls{i}, lyr_mat{i}, hv, theta, phi, P, formalism_xsect, formalism_imfp);
    % - Determine the value of I0
    lyr_I0{i}   = lyr_n{i} .* lyr_sf{i};
end
%% 2    :   Determination of the IMFP at all kinetic energy for each layer material
for i = 1:Nlyrs
    lyr_be{i}       = calc_be(lyr_ele{i}, lyr_cls{i}); 
    lyr_ke{i}       = hv - lyr_be{i};
    lyr_imfp{i}     = calc_imfp(lyr_ke{i}, formalism_imfp, lyr_mat{i});
    lyr_imfp{i}     = 0.1 .* lyr_imfp{i}; % Convert the IMFP into nm
end
%% 3    :   Use Beer-Lambert law to determine the total photoelectron intensity for each layer
% Defining the Beer-Lambert Law of Attenuation for each layer
Theta = deg2rad(theta);
lyr_I = {}; lyr_Inorm = {};
% -- For one layer (Bulk)
if Nlyrs == 1;  lyr_I{1} = lyr_I0{1} .* exp(-lyr_thick{1} ./ (lyr_imfp{1} .* cos(Theta)));
% -- For multi-layer samples
else
    for i = 1:Nlyrs
        % -- For the uppermost layer, there is only attenuation through 1 layer
        if i == 1; lyr_I{i} = lyr_I0{i} .* (1 - exp(-lyr_thick{i} ./ (lyr_imfp{i} .* cos(Theta))));
        % -- For the lowermost bulk layer, there is attenuation through all layers above it
        elseif i == Nlyrs
            lyr_I{i} = lyr_I0{i};
            for j = 1:Nlyrs - 1; lyr_I{i} = lyr_I{i} .* (exp(-lyr_thick{j} ./ (lyr_imfp{j} .* cos(Theta)))); end
        % -- For intermediate layers
        else
            lyr_I{i} = lyr_I0{i} .* (1 - exp(-lyr_thick{i} ./ (lyr_imfp{i} .* cos(Theta))));
            for j = 1:i-1; lyr_I{i} = lyr_I{i} .* (exp(-lyr_thick{j} ./ (lyr_imfp{j} .* cos(Theta)))); end
        end
    end
end
% Normalising the intensities to the total sums
lyr_norm = lyr_I{1};
for i = 2:Nlyrs; lyr_norm = lyr_norm + lyr_I{i}; end
for i = 1:Nlyrs; lyr_Inorm{i} = lyr_I{i} ./ lyr_norm; end
%% 4    :   Creating a MATLAB data structure to store all of the results
% - Defining the structure
pes_model               = struct();
% - Defining the input arguments
pes_model.lyr_mat       = lyr_mat;
pes_model.lyr_thick     = lyr_thick;
pes_model.lyr_ele       = lyr_ele;
pes_model.lyr_cls       = lyr_cls;
pes_model.hv            = hv;
pes_model.theta         = theta;
pes_model.phi           = phi;
pes_model.P             = P;
pes_model.formalism_xsect   = formalism_xsect;
pes_model.formalism_imfp    = formalism_imfp;
% - Determination of the Initial Photoelectron Intensity
pes_model.Nlyrs             = Nlyrs;
pes_model.lyr_density       = lyr_density;
pes_model.lyr_atmass        = lyr_atmass;
pes_model.lyr_atznum        = lyr_atznum;
pes_model.lyr_atstoic       = lyr_atstoic;
pes_model.lyr_n             = lyr_n;
pes_model.lyr_sf            = lyr_sf;
pes_model.lyr_sf_params     = lyr_sf_params;
pes_model.lyr_I0            = lyr_I0;
% - Defining the properties for extracting the IMFP
pes_model.lyr_be        = lyr_be;
pes_model.lyr_ke	    = lyr_ke;
pes_model.lyr_imfp      = lyr_imfp;
% - Defining the properties for extracting value of I0
pes_model.lyr_I0        = lyr_I0;
% - Defining the final intensities
pes_model.lyr_norm    	= lyr_norm;
pes_model.lyr_I         = lyr_I;
pes_model.lyr_Inorm  	= lyr_Inorm;
%% 5    :   Plotting the multilayer model stack and solutions
if plot_result == 1; nlayer_model01_view(pes_model); end
end