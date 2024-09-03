function pes_model = nlayer_model02_run(lyr_mat, lyr_conc, lyr_ele, lyr_cls, lyr_thick, hv, theta, phi, P, formalism_xsect, formalism_imfp, plot_result)
% pes_model = nlayer_model02_run(lyr_mat, lyr_conc, lyr_ele, lyr_thick, lyr_cls, hv, theta, phi, P, formalism_xsect, formalism_imfp, plot_result)
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
%   -   lyr_ele:            Mx1 cell-vector of strings of the element to be probed; e.g. "H", "He", "Si", "In"...
%   -   lyr_cls:  	        Mx1 cell-vector of strings of the core-level to be probed; e.g. "5s1", "5p1", "5p3", "5d3", "5d5", "5f5', "5f7"...
%   -   lyr_thick:	        Mx1 cell-vector of the thickness of each layer in the stack (nm)
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
if nargin < 12; plot_result = 0; end
if nargin < 11; formalism_imfp = "S3"; end
if nargin < 10; formalism_xsect = "Cant2022"; end
if nargin < 9; P = 0.5;  end
if nargin < 8; phi = 0;  end
if nargin < 7; theta = 0;  end
if isempty(plot_result); plot_result = 0; end
if isempty(formalism_imfp); formalism_imfp = "S3"; end
if isempty(formalism_xsect); formalism_xsect = "Cant2022"; end
if isempty(P); P = 0.5; end
if isempty(phi); phi = 0; end
if isempty(theta); theta = 0; end
%% Validity check on inputs
% -- Ensuring the user-defined material input is a cell-array
if ~iscell(lyr_mat);    lyr_mat = cellstr(lyr_mat); end
if ~iscell(lyr_conc);   lyr_conc = cellstr(lyr_conc); end
if ~iscell(lyr_ele);    lyr_ele = cellstr(lyr_ele); end
if ~iscell(lyr_cls);    lyr_cls = cellstr(lyr_cls); end
if ~iscell(lyr_thick);  lyr_thick = num2cell(lyr_thick); end
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
    for j = 1:length(lyr_mat{i})
        % - Extracting all relevant material parameters
        lyr_mpd{i}{j}	    = get_mpd_props(lyr_mat{i}{j});
        lyr_density{i}{j}   = lyr_mpd{i}{j}.DENSITY;
        lyr_atmass{i}{j}    = lyr_mpd{i}{j}.ATOM_MASS;
        lyr_atznum{i}{j}    = lyr_mpd{i}{j}.ATOM_ZNUM;
        % - Determination of the atomic concentration (n)
        pc = physics_constants();
        lyr_atstoic{i}{j}   = get_ele_ratio_from_mat(lyr_mat{i}{j}, lyr_ele{i}{j});                     % number of atoms per molecule (scales the intensity according the ratio of elements in the compound)
        lyr_n{i}{j}         = lyr_conc{i}{j} / 100 .* lyr_atstoic{i}{j} .* pc.NA.*(lyr_density{i}{j} ./ lyr_atmass{i}{j});   % number of atoms per unit volume
        % - Determination of the Sensitivity Factor (SF)
        [lyr_sf{i}{j}, lyr_sf_params{i}{j}] = calc_sf(lyr_ele{i}{j}, lyr_cls{i}{j}, lyr_mat{i}{j}, hv, theta, phi, P, formalism_xsect, formalism_imfp);
        % - Determine the value of I0
        lyr_I0{i}{j}   = lyr_n{i}{j} .* lyr_sf{i}{j};
    end
end
%% 2    :   Determination of the IMFP at all kinetic energy for each layer material
for i = 1:Nlyrs
    for j = 1:length(lyr_mat{i})
        lyr_be{i}{j}       = calc_be(lyr_ele{i}{j}, lyr_cls{i}{j}); 
        lyr_ke{i}{j}       = hv - lyr_be{i}{j};
        % - Extracting the material properties required for the TPP-2M formalism
        rho{i}{j}   = (lyr_conc{i}{j} / 100) .*lyr_mpd{i}{j}.DENSITY;
        Nv{i}{j}    = lyr_mpd{i}{j}.ELECT_VALENCY;
        M{i}{j}     = lyr_mpd{i}{j}.ATOM_MASS;
        Egap{i}{j}  = lyr_mpd{i}{j}.ELE_BGAP;
        Z{i}{j}     = lyr_mpd{i}{j}.ATOM_ZNUM;
        stoic{i}{j} = lyr_mpd{i}{j}.STOICHIOMETRY;
        if strcmpi(formalism_imfp, "TPP2M") || strcmpi(formalism_imfp, "TPP-2M") || strcmpi(formalism_imfp, "1994Tanuma") || strcmpi(formalism_imfp, "TPP2M Formalism")
            lyr_imfp{i}{j} = imfp_tpp2m(lyr_ke{i}{j}, rho{i}{j}, Nv{i}{j}, M{i}{j}, Egap{i}{j});
        elseif strcmpi(formalism_imfp, "S1") || strcmpi(formalism_imfp, "2011Seah-S1")  || strcmpi(formalism_imfp, "S1 Formalism")
            lyr_imfp{i}{j} = imfp_S1(lyr_ke{i}{j}, rho{i}{j}, M{i}{j}, Egap{i}{j}, Z{i}{j}, stoic{i}{j});
        elseif strcmpi(formalism_imfp, "S3") || strcmpi(formalism_imfp, "2012Seah-S3")
            lyr_imfp{i}{j} = eal_S3(lyr_ke{i}{j},rho{i}{j}, M{i}{j}, Z{i}{j}, stoic{i}{j});
        elseif strcmpi(formalism_imfp, "JTP") || strcmpi(formalism_imfp, "2023JTP") || strcmpi(formalism_imfp, "JTP2023")
            lyr_imfp{i}{j} = imfp_jtp(lyr_ke{i}{j}, rho{i}{j}, Nv{i}{j}, M{i}{j}, Egap{i}{j});
        end
        lyr_imfp{i}{j}     = 0.1 .* lyr_imfp{i}{j}; % Convert the IMFP into nm
    end
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
        if i == 1; for k = 1:length(lyr_mat{i}); lyr_I{i}{k} = lyr_I0{i}{k} .* (1 - exp(-lyr_thick{i} ./ (lyr_imfp{i}{k} .* cos(Theta)))); end
        % -- For the lowermost bulk layer, there is attenuation through all layers above it
        elseif i == Nlyrs
            for k = 1:length(lyr_mat{i})
                lyr_I{i}{k} = lyr_I0{i}{k};
                for j = 1:Nlyrs - 1; for l = 1:length(lyr_mat{j}); lyr_I{i}{k} = lyr_I{i}{k} .* (exp(-lyr_thick{j} ./ (lyr_imfp{j}{l} .* cos(Theta)))); end; end
            end
        % -- For intermediate layers
        else
            for k = 1:length(lyr_mat{i})
                lyr_I{i}{k} = lyr_I0{i}{k} .* (1 - exp(-lyr_thick{i} ./ (lyr_imfp{i}{k} .* cos(Theta))));
                for j = 1:i-1; for l = 1:length(lyr_mat{j}); lyr_I{i}{k} = lyr_I{i}{k} .* (exp(-lyr_thick{j} ./ (lyr_imfp{j}{l} .* cos(Theta)))); end; end
            end
        end
    end
end
% Normalising the intensities to the total sums
lyr_norm = zeros(size(lyr_I{1}{1}));
for i = 1:Nlyrs; for j = 1:length(lyr_mat{i}); lyr_norm =  lyr_norm + lyr_I{i}{j}; end; end
for i = 1:Nlyrs; for j = 1:length(lyr_mat{i}); lyr_Inorm{i}{j} = lyr_I{i}{j} ./ lyr_norm; end; end

%% 4    :   Creating a MATLAB data structure to store all of the results
% - Defining the structure
pes_model               = struct();
% - Defining the input arguments
pes_model.lyr_mat       = lyr_mat;
pes_model.lyr_conc      = lyr_conc;
pes_model.lyr_ele       = lyr_ele;
pes_model.lyr_cls       = lyr_cls;
pes_model.lyr_thick     = lyr_thick;
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
if plot_result == 1; view_nlayer_pes_model(pes_model); end
end