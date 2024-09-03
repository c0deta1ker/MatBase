function pes_model = nlayer_model03_run(lyr_mat, lyr_ele, lyr_cls, lyr_type, lyr_conc, lyr_thick, lyr_cdl, hv, theta, phi, P, formalism_xsect, formalism_imfp, plot_result)
% pes_model = nlayer_model03_run(lyr_mat, lyr_ele, lyr_cls, lyr_type, lyr_conc, lyr_thick, lyr_cdl, hv, theta, phi, P, formalism_xsect, formalism_imfp, plot_result)
%   This function calculates the total photoelectron intensity originating
%   from N independent layers in a sample. The sample is composed of various


%% Default parameters
% -- Default Inputs
if nargin < 14; plot_result = 0; end
if nargin < 13; formalism_imfp = "JTP"; end
if nargin < 12; formalism_xsect = "Cant2022"; end
if nargin < 11; P = 0.5;  end
if nargin < 10; phi = 0;  end
if nargin < 9; theta = 0;  end
if isempty(plot_result); plot_result = 0; end
if isempty(formalism_imfp); formalism_imfp = "JTP"; end
if isempty(formalism_xsect); formalism_xsect = "Cant2022"; end
if isempty(P); P = 0.5; end
if isempty(phi); phi = 0; end
if isempty(theta); theta = 0; end
%% Validity check on inputs
% -- Ensuring the user-defined material input is a cell-array
if ~iscell(lyr_mat);    lyr_mat = cellstr(lyr_mat); end
if ~iscell(lyr_type);   lyr_type = cellstr(lyr_type); end
if ~iscell(lyr_conc);   lyr_conc = num2cell(lyr_conc); end
if ~iscell(lyr_thick);  lyr_thick = num2cell(lyr_thick); end
if ~iscell(lyr_cdl);    lyr_cdl = num2cell(lyr_cdl); end
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
%% 1    :   Determination of the IMFP at all kinetic energy for each layer material
for i = 1:Nlyrs
    lyr_be{i}       = calc_be(lyr_ele{i}, lyr_cls{i}); 
    lyr_ke{i}       = hv - lyr_be{i};
    % - Extracting all relevant material parameters
    lyr_mpd{i}	    = get_mpd_props(lyr_mat{i});
    rho{i}          = (lyr_conc{i} / 100) .*lyr_mpd{i}.DENSITY;
    Nv{i}           = lyr_mpd{i}.ELECT_VALENCY;
    M{i}            = lyr_mpd{i}.ATOM_MASS;
    Egap{i}         = lyr_mpd{i}.ELE_BGAP;
    Z{i}            = lyr_mpd{i}.ATOM_ZNUM;
    stoic{i}        = lyr_mpd{i}.STOICHIOMETRY;
    if strcmpi(formalism_imfp, "TPP2M") || strcmpi(formalism_imfp, "TPP-2M") || strcmpi(formalism_imfp, "1994Tanuma") || strcmpi(formalism_imfp, "TPP2M Formalism")
        lyr_imfp{i} = imfp_tpp2m(lyr_ke{i}, rho{i}, Nv{i}, M{i}, Egap{i});
    elseif strcmpi(formalism_imfp, "S1") || strcmpi(formalism_imfp, "2011Seah-S1")  || strcmpi(formalism_imfp, "S1 Formalism")
        lyr_imfp{i} = imfp_S1(lyr_ke{i}, rho{i}, M{i}, Egap{i}, Z{i}, stoic{i});
    elseif strcmpi(formalism_imfp, "S3") || strcmpi(formalism_imfp, "2012Seah-S3")
        lyr_imfp{i} = eal_S3(lyr_ke{i},rho{i}, M{i}, Z{i}, stoic{i});
    elseif strcmpi(formalism_imfp, "JTP") || strcmpi(formalism_imfp, "2023JTP") || strcmpi(formalism_imfp, "JTP2023")
        lyr_imfp{i} = imfp_jtp(lyr_ke{i}, rho{i}, Nv{i}, M{i}, Egap{i});
    end
    lyr_imfp{i}     = 0.1 .* lyr_imfp{i}; % Convert the IMFP into nm
end
%% 2    :   Determination of the atomic concentration profile for each layer material
for i = 1:Nlyrs
    if i == 1;          lyr_z0{i} = 0.5 * sum(cell2mat(lyr_thick(i)));
    elseif i == Nlyrs;  lyr_z0{i} = sum(cell2mat(lyr_thick(1:i-1)));
    else;               lyr_z0{i} = sum(cell2mat(lyr_thick(1:i-1))) + 0.5 * sum(cell2mat(lyr_thick(i)));
    end
end
zpad = 100;
lyr_z = 0:0.001:sum(cell2mat(lyr_thick(1:end-1))) + zpad; lyr_z = lyr_z';
lyr_atmconc = {}; 
for i = 1:Nlyrs
    % -- For the bulk layer
    if i == Nlyrs
        if strcmpi(lyr_type{i}, "StLHS");               lyr_atmconc{i} = StepModel_LHS(lyr_z, lyr_z0{i}, lyr_conc{i});
        elseif strcmpi(lyr_type{i}, "StExLHS");         lyr_atmconc{i} = StepModel_Exp_LHS(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "StGaLHS");         lyr_atmconc{i} = StepModel_Gauss_LHS(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "ToHaExTrLHS");     lyr_atmconc{i} = StepModel_Exp_LHS_trunc(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i}, lyr_thick{i-1});
        elseif strcmpi(lyr_type{i}, "ToHaGaTrLHS");     lyr_atmconc{i} = StepModel_Gauss_LHS_trunc(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i}, lyr_thick{i-1});
        end
    % -- For each layer
    else
        if strcmpi(lyr_type{i}, "StLHS");           lyr_atmconc{i} = StepModel_LHS(lyr_z, lyr_z0{i}, lyr_conc{i});
        elseif strcmpi(lyr_type{i}, "StRHS");       lyr_atmconc{i} = StepModel_RHS(lyr_z, lyr_z0{i}, lyr_conc{i});
        elseif strcmpi(lyr_type{i}, "StExLHS");     lyr_atmconc{i} = StepModel_Exp_LHS(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "StExRHS");     lyr_atmconc{i} = StepModel_Exp_RHS(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "StGaLHS");     lyr_atmconc{i} = StepModel_Gauss_LHS(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "StGaRHS");     lyr_atmconc{i} = StepModel_Gauss_RHS(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "StExTrLHS");   lyr_atmconc{i} = StepModel_Exp_LHS_trunc(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i}, lyr_thick{i-1});
        elseif strcmpi(lyr_type{i}, "StExTrRHS");   lyr_atmconc{i} = StepModel_Exp_RHS_trunc(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i}, lyr_thick{i+1});
        elseif strcmpi(lyr_type{i}, "StGaTrLHS");   lyr_atmconc{i} = StepModel_Gauss_LHS_trunc(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i}, lyr_thick{i-1});
        elseif strcmpi(lyr_type{i}, "StGaTrRHS");   lyr_atmconc{i} = StepModel_Gauss_RHS_trunc(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_cdl{i}, lyr_thick{i+1});
        elseif strcmpi(lyr_type{i}, "ToHa");        lyr_atmconc{i} = TopHatModel(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i});
        elseif strcmpi(lyr_type{i}, "ToHaEx");      lyr_atmconc{i} = TopHatModel_Exp(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "ToHaExLHS");   lyr_atmconc{i} = TopHatModel_Exp_LHS(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "ToHaExRHS");   lyr_atmconc{i} = TopHatModel_Exp_RHS(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "ToHaGa");      lyr_atmconc{i} = TopHatModel_Gauss(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "ToHaGaLHS");   lyr_atmconc{i} = TopHatModel_Gauss_LHS(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "ToHaGaRHS");   lyr_atmconc{i} = TopHatModel_Gauss_LHS(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i}, lyr_cdl{i});
        elseif strcmpi(lyr_type{i}, "ToHaExTrLHS"); lyr_atmconc{i} = TopHatModel_Exp_LHS_trunc(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i}, lyr_cdl{i}, lyr_thick{i-1});
        elseif strcmpi(lyr_type{i}, "ToHaExTrRHS"); lyr_atmconc{i} = TopHatModel_Exp_RHS_trunc(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i}, lyr_cdl{i}, lyr_thick{i+1});
        elseif strcmpi(lyr_type{i}, "ToHaGaTrLHS"); lyr_atmconc{i} = TopHatModel_Gauss_LHS_trunc(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i}, lyr_cdl{i}, lyr_thick{i-1});
        elseif strcmpi(lyr_type{i}, "ToHaGaTrRHS"); lyr_atmconc{i} = TopHatModel_Gauss_RHS_trunc(lyr_z, lyr_z0{i}, lyr_conc{i}, lyr_thick{i}, lyr_cdl{i}, lyr_thick{i+1});
        end
    end
end
%% 3    :   Determine the maximum intensity I0 for each layer material
for i = 1:Nlyrs
    % - Extracting all relevant material parameters
    lyr_mpd{i}	    = get_mpd_props(lyr_mat{i});
    lyr_density{i}	= lyr_mpd{i}.DENSITY;
    lyr_atmass{i}   = lyr_mpd{i}.ATOM_MASS;
    lyr_atznum{i}   = lyr_mpd{i}.ATOM_ZNUM;
    % - Determination of the atomic concentration (n)
    pc = physics_constants();
    lyr_atstoic{i}  = get_ele_ratio_from_mat(lyr_mat{i}, lyr_ele{i});              % number of atoms per molecule (scales the intensity according the ratio of elements in the compound)
    lyr_n{i}        = lyr_atmconc{i} / 100 .* lyr_atstoic{i} .* pc.NA.*(lyr_density{i} ./ lyr_atmass{i});   % number of atoms per unit volume
    % - Determination of the Sensitivity Factor (SF)
    [lyr_sf{i}, lyr_sf_params{i}] = calc_sf(lyr_ele{i}, lyr_cls{i}, lyr_mat{i}, hv, theta, phi, P, formalism_xsect, formalism_imfp);
    % - Determine the value of I0
    if i == Nlyrs;          lyr_I0{i}   = lyr_n{i} .* lyr_sf{i} .* exp(-(lyr_z-lyr_z0{i}) ./ (lyr_imfp{i} .* cos(deg2rad(theta))));
    else;                   lyr_I0{i}   = lyr_n{i} .* lyr_sf{i} .* exp(-(lyr_z-lyr_z0{i}+0.5.*lyr_thick{i}) ./ (lyr_imfp{i} .* cos(deg2rad(theta))));
    end
end
if plot_result == 1
    figure(); hold on;
    cols = jet(length(theta));
    for i = 1:Nlyrs
        for j = 1:length(theta)
            plot(lyr_z, lyr_I0{i}(:,j), 'color', cols(j,:));
        end
    end
end
% -- final value of I0
for i = 1:Nlyrs
    lyr_A0{i} = trapz(lyr_z, lyr_I0{i});
end
if plot_result == 1
    figure(); hold on;
    for i = 1:Nlyrs
        plot(theta, lyr_A0{i});
    end
end

%% 4    :   Use Beer-Lambert law to determine the total intensity originating from each layer
% Defining the Beer-Lambert Law of Attenuation for each layer
Theta = deg2rad(theta);
% Attenuation through overlayers
lyr_I = {}; lyr_Inorm = {};
for i = 1:Nlyrs
    lyr_I{i} = lyr_A0{i};
    if i == 1; lyr_I{i} = lyr_I{i};
    else
        % -- For the intermediate layers, there is attenuation through all layers above it
        for j = 1:i-1; lyr_I{i} = lyr_I{i} .* exp(-lyr_thick{j} ./ (lyr_imfp{j} .* cos(Theta)));  end
    end
    lyr_I{i}(isnan(lyr_I{i})) = 0;
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
pes_model.lyr_x0        = lyr_z0;
% - Defining the properties for extracting value of I0
pes_model.lyr_I0        = lyr_I0;
pes_model.Z             = lyr_z;
pes_model.lyr_atmconc     = lyr_atmconc;
% - Defining the final intensities
pes_model.lyr_norm    	= lyr_norm;
pes_model.lyr_I         = lyr_I;
pes_model.lyr_Inorm  	= lyr_Inorm;
%% 5    :   Plotting the multilayer model stack and solutions
% if plot_result == 1; nlayer_model01_view(pes_model); end
end