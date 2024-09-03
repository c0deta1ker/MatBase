%% MATLAB Digitisation of Cant photoionisation cross-sections
close all; clear all;
path_matbase    = what('MatBase'); path_matbase = string(path_matbase.path);
path_coeffs                 = path_matbase + "\MatBase-v1.2\Photoionisation Cross-Section and Asymmetry Database (PIXSAD)\2022Cant_AMRSF\coeff_tables\";
path_data_save              = path_matbase + "\MatBase-v1.2\Photoionisation Cross-Section and Asymmetry Database (PIXSAD)\2022Cant_AMRSF\";
path_data_save_sigma        = path_matbase + "\MatBase-v1.2\Photoionisation Cross-Section and Asymmetry Database (PIXSAD)\2022Cant_AMRSF\data_sigma\";
path_data_save_beta         = path_matbase + "\MatBase-v1.2\Photoionisation Cross-Section and Asymmetry Database (PIXSAD)\2022Cant_AMRSF\data_beta\";
path_data_save_gamma        = path_matbase + "\MatBase-v1.2\Photoionisation Cross-Section and Asymmetry Database (PIXSAD)\2022Cant_AMRSF\data_gamma\";
path_data_save_delta        = path_matbase + "\MatBase-v1.2\Photoionisation Cross-Section and Asymmetry Database (PIXSAD)\2022Cant_AMRSF\data_delta\";
path_fig_save_sigma         = path_matbase + "\MatBase-v1.2\Photoionisation Cross-Section and Asymmetry Database (PIXSAD)\2022Cant_AMRSF\figs_sigma\";
path_fig_save_beta          = path_matbase + "\MatBase-v1.2\Photoionisation Cross-Section and Asymmetry Database (PIXSAD)\2022Cant_AMRSF\figs_beta\";
path_fig_save_gamma         = path_matbase + "\MatBase-v1.2\Photoionisation Cross-Section and Asymmetry Database (PIXSAD)\2022Cant_AMRSF\figs_gamma\";
path_fig_save_delta         = path_matbase + "\MatBase-v1.2\Photoionisation Cross-Section and Asymmetry Database (PIXSAD)\2022Cant_AMRSF\figs_delta\";
%% (1)     :    Defining all of the variables
% -- Defining the elements whose photoionisation cross-sections we have
ATOM_SYMB = {...
    'H','He',...
    'Li','Be','B','C','N','O','F','Ne',...
    'Na','Mg','Al','Si','P','S','Cl','Ar',...
    'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',...
    'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',...
    'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',...
    'Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf'};
% -- Defining the Z number for each one of these elements
ATOM_ZNUM = 1:length(ATOM_SYMB);
% -- Defining all of the core-levels that are probed
% --- Defining in terms of orbital notation (for sigma)
ATOM_CL_01 = {...
    '1s1', '2s1', '3s1', '4s1', '5s1',...
    '2p1', '2p3','3p1', '3p3','4p1', '4p3','5p1', '5p3'...
    '3d3', '3d5','4d3', '4d5',...
    '4f5', '4f7',...
    };
% --- Defining in terms of orbital notation (for beta, gamma, delta)
ATOM_CL_02 = {...
    '1s1', '2s1', '3s1', '4s1', '5s1',...
    '2p1', '2p3','3p1', '3p3','4p1', '4p3','5p1', '5p3'...
    '3d3', '3d5','4d3', '4d5','5d3', '5d5',...
    '4f5', '4f7',...
    };
% --- Defining the photon energy domain
hv_dom = 1e3:50:1.5e4; hv_dom = hv_dom';
HV = table(); HV.('PhotonEnergy_eV') = hv_dom;
%% (2)     :    Calculating the photoionization cross-sections
XSECT_SIGMA = {};XSECT_BETA = {};XSECT_GAMMA = {};XSECT_DELTA = {};
for Z = ATOM_ZNUM
    XSECT_SIGMA{Z}  = table();  XSECT_SIGMA{Z}.('PhotonEnergy_eV')  = hv_dom;
    XSECT_BETA{Z}   = table();  XSECT_BETA{Z}.('PhotonEnergy_eV')   = hv_dom;
    XSECT_GAMMA{Z}  = table();  XSECT_GAMMA{Z}.('PhotonEnergy_eV')  = hv_dom;
    XSECT_DELTA{Z}  = table();  XSECT_DELTA{Z}.('PhotonEnergy_eV')  = hv_dom;
    % -- Evaluating sigma
    for j = 1:length(ATOM_CL_01)
        XSECT_SIGMA{Z}.(ATOM_CL_01{j}) = xsect_sigma_cant2022(ATOM_SYMB{Z}, ATOM_CL_01{j}, hv_dom);
    end
    writetable(XSECT_SIGMA{Z}, path_data_save_sigma + sprintf("Z_%i_%s.csv", Z, ATOM_SYMB{Z}));
    XSECT_SIGMA{Z}(:,1) = [];   % delete the photon energy from the data structure
    % -- Evaluating beta, gamma & delta
    for j = 1:length(ATOM_CL_02)
        XSECT_BETA{Z}.(ATOM_CL_02{j})   = xsect_beta_cant2022(ATOM_SYMB{Z}, ATOM_CL_02{j}, hv_dom);
        XSECT_GAMMA{Z}.(ATOM_CL_02{j})  = xsect_gamma_cant2022(ATOM_SYMB{Z}, ATOM_CL_02{j}, hv_dom);
        XSECT_DELTA{Z}.(ATOM_CL_02{j})  = xsect_delta_cant2022(ATOM_SYMB{Z}, ATOM_CL_02{j}, hv_dom);
    end
    % ---- Saving a .txt file for the dipole parameter, beta
    writetable(XSECT_BETA{Z}, path_data_save_beta + sprintf("Z_%i_%s.csv", Z, ATOM_SYMB{Z}));
    XSECT_BETA{Z}(:,1) = [];   % delete the photon energy from the data structure
    % ---- Saving a .txt file for the non-dipole parameter, gamma
    writetable(XSECT_GAMMA{Z}, path_data_save_gamma + sprintf("Z_%i_%s.csv", Z, ATOM_SYMB{Z}));
    XSECT_GAMMA{Z}(:,1) = [];   % delete the photon energy from the data structure
    % ---- Saving a .txt file for the non-dipole parameter, delta
    writetable(XSECT_DELTA{Z}, path_data_save_delta + sprintf("Z_%i_%s.csv", Z, ATOM_SYMB{Z}));
    XSECT_DELTA{Z}(:,1) = [];   % delete the photon energy from the data structure
end
%% (3)     :    Appending the data into a MATLAB structure file
XS_DB_Cant2022                  = struct();
XS_DB_Cant2022.ATOM_SYMB        = ATOM_SYMB;
XS_DB_Cant2022.ATOM_ZNUM        = ATOM_ZNUM;
XS_DB_Cant2022.HV               = HV;
XS_DB_Cant2022.XSECT_SIGMA      = XSECT_SIGMA;
XS_DB_Cant2022.XSECT_BETA       = XSECT_BETA;
XS_DB_Cant2022.XSECT_GAMMA      = XSECT_GAMMA;
XS_DB_Cant2022.XSECT_DELTA      = XSECT_DELTA;
save(char(path_data_save + "XS_DB_Cant2022"), 'XS_DB_Cant2022', '-v7.3');
XS_DB_Cant2022
%% (7)     :    Running through all elements and plotting the cross-sections (sigma)
close all;
for i = ATOM_ZNUM
    xsect_sigma_cant2022(ATOM_SYMB{i}, [], [], 1);      print(path_fig_save_sigma + sprintf("Z%i_%s_xsect_sigma", ATOM_ZNUM(i), ATOM_SYMB{i}),'-dpng', '-r500');
    xsect_beta_cant2022(ATOM_SYMB{i}, [], [], 1);       print(path_fig_save_beta + sprintf("Z%i_%s_xsect_beta", ATOM_ZNUM(i), ATOM_SYMB{i}),'-dpng', '-r500');
    xsect_gamma_cant2022(ATOM_SYMB{i}, [], [], 1);      print(path_fig_save_gamma + sprintf("Z%i_%s_xsect_gamma", ATOM_ZNUM(i), ATOM_SYMB{i}),'-dpng', '-r500');
    xsect_delta_cant2022(ATOM_SYMB{i}, [], [], 1);      print(path_fig_save_delta + sprintf("Z%i_%s_xsect_delta", ATOM_ZNUM(i), ATOM_SYMB{i}),'-dpng', '-r500');
    close all;
end