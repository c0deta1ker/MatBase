%% MATLAB Digitisation of Trzhaskovskaya photoionisation cross-sections
close all; clear all;
path_matbase    = what('MatBase'); path_matbase = string(path_matbase.path);
path_data           = path_matbase + "\MatBase-v1.2\Photoionisation Cross-Section and Asymmetry Database (PIXSAD)\2018Trzhaskovskaya_HAXPES_parameters\data\";
path_data_save      = path_matbase + "\MatBase-v1.2\Photoionisation Cross-Section and Asymmetry Database (PIXSAD)\2018Trzhaskovskaya_HAXPES_parameters\";
path_fig_save_sigma = path_matbase + "\MatBase-v1.2\Photoionisation Cross-Section and Asymmetry Database (PIXSAD)\2018Trzhaskovskaya_HAXPES_parameters\figs_sigma\";
path_fig_save_beta  = path_matbase + "\MatBase-v1.2\Photoionisation Cross-Section and Asymmetry Database (PIXSAD)\2018Trzhaskovskaya_HAXPES_parameters\figs_beta\";
path_fig_save_gamma = path_matbase + "\MatBase-v1.2\Photoionisation Cross-Section and Asymmetry Database (PIXSAD)\2018Trzhaskovskaya_HAXPES_parameters\figs_gamma\";
path_fig_save_delta = path_matbase + "\MatBase-v1.2\Photoionisation Cross-Section and Asymmetry Database (PIXSAD)\2018Trzhaskovskaya_HAXPES_parameters\figs_delta\";
%% (1)     :    Defining all of the variables
% -- Defining the elements whose photoionisation cross-sections we have
ATOM_SYMB = {...
    'H','He',...
    'Li','Be','B','C','N','O','F','Ne',...
    'Na','Mg','Al','Si','P','S','Cl','Ar',...
    'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',...
    'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',...
    'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',...
    'Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm'};
% -- Defining the Z number for each one of these elements
ATOM_ZNUM = 1:length(ATOM_SYMB);
% -- Defining all of the core-levels that are probed
% --- Defining in terms of orbital notation
ATOM_CL = {...
    '1s1',...
    '2s1', '2p1', '2p3',...
    '3s1', '3p1', '3p3', '3d3', '3d5',...
    '4s1', '4p1', '4p3', '4d3', '4d5', '4f5', '4f7',...
    '5s1', '5p1', '5p3', '5d3', '5d5', '5f5', '5f7',...
    '6s1', '6p1', '6p3',...
    '7s1',...
    };
%% (2)     :    Iterating through all the data-files and storing the data
% -- Extracting all the filenames for the photoionisation data
file_names = dir(path_data);
file_names = {file_names(~[file_names.isdir]).name};
% -- Filing through all the elements to isolate each element
file_names_all = strings(length(ATOM_ZNUM),1);
for i = ATOM_ZNUM
    % --- Finding the index number for the element
    findx           = find(contains(file_names, "_"+string(ATOM_SYMB{i}+".csv")));
    selected_names  = file_names(findx);
    % --- Appending the filename to a cell array
    file_names_all(i) = string(selected_names{1});
end
% -- Loading in all the data for each data file
HV = table(); HV.PhotonEnergy_eV = [1500, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]';
for i = ATOM_ZNUM
    raw_table   = readtable(string(path_data) + string(file_names_all(i)), 'VariableNamingRule','preserve');
    raw_table(:,2) = [];
    DATA_TABLE{1,i} = raw_table;
    % -- Reading in all the values from the table
    i_hv      = table2array(raw_table(1,2:end))';
    i_shell   = table2array(raw_table(2:4:end,1));
    for k = 1:length(i_shell); i_shell{k} = i_shell{k}(1:3); end
    i_sigma   = 1000.*table2array(raw_table(2:4:end,2:end));
    i_beta    = table2array(raw_table(3:4:end,2:end));
    i_gamma   = table2array(raw_table(4:4:end,2:end));
    i_delta   = table2array(raw_table(5:4:end,2:end));
    % - Assinging the values to a nice table
    % --- Initialising the tables
    XSECT_SIGMA{1,i}      = table();
    XSECT_BETA{1,i}      = table();
    XSECT_GAMMA{1,i}      = table();
    XSECT_DELTA{1,i}      = table();
    % --- Appending the photoionization parameters to tables
    for j = 1:length(i_shell)
        XSECT_SIGMA{1,i}.(i_shell{j})  = i_sigma(j,:)';
        XSECT_BETA{1,i}.(i_shell{j})   = i_beta(j,:)';
        XSECT_GAMMA{1,i}.(i_shell{j})  = i_gamma(j,:)';
        XSECT_DELTA{1,i}.(i_shell{j})  = i_delta(j,:)';
    end
end
XS_DB_Trzh2018                   = struct();
XS_DB_Trzh2018.file_names        = file_names_all;
XS_DB_Trzh2018.ATOM_SYMB         = ATOM_SYMB;
XS_DB_Trzh2018.ATOM_ZNUM         = ATOM_ZNUM;
XS_DB_Trzh2018.ATOM_CL           = ATOM_CL;
XS_DB_Trzh2018.DATA_TABLE        = DATA_TABLE;
XS_DB_Trzh2018.HV                = HV;
XS_DB_Trzh2018.XSECT_SIGMA       = XSECT_SIGMA;
XS_DB_Trzh2018.XSECT_BETA        = XSECT_BETA;
XS_DB_Trzh2018.XSECT_GAMMA       = XSECT_GAMMA;
XS_DB_Trzh2018.XSECT_DELTA       = XSECT_DELTA;
save(char(path_data_save + "XS_DB_Trzh2018"), 'XS_DB_Trzh2018', '-v7.3');
XS_DB_Trzh2018
%% (3)     :    Running through all elements and plotting the cross-sections (sigma)
close all;
hv_domain = logspace(3,5,50);
for i = ATOM_ZNUM
    xsect_sigma_trzh2018(ATOM_SYMB{i}, [], hv_domain, 1);  print(path_fig_save_sigma + sprintf("Z%i_%s_xsect_sigma", ATOM_ZNUM(i), ATOM_SYMB{i}),'-dpng', '-r500');
    close all;
end
%% (4)     :    Running through all elements and plotting the dipole & non-dipole parameters (beta, gamma and delta)
close all;
hv_domain = 1000:200:11000;
for i = ATOM_ZNUM
    xsect_beta_trzh2018(ATOM_SYMB{i}, [], hv_domain, 1);   print(path_fig_save_beta  + sprintf("Z%i_%s_xsect_beta",  ATOM_ZNUM(i), ATOM_SYMB{i}),'-dpng', '-r500');
    xsect_gamma_trzh2018(ATOM_SYMB{i}, [], hv_domain, 1);  print(path_fig_save_gamma + sprintf("Z%i_%s_xsect_gamma", ATOM_ZNUM(i), ATOM_SYMB{i}),'-dpng', '-r500');
    xsect_delta_trzh2018(ATOM_SYMB{i}, [], hv_domain, 1);  print(path_fig_save_delta + sprintf("Z%i_%s_xsect_delta", ATOM_ZNUM(i), ATOM_SYMB{i}),'-dpng', '-r500');
    close all;
end