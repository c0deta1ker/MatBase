%% MATLAB Digitisation of NIST IMFP database
close all; clear all;
path_matbase    = what('MatBase'); path_matbase = string(path_matbase.path);
path_data       = path_matbase + "\MatBase-v1.2\Electron Inelastic Mean Free Path Database (IMFPD)\1999NIST_IMFP_NIST_Optical_Experiments\data_optical\";
path_data_save  = path_matbase + "\MatBase-v1.2\Electron Inelastic Mean Free Path Database (IMFPD)\1999NIST_IMFP_NIST_Optical_Experiments\";
path_fig_save   = path_matbase + "\MatBase-v1.2\Electron Inelastic Mean Free Path Database (IMFPD)\0_figs\";
%% (1)     :    Defining all of the variables
% -- Defining the elements whose photoionisation cross-sections we have
ATOM_SYMB = {...
    'H','He',...
    'Li','Be','B','C','N','O','F','Ne',...
    'Na','Mg','Al','Si','P','S','Cl','Ar',...
    'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',...
    'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',...
    'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',...
    'Fr','Ra','Ac','Th','Pa','U'};
% -- Defining the Z number for each one of these elements
ATOM_ZNUM = 1:length(ATOM_SYMB);
%% (2)     :    Iterating through all the data-files and extracting the IMFPs
% -- Extracting all the filenames for the photoionisation data
file_names = dir(path_data);
file_names = {file_names(~[file_names.isdir]).name};
% -- Filing through all the elements to isolate each element
for i = ATOM_ZNUM
    % --- Finding the index number for the element
    findx           = find(contains(file_names, "_"+string(ATOM_SYMB{i})+"_"));
    selected_names  = file_names(findx);
    eIMFP_FileNames{i} = selected_names;
end
eIMFP_Data      = {};
% -- Loading in all the data for each text file
for i = ATOM_ZNUM
    ke_dat      = [];
    imfp_dat    = [];
    for j = 1:length(eIMFP_FileNames{i})
        % -- Loading in all of the data
        data_00     = readtable(string(path_data) + string(eIMFP_FileNames{i}{j}));
        % -- Modify variable names
        data_00.Properties.VariableNames{1} = 'Energy';
        data_00.Properties.VariableNames{2} = 'IMFP';
        % -- Finding the first NaN value to crop to
        indx    = max(find(isnan(data_00.Energy)));
        ke_dat      = data_00.Energy(indx+1:end);
        imfp_dat    = data_00.IMFP(indx+1:end);
        % -- Storing the data into a cell array
        eIMFP_Data{i}{j} = [ke_dat, imfp_dat];
    end
end
% -- Finding the total number of text-files laoded for each element
eIMFP_Length    = [];
for i = ATOM_ZNUM; eIMFP_Length(i) = length(eIMFP_Data{i}); end
%% (3)     :    Interpolating the IMFPs
eIMFP_DataInterp = {};
for i = 1:length(eIMFP_Data)
    ke_dat      = [];
    imfp_dat    = [];
    for j = 1:length(eIMFP_Data{i})
        % -- Loading in the original data
        ke_dat      = eIMFP_Data{i}{j}(:,1);   % -- Electron Kinetic Energy
        imfp_dat    = eIMFP_Data{i}{j}(:,2);   % -- Inelastic Mean Free Path
        % -- Storing the data into a cell array
        ek_V    = linspace(min(ke_dat(:)), max(ke_dat(:)), 2e3);
        imfp_V  = interp1(ke_dat, imfp_dat, ek_V, 'pchip');
        % -- Appending the interpolated data into a cell array
        eIMFP_DataInterp{i}{j} = [ek_V', imfp_V'];
    end
end
%% (4)     :    Create a table for each element that sumamrises all of the IMFP data
for i = 1:length(eIMFP_Data)
    eIMFP_Table{i}  = table();
    for j = 1:length(eIMFP_Data{i})
        ek_name     = "ek_" + string(j);
        imfp_name   = "imfp_" + string(j);
        eIMFP_Table{i}.(ek_name)    = eIMFP_DataInterp{i}{j}(:,1);
        eIMFP_Table{i}.(imfp_name)  = eIMFP_DataInterp{i}{j}(:,2);
    end
end
IMFPD_NIST1999                  = struct();
IMFPD_NIST1999.ATOM_SYMB        = ATOM_SYMB;             % Atomic symbols used
IMFPD_NIST1999.ATOM_ZNUM        = ATOM_ZNUM;             % Atomic Z number for each one of these elements
IMFPD_NIST1999.FileNames        = eIMFP_FileNames;       % Cell-array of all the text file names
IMFPD_NIST1999.Length           = eIMFP_Length;          % Array where; {Atom Z number}(total number of files used)
IMFPD_NIST1999.DataRaw          = eIMFP_Data;            % Cell-array where; {Atom Z number}{File Index}(electron kinetic energy [eV], imfp [Ang])
IMFPD_NIST1999.DataInterp       = eIMFP_DataInterp;      % Cell-array where; {Atom Z number}{Core-level}(electron kinetic energy [eV], imfp [Ang])
IMFPD_NIST1999.Data             = eIMFP_Table;           % MATLAB Table where; {Atom Z number}
save(char(path_data_save + "IMFPD_NIST1999"), 'IMFPD_NIST1999', '-v7.3');
%% (5)     :    Running through all elements and plotting the imfps
close all;
for i = ATOM_ZNUM
    figX = view_imfp(ATOM_SYMB{i});  
    print(path_fig_save + sprintf("Z%i_%s_imfp", ATOM_ZNUM(i), ATOM_SYMB{i}),'-dpng', '-r500');
    saveas(figX, path_fig_save + sprintf("Z%i_%s_imfp", ATOM_ZNUM(i), ATOM_SYMB{i}), 'fig');
    close all;
end