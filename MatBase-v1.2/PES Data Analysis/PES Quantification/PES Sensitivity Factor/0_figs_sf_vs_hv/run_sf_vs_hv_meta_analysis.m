close all; clear all;
%% MATLAB Digitisation of SFs
close all; clear all;
path_matbase    = what('MatBase'); path_matbase = string(path_matbase.path);
path_figs       = path_matbase + "\MatBase-v1.2\PES Data Analysis\PES Quantification\PES Sensitivity Factor\0_figs_sf_vs_hv\";
%% 1    :   Calculating the SFs From My Library
% -- Defining the elements
ATOM_SYMB = {...
    'H','He',...
    'Li','Be','B','C','N','O','F','Ne',...
    'Na','Mg','Al','Si','P','S','Cl','Ar',...
    'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',...
    'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',...
    'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',...
    'Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf'};
% -- Defining the core-levels to be probed
ATOM_ZNUM = 1:length(ATOM_SYMB);
XS_DB_Trzh2018	= load('XS_DB_Trzh2018.mat'); XS_DB_Trzh2018 = XS_DB_Trzh2018.XS_DB_Trzh2018;
% -- Assuming normal emission geometry
theta    = 0;
phi      = 0;
P        = 0.5;
hv_SX    = 250:10:1500;
hv_HX    = 500:10:100000;
% -- Filing through all elements and core-levels
for i = 1:length(ATOM_SYMB)
    close all;
    DataTable   = XS_DB_Trzh2018.XSECT_SIGMA{i};
    CoreLevels  = DataTable.Properties.VariableNames;
    SF_SX = {}; SF_HX = {};
    for j = 1:length(CoreLevels)
        SF_SX{j}   = calc_sf(ATOM_SYMB{i}, CoreLevels{j}, ATOM_SYMB{i}, hv_SX, theta, phi, P, 'yl');
        SF_HX{j}   = calc_sf(ATOM_SYMB{i}, CoreLevels{j}, ATOM_SYMB{i}, hv_HX, theta, phi, P, 'cant2022');
    end
    % -- Plotting the RSF
    fig = figure(); hold on; grid on; grid minor;
    fig.Position(1) = 100; fig.Position(2) = 100;
    fig.Position(3) = 600;fig.Position(4) = 500;
    cols = lines(length(CoreLevels));
    for j = 1:length(CoreLevels); plot(hv_HX, SF_HX{j}, 'k-', 'LineWidth',1.5, 'color', cols(j,:)); end
    % for j = 1:length(CoreLevels); plot(hv_SX, SF_SX{j}, 'k:', 'LineWidth',1.5, 'color', cols(j,:)); end
    title(sprintf("Z = %i, %s", i, ATOM_SYMB{i}),'FontWeight','bold','interpreter', 'none', 'fontsize', 16);
    legend(CoreLevels, 'fontsize', 7, 'location', 'eastoutside');
    xlabel(' Photon Energy [eV] ', 'FontWeight','bold');
    ylabel(' Sensitivity Factor (SF) ', 'FontWeight','bold');
    ax = gca; ax.YScale = 'log'; ax.XScale = 'log';
    axis([1e2, 1e5, 1e-3, 1e9]);
    % -- Adding annotation with the experimental variables
    txt_str = sprintf("theta = %.2f deg. \nphi = %.2f deg. \nP = %.2f", theta, phi, P);
    hAnnotation = annotation('textbox',[0.15, 0.15, 0.1, 0.1],'string', txt_str, 'FontSize', 8, 'BackgroundColor','w');
    print(path_figs + sprintf("Z%i_%s_SF_vs_HV", ATOM_ZNUM(i), ATOM_SYMB{i}),'-dpng', '-r300');
    saveas(fig, path_figs + sprintf("Z%i_%s_SF_vs_HV", ATOM_ZNUM(i), ATOM_SYMB{i}), 'fig');
end