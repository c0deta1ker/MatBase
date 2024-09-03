%% MATLAB Digitisation of Cant photoionisation cross-sections
close all; clear all;
path_matbase    = what('MatBase'); path_matbase = string(path_matbase.path);
path_figs       = path_matbase + "\MatBase-v1.2\Photoionisation Cross-Section and Asymmetry Database (PIXSAD)\0_figs\";
%% 1    :   COMPARING THE CROSS-SECTIONS OF ALL THE FORMALISMS
% -- Defining the elements
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
XS_DB_Trzh2018	= load('XS_DB_Trzh2018.mat'); XS_DB_Trzh2018 = XS_DB_Trzh2018.XS_DB_Trzh2018;
% -- Filing through all elements and core-levels
for i = 1:length(ATOM_SYMB)
    DataTable = XS_DB_Trzh2018.XSECT_SIGMA{i};
    CoreLevels = DataTable.Properties.VariableNames;
    for j = 1:length(CoreLevels)
        close all;
        figX = view_xsect(ATOM_SYMB{i}, CoreLevels{j});
        print(path_figs + sprintf("Z%i_%s_%i_%s_XSECT", ATOM_ZNUM(i), ATOM_SYMB{i}, j, CoreLevels{j}),'-dpng', '-r300');
        saveas(figX, path_figs + sprintf("Z%i_%s_%i_%s_XSECT", ATOM_ZNUM(i), ATOM_SYMB{i}, j, CoreLevels{j}), 'fig');
    end
end

%% -- Comparing the known cross-sectios from literature
close all;
labels = {'C1s','Si1s','Si2p3','Ag3d5','Au4s1','Au4f7'};
hv = 1500:20:11000;
%% - 1 - SIGMA, CROSS-SECTIONS
xsect_C1s   = calc_xsect_sigma('C','1s1',hv);
xsect_Si1s  = calc_xsect_sigma('Si','1s1',hv);
xsect_Si2p3 = calc_xsect_sigma('Si','2p3',hv);
xsect_Ag3d5 = calc_xsect_sigma('Ag','3d5',hv);
xsect_Au4s1 = calc_xsect_sigma('Au','4s1',hv);
xsect_Au4f7 = calc_xsect_sigma('Au','4f7',hv);
% -- Plotting a comparison of all photon energies
barn2cm2 = 1e-24;
figure(); hold on;
plot(hv, xsect_C1s.*barn2cm2, 'g-', 'linewidth', 1.5);
plot(hv, xsect_Si1s.*barn2cm2, 'r-', 'linewidth', 1.5);
plot(hv, xsect_Si2p3.*barn2cm2, 'b-', 'linewidth', 1.5);
plot(hv, xsect_Ag3d5.*barn2cm2, 'm-','linewidth', 1.5);
plot(hv, xsect_Au4s1.*barn2cm2, 'y-','linewidth', 1.5, 'color', [0.9290 0.6940 0.1250]);
plot(hv, xsect_Au4f7.*barn2cm2, 'k-','linewidth', 1.5);
legend(labels, 'Location','best');
xlabel('$$ \bf Photon\ Energy\ [eV] $$', 'interpreter', 'latex');
ylabel('$$ \bf Cross\ Section\ [cm^2] $$', 'interpreter', 'latex');
ax = gca; ax.YScale = 'log'; ax.XScale = 'log';
axis([1000, 11000, 1e-24, 1e-18]); box on; grid on;
title('Change of \sigma of selected peaks with photon energy');
print(path_figs + "fig_Ref_XSect_Sigma",'-dpng', '-r300');

%% - 2 - BETA, ASYMMETRY PARAMETER
beta_C1s   = calc_xsect_beta('C','1s1',hv);
beta_Si1s  = calc_xsect_beta('Si','1s1',hv);
beta_Si2p3 = calc_xsect_beta('Si','2p3',hv);
beta_Ag3d5 = calc_xsect_beta('Ag','3d5',hv);
beta_Au4s1 = calc_xsect_beta('Au','4s1',hv);
beta_Au4f7 = calc_xsect_beta('Au','4f7',hv);
% -- Plotting a comparison of all photon energies
figure(); hold on;
plot(hv, beta_C1s, 'g-', 'linewidth', 1.5);
plot(hv, beta_Si1s, 'r-', 'linewidth', 1.5);
plot(hv, beta_Si2p3, 'b-', 'linewidth', 1.5);
plot(hv, beta_Ag3d5, 'm-', 'linewidth', 1.5);
plot(hv, beta_Au4s1, 'y-', 'linewidth', 1.5, 'color', [0.9290 0.6940 0.1250]);
plot(hv, beta_Au4f7, 'k-', 'linewidth', 1.5);
box on; legend(labels, 'Location','best');
xlabel('$$ \bf Photon\ Energy\ [eV] $$', 'interpreter', 'latex');
ylabel('$$ \bf \beta $$', 'interpreter', 'latex');
axis([0, 11000, 0, 2.5]); box on; grid on;
title('Change of \beta of selected peaks with photon energy');
print(path_figs + "fig_Ref_XSect_Beta",'-dpng', '-r300');

%% - 3 - GAMMA, ASYMMETRY PARAMETER
gamma_C1s   = calc_xsect_gamma('C','1s1',hv);
gamma_Si1s  = calc_xsect_gamma('Si','1s1',hv);
gamma_Si2p3 = calc_xsect_gamma('Si','2p3',hv);
gamma_Ag3d5 = calc_xsect_gamma('Ag','3d5',hv);
gamma_Au4s1 = calc_xsect_gamma('Au','4s1',hv);
gamma_Au4f7 = calc_xsect_gamma('Au','4f7',hv);
% -- Plotting a comparison of all photon energies
figure(); hold on;
plot(hv, gamma_C1s, 'g-', 'linewidth', 1.5);
plot(hv, gamma_Si1s, 'r-', 'linewidth', 1.5);
plot(hv, gamma_Si2p3, 'b-', 'linewidth', 1.5);
plot(hv, gamma_Ag3d5, 'm-', 'linewidth', 1.5);
plot(hv, gamma_Au4s1, 'y-', 'linewidth', 1.5, 'color', [0.9290 0.6940 0.1250]);
plot(hv, gamma_Au4f7, 'k-', 'linewidth', 1.5);
box on; legend(labels, 'Location','best');
xlabel('$$ \bf Photon\ Energy\ [eV] $$', 'interpreter', 'latex');
ylabel('$$ \bf \gamma $$', 'interpreter', 'latex');
axis([0, 11000, -0.5, 2.5]); box on; grid on;
title('Change of \gamma of selected peaks with photon energy');
print(path_figs + "fig_Ref_XSect_Gamma",'-dpng', '-r300');

%% - 4 - DELTA, ASYMMETRY PARAMETER
delta_C1s   = calc_xsect_delta('C','1s1',hv);
delta_Si1s  = calc_xsect_delta('Si','1s1',hv);
delta_Si2p3 = calc_xsect_delta('Si','2p3',hv);
delta_Ag3d5 = calc_xsect_delta('Ag','3d5',hv);
delta_Au4s1 = calc_xsect_delta('Au','4s1',hv);
delta_Au4f7 = calc_xsect_delta('Au','4f7',hv);
% -- Plotting a comparison of all photon energies
figure(); hold on;
plot(hv, delta_C1s, 'g-', 'linewidth', 1.5);
plot(hv, delta_Si1s, 'r-', 'linewidth', 1.5);
plot(hv, delta_Si2p3, 'b-', 'linewidth', 1.5);
plot(hv, delta_Ag3d5, 'm-', 'linewidth', 1.5);
plot(hv, delta_Au4s1, 'y-', 'linewidth', 1.5, 'color', [0.9290 0.6940 0.1250]);
plot(hv, delta_Au4f7, 'k-', 'linewidth', 1.5);
box on; legend(labels, 'Location','best');
xlabel('$$ \bf Photon\ Energy\ [eV] $$', 'interpreter', 'latex');
ylabel('$$ \bf \delta $$', 'interpreter', 'latex');
axis([1000, 10000, 0, 0.45]); box on; grid on;
title('Change of \delta of selected peaks with photon energy');
print(path_figs + "fig_Ref_XSect_Delta",'-dpng', '-r300');
