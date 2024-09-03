%% MATLAB Digitisation of Cant photoionisation cross-sections
close all; clear all;
path_matbase    = what('MatBase'); path_matbase = string(path_matbase.path);
path_figs       = path_matbase + "\MatBase-v1.2\Photoionisation Energy and Fluorescence Database (PIEFD)\0_figs\";
%% 1    :   COMPARING THE BINDING ENERGY OF ALL THE FORMALISMS
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
% -- Filing through all formalisms of the BE spectra
formalism = {"Moulder1993", "Trzhaskovskaya2018", "Cant2022", "Constantinou2023"};
for i = 1:length(ATOM_SYMB)
    % -- Extracting the BE and CLs of the i'th element for all formalisms
    be = []; cls = [];
    for j = 1:length(formalism); [be{j}, cls{j}] = calc_be(ATOM_SYMB{i}, [], formalism{j}, 0); end
    % -- Plotting the BE and CLs for all formalisms
    figX = figure(); figX.Position(3) = 850; figX.Position(4) = 400;
    hold on; grid on;
    for j = 1:length(formalism)
        cols    = lines(length(be));
        for k = 1:length(be{j})
            if k == 1;  stem(be{j}(k), (5-j)/length(formalism), '-', 'linewidth', 1.0, 'marker', 'none', 'color', cols(j,:));
            else;       stem(be{j}(k), (5-j)/length(formalism), '-', 'linewidth', 1.0, 'marker', 'none', 'color', cols(j,:), 'HandleVisibility','off');
            end
            if j == 4; text(be{j}(k), (5-j)/length(formalism), sprintf('%s-%s(%.2f)', ATOM_SYMB{i}, cls{j}{k}, be{j}(k)), 'Rotation',45, 'FontWeight','bold', 'FontSize',8); end
        end
    end
    legend(formalism, 'location', 'northeast', 'FontSize', 6);
    title(sprintf("Z%i_%s_Eb_Spectrum", ATOM_ZNUM(i), ATOM_SYMB{i}), 'Interpreter','none');
    xlabel('$$ \bf Binding\ Energy\ [eV] $$', 'interpreter', 'latex');
    ylabel('$$ \bf Intensity\ [arb.] $$', 'interpreter', 'latex');
    ax = gca; ax.YScale = 'linear'; ax.XScale = 'log';
    axis([0.1, 130000, 0, 1.40]);
    print(path_figs + sprintf("Z%i_%s_Eb", ATOM_ZNUM(i), ATOM_SYMB{i}),'-dpng', '-r300');
    saveas(figX, path_figs + sprintf("Z%i_%s_Eb", ATOM_ZNUM(i), ATOM_SYMB{i}), 'fig');
    close all;
end
%% 3   :   TESTING THAT THE BE OVERLAY WORKS
view_be(["Si", "Au"]);
%% 2   :   TESTING THAT THE BE OVERLAY WORKS
figure(); hold on;
overlay_be(["Si", "Au"], 1);
