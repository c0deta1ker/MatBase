function [fig, ele_imfp] = view_imfp(material)
% [fig, ele_imfp] = view_imfp(material)
%   This is a function that plots the electron IMFP curves of a particular material
%   defined by the user using all available calculators. This can be used to
%   quickly view all the available IMFP data for a particular element /
%   material.
%
%   IN:
%   -   material:  	string of the material whose imfp is to be shown; e.g. "Si", "SiO2", "Al2O3"...
%
%   OUT: 
%   -   fig:        figure output
%   -   ele_imfp:   data structure containing all of the data within the figure plot

%% Default parameters (Parameters for Silicon)
if nargin < 1; material = "Si";  end
if isempty(material); material = "Si"; end
%% - 1 - Extracting the material parameters from the materials database
% - Extracting the material properties
material_props = get_mpd_props(material);
% - Extracting the material properties required for the S1 formalism
rho     = material_props.DENSITY;
Nv      = material_props.ELECT_VALENCY;
M       = material_props.ATOM_MASS;
Egap    = material_props.ELE_BGAP;
Z       = material_props.ATOM_ZNUM;
%% - 2 - Extracting the eIMFP properties from each calculator
ele_imfp                 = struct();
% -- Defining the kinetic energy range
ele_imfp.Ek              = logspace(0,6,100);
% -- Unversal curve from Seah (1979)
ele_imfp.imfp_universal  = imfp_universal(ele_imfp.Ek);
% -- TPP2M predictive equations from Tanuma (1994)
ele_imfp.imfp_tpp2m      = imfp_tpp2m_mpd(ele_imfp.Ek, material);
ele_imfp.imfp_tpp2m_avg  = imfp_tpp2m_avg(ele_imfp.Ek);
% -- Optical data from NIST (1999)
[ele_imfp.imfp_opt, ele_imfp.dimfp_opt] = imfp_optical(ele_imfp.Ek, material);
% -- S1 & S2 predictive equations from Seah (2011)
ele_imfp.imfp_S1        = imfp_S1_mpd(ele_imfp.Ek, material);
ele_imfp.imfp_S2        = imfp_S2_mpd(ele_imfp.Ek, material);
% -- S3 & S4 predictive equations from Seah (2012)
ele_imfp.eal_S3         = eal_S3_mpd(ele_imfp.Ek, material);
ele_imfp.eal_S4         = eal_S4_mpd(ele_imfp.Ek, material);
% -- JTP (2023)
ele_imfp.imfp_jtp       = imfp_jtp_mpd(ele_imfp.Ek, material);
%% - 3 - Plotting a summary of the final IMFP figure
% -- Initialising the figure
fig = figure(); hold on;
grid on; grid minor;
cols = lines(10);
% - PLOTTING UNIVERSAL CURVE
loglog(ele_imfp.Ek, ele_imfp.imfp_universal, '-', 'color', 'k', 'LineWidth', 1.5);
% - PLOTTING OPTICAL/EXPERIMENTAL DATA
errorbar(ele_imfp.Ek, ele_imfp.imfp_opt, ele_imfp.dimfp_opt, ele_imfp.dimfp_opt, 'rx-', 'color', cols(1,:));
% - PLOTTING TPP2M-CURVE
loglog(ele_imfp.Ek, ele_imfp.imfp_tpp2m, ':', 'color', 'b', 'LineWidth', 1.5, 'color', cols(4,:));
% - PLOTTING TPP2M-AVERAGE-CURVE
loglog(ele_imfp.Ek, ele_imfp.imfp_tpp2m_avg, '--', 'color', 'b', 'LineWidth', 1.5, 'color', cols(4,:));
% - PLOTTING S1- & S2-CURVE
loglog(ele_imfp.Ek, ele_imfp.imfp_S1, ':', 'LineWidth', 1.5, 'color', cols(5,:));
loglog(ele_imfp.Ek, ele_imfp.imfp_S2, '--', 'LineWidth', 1.5, 'color', cols(5,:));
% - PLOTTING S3- & S4-CURVE
loglog(ele_imfp.Ek, ele_imfp.eal_S3, ':', 'LineWidth', 1.5, 'color', cols(6,:));
loglog(ele_imfp.Ek, ele_imfp.eal_S4, '--', 'LineWidth', 1.5, 'color', cols(6,:));
% - PLOTTING JTP CURVE
loglog(ele_imfp.Ek, ele_imfp.imfp_jtp, '-', 'LineWidth', 1.5, 'color', cols(7,:));
% - FORMATTING THE FIGURE
% -- Plotting the x- and y-axes
a = line([0 0], [-1e5, 1e5], 'Color', [0 0 0], 'LineWidth', 0.75, 'Linestyle', '--');
b = line([-1e5, 1e5], [0 0], 'Color', [0 0 0], 'LineWidth', 0.75, 'Linestyle', '--');
a.Annotation.LegendInformation.IconDisplayStyle = 'off';
b.Annotation.LegendInformation.IconDisplayStyle = 'off';
% -- Labelling the x- and y-axes
xlabel('$$ \bf Electron\ Kinetic\ Energy\ [eV] $$', 'interpreter', 'latex');
ylabel('$$ \bf IMFP\ [Angstrom] $$', 'interpreter', 'latex');
ax = gca;
ax.FontSize = 12;
ax.XLim = [1, 1e5]; ax.YLim = [1, 1e3];
ax.XScale = 'log'; ax.YScale = 'log';
text(0.75, 0.07, sprintf("%s (Z=%i)", material,Z), 'FontSize', 12, 'color', 'k', 'Units','normalized', 'FontWeight', 'bold');
% -- Add a legend
h       = zeros(9, 1);
h(1)    = plot(NaN,NaN,'k-', 'LineWidth', 2.0);
h(2)    = plot(NaN,NaN,'rx-', 'color',cols(1,:));
h(3)    = plot(NaN,NaN,'b:', 'LineWidth', 2.0, 'color',cols(4,:));
h(4)    = plot(NaN,NaN,'b--', 'LineWidth', 2.0, 'color',cols(4,:));
h(5)    = plot(NaN,NaN,'g:', 'LineWidth', 2.0, 'color',cols(5,:));
h(6)    = plot(NaN,NaN,'g--', 'LineWidth', 2.0, 'color',cols(5,:));
h(7)    = plot(NaN,NaN,'r:', 'LineWidth', 2.0, 'color',cols(6,:));
h(8)    = plot(NaN,NaN,'r--', 'LineWidth', 2.0, 'color',cols(6,:));
h(9)    = plot(NaN,NaN,'r-', 'LineWidth', 2.0, 'color',cols(7,:));
legend(h, 'Universal(Seah-1979)','Optical(NIST-1999)','TPP-2M(Tanuma-1994)','TPP-2M(avg)(Tanuma-1994)','S1(Seah-2011)','S2(Seah-2011)','S3[EAL](Seah-2012)','S4[EAL](Seah-2012)', 'JTP(2023)', 'fontsize', 7, 'location', 'northwest');
end