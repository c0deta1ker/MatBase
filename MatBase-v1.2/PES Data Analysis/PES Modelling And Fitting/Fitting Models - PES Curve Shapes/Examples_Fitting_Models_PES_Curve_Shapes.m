close all; clear all; clc;
%% 1    :   Atomic Concentration Curves
close all;
% Defining constants
xdat = linspace(-6, 6, 1e3);
center      = 0;
amplitude   = 3.5;
width       = 6;
fwhm        = 2;
cdl         = 1;
cutoff      = 1;
%% 1.1    :   Step
ydat = PES_AtmConcCurve(xdat, "StLHS", center, amplitude);
ydat = PES_AtmConcCurve(xdat, "StRHS", center, amplitude);
%% 1.2    :   Step Erf Broadened
ydat = PES_AtmConcCurve(xdat, "StErfLHS", center, amplitude, fwhm);
ydat = PES_AtmConcCurve(xdat, "StErfRHS", center, amplitude, fwhm);
%% 1.3    :   Step Gaussian Broadened
ydat = PES_AtmConcCurve(xdat, "StGaLHS", center, amplitude, fwhm);
ydat = PES_AtmConcCurve(xdat, "StGaRHS", center, amplitude, fwhm);
%% 1.4    :   Step Gaussian Broadened & Truncated
ydat = PES_AtmConcCurve(xdat, "StGaTrLHS", center, amplitude, fwhm, cutoff);
ydat = PES_AtmConcCurve(xdat, "StGaTrRHS", center, amplitude, fwhm, cutoff);
%% 1.5    :   Step Exponential Broadened
ydat = PES_AtmConcCurve(xdat, "StExLHS", center, amplitude, cdl);
ydat = PES_AtmConcCurve(xdat, "StExRHS", center, amplitude, cdl);
%% 1.6    :   Step Exponential Broadened & Truncated
ydat = PES_AtmConcCurve(xdat, "StExTrLHS", center, amplitude, cdl, cutoff);
ydat = PES_AtmConcCurve(xdat, "StExTrRHS", center, amplitude, cdl, cutoff);
%% 1.7    :   TopHat
ydat = PES_AtmConcCurve(xdat, "ToHa", center, amplitude, width);
%% 1.8    :   TopHat Erf Broadened
ydat = PES_AtmConcCurve(xdat, "ToHaErf", center, amplitude, width, fwhm);
ydat = PES_AtmConcCurve(xdat, "ToHaErfLHS", center, amplitude, width, fwhm);
ydat = PES_AtmConcCurve(xdat, "ToHaErfRHS", center, amplitude, width, fwhm);
%% 1.9    :   TopHat Gaussian Broadened
ydat = PES_AtmConcCurve(xdat, "ToHaGa", center, amplitude, width, fwhm);
ydat = PES_AtmConcCurve(xdat, "ToHaGaLHS", center, amplitude, width, fwhm);
ydat = PES_AtmConcCurve(xdat, "ToHaGaRHS", center, amplitude, width, fwhm);
%% 1.10    :   TopHat Gaussian Broadened & Truncated
ydat = PES_AtmConcCurve(xdat, "ToHaGaTrLHS", center, amplitude, width, fwhm, cutoff);
ydat = PES_AtmConcCurve(xdat, "ToHaGaTrRHS", center, amplitude, width, fwhm, cutoff);
%% 1.11    :   TopHat Exponential Broadened
ydat = PES_AtmConcCurve(xdat, "ToHaEx", center, amplitude, width, cdl);
ydat = PES_AtmConcCurve(xdat, "ToHaExLHS", center, amplitude, width, cdl);
ydat = PES_AtmConcCurve(xdat, "ToHaExRHS", center, amplitude, width, cdl);
%% 1.12    :   TopHat Exponential Broadened & Truncated
ydat = PES_AtmConcCurve(xdat, "ToHaExTrLHS", center, amplitude, width, cdl, cutoff);
ydat = PES_AtmConcCurve(xdat, "ToHaExTrRHS", center, amplitude, width, cdl, cutoff);

%% 2    :   General XPS Lineshapes
close all;
%% 2.1    :   Generic, vanilla PES curve
% 1 - Defining the input parameters
xdat = linspace(-3, 3, 1e4);
TYPE    = "sGLA";   % type of curve to use for fitting. Default: "sGLA" ("G","L","V","DS","sGL","sGLA","pGL","pGLA")
BE      =  0.0;      % scalar of the binding energy of PE curve.
INT     =  1.0;      % scalar of the peak intensity of PE curve.
FWHM    =  0.1;     % scalar of the FWHM of the PE curve.
MR      =  0.5;     % scalar of the Mixing Ratio: 0 for pure Gaussian, 1 is for pure Lorentzian.
LSE     = -1.0;     % scalar of the binding energy of spin-orbit split PE curve.
LSI     =  0.5;     % calar of the branching ratio of spin-orbit split PE curve.
LSW     =  0.0;     % scalar of the additional lorentzian width of spin-orbit split PE curve.
ASY     =  0.;     % scalar of the PE curve asymmetry parameter (usually for metallic systems).
% 2 - Calculating the lineshape
ydat    = PES_SpecIntCurve(xdat, TYPE, BE, INT, FWHM, MR, LSE, LSI, LSW, ASY);
%% 2.2    :   PES curve convolved with a Fermi-Dirac Distribution (FDD) for near Fermi-edge fitting problems
% 1 - Defining the input parameters
xdat = linspace(-2, 2, 1e4);
TYPE    = "sGLA";   % type of curve to use for fitting. Default: "sGLA" ("G","L","V","DS","sGL","sGLA","pGL","pGLA")
BE      =  0.0;     % scalar of the binding energy of PE curve.
INT     =  1.0;     % scalar of the peak intensity of PE curve.
FWHM    =  0.1;     % scalar of the FWHM of the PE curve.
MR      =  0.5;     % scalar of the Mixing Ratio: 0 for pure Gaussian, 1 is for pure Lorentzian.
LSE     = -1.0;     % scalar of the binding energy of spin-orbit split PE curve.
LSI     =  0.5;     % calar of the branching ratio of spin-orbit split PE curve.
LSW     =  0.0;     % scalar of the additional lorentzian width of spin-orbit split PE curve.
ASY     =  0.0;     % scalar of the PE curve asymmetry parameter (usually for metallic systems).
FDEF    =  0.0;     % scalar of the FDD Fermi-Level position.
FDT     =  12.0;    % scalar of the FDD temperature.
FDW     =  0.10;    % scalar of the FDD Gaussian width after convolution.
% 2 - Calculating the lineshape
ydat     = PES_SpecIntCurve_FDD(xdat, TYPE, BE, INT, FWHM, MR, LSE, LSI, LSW, ASY, FDEF, FDT, FDW);
%% 2.3    :   PES curve integrated over a pre-defined/theoretical potential profile
% 1 - Defining the input parameters
xdat    = linspace(-3.0, 3.0, 5e3);
TYPE    = "G";   % type of curve to use for fitting. Default: "sGLA" ("pGLA", "sGL", "pGL", "G", "L", "DS")
BE      =  -1.0;    % scalar of the binding energy of PE curve.
INT     =  1.05;    % scalar of the peak intensity of PE curve.
FWHM    =  0.1;     % scalar of the FWHM of the PE curve.
MR      =  0.5;     % scalar of the Mixing Ratio: 0 for pure Gaussian, 1 is for pure Lorentzian.
LSE     = -1.0;     % scalar of the binding energy of spin-orbit split PE curve.
LSI     =  0.5;     % scalar of the branching ratio of spin-orbit split PE curve.
LSW     =  0.0;     % scalar of the additional lorentzian width of spin-orbit split PE curve.
ASY     =  0.0;     % scalar of the PE curve asymmetry parameter (usually for metallic systems).
MFP     =  5.0;     % scalar of the mean-free path of the emitted photoelectrons (either from optical, TPP-2M or fits) [nm]
ZPOT    = linspace(0, 25, 100);                           % 1xN array of the z-domain (depth) of the potential profile [eV]
EPOT    = 1./(1 + exp(-ZPOT*0.05)) - 0.01;      % 1xN array of the potential energy relative to the Fermi-level [eV]
% 2 - Calculating the lineshape
[ydat, PES_int]     = PES_SpecIntCurve_POT(xdat, TYPE, BE, INT, FWHM, MR, LSE, LSI, LSW, ASY, MFP, ZPOT, EPOT);

%% 3    :   Background Subtraction Methods
close all;
LHS         = -7;
RHS         = 2;
Win         = 0.5;
Bgr         = 0;
%% 3.1    :   Artificial XPS spectrum
xdat        = linspace(-10, 10, 1e3);
ydat        = PES_SpecIntCurve(xdat, "sGLA", -2, 1, 0.1, 0.5, -1, 0.5, 0, 0.1);
ydat        = ydat - 0.01.*xdat + 0.1;
ydat        = ydat .* (1 + 0.5*rand(1, length(ydat)));
figure(); hold on; plot(xdat, ydat);
%% 3.2    :   Polynomial Background
[roi_xdat, roi_ydat, roi_bgrnd] = PES_BgrndSubCurve(xdat, ydat, "L", LHS, RHS, Win, Bgr);
[roi_xdat, roi_ydat, roi_bgrnd] = PES_BgrndSubCurve(xdat, ydat, "P", LHS, RHS, Win, Bgr, 4);
%% 3.3    :   Shirley Background
[roi_xdat, roi_ydat, roi_bgrnd] = PES_BgrndSubCurve(xdat, ydat, "S", LHS, RHS, Win, Bgr);
[roi_xdat, roi_ydat, roi_bgrnd] = PES_BgrndSubCurve(xdat, ydat, "SO", LHS, RHS, Win, Bgr);
[roi_xdat, roi_ydat, roi_bgrnd] = PES_BgrndSubCurve(xdat, ydat, "SS", LHS, RHS, Win, Bgr);
%% 3.4    :   Tougaard Background
[roi_xdat, roi_ydat, roi_bgrnd] = PES_BgrndSubCurve(xdat, ydat, "T1", LHS, RHS, Win, Bgr);
[roi_xdat, roi_ydat, roi_bgrnd] = PES_BgrndSubCurve(xdat, ydat, "T2", LHS, RHS, Win, Bgr);
%% 3.5    :   Step Fermi-Edge Background
[roi_xdat, roi_ydat, roi_bgrnd] = PES_BgrndSubCurve(xdat, ydat, "FDDGpL", LHS, RHS, Win);
[roi_xdat, roi_ydat, roi_bgrnd] = PES_BgrndSubCurve(xdat, ydat, "FDDGsL", LHS, RHS, Win);
%% 3.6    :   No Background
[roi_xdat, roi_ydat, roi_bgrnd] = PES_BgrndSubCurve(xdat, ydat, "R", LHS, RHS, Win);

%% 4    :   Fitting an FDD Spectrum
help fdd2fit_solver
% 1 - Arbitrary 1D Fermi Edge
xdat    = linspace(-3, 3, 1e3);
ydat    = ThermalModel_FDDGsL(xdat, -0.25, 12, 0.3, -0.2, 0, 0.05) + 0.15*rand(size(xdat));
[fitStr, Ef] = fdd2fit_solver([], xdat, ydat, [], [], "FDDGsL");

%% APPENDIX
close all;
%% Appendix I - Extracting peak parameters
% 1 - Defining a curve to be used
xdat = linspace(-3, 3, 1e3);
TYPE    = "sGLA";   % type of curve to use for fitting. Default: "sGLA" ("pGLA", "sGL", "pGL", "G", "L", "DS")
BE      = 0.0;      % scalar of the binding energy of PE curve.
INT     = 1.0;      % scalar of the peak intensity of PE curve.
FWHM    =  0.1;     % scalar of the FWHM of the PE curve.
MR      =  0.5;     % scalar of the Mixing Ratio: 0 for pure Gaussian, 1 is for pure Lorentzian.
LSE     = -1.0;     % scalar of the binding energy of spin-orbit split PE curve.
LSI     =  0.5;     % calar of the branching ratio of spin-orbit split PE curve.
LSW     =  0.0;     % scalar of the additional lorentzian width of spin-orbit split PE curve.
ASY     =  0.0;     % scalar of the PE curve asymmetry parameter (usually for metallic systems).
% 2 - Calculating the lineshape
ydat    = PES_SpecIntCurve(xdat, TYPE, BE, INT, FWHM, MR, LSE, LSI, LSW, ASY);
% 3 - Plotting the lineshapes & verifying the FWHM values
figure(); plot(xdat, ydat, 'b-', 'linewidth', 2);
title('PESCurve()'); 
xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
ylabel('$$ \bf  Intensity $$', 'Interpreter', 'latex');
%% Appendix II - Extracting peak positions & heights
close all;
help find_peak_loc; 
Win = 0.02 .* [-1, 1];
[xval_R, yval_R]            = find_peak_loc(xdat, ydat, Win, "R", 1);
[xval_S, yval_S]            = find_peak_loc(xdat, ydat, Win, "S", 1);
[xval_G1, yval_G1]          = find_peak_loc(xdat, ydat, Win, "G1", 1);
[xval_D1, yval_D1]          = find_peak_loc(xdat, ydat, Win, "D1", 1);
[xval_D2, yval_D2]          = find_peak_loc(xdat, ydat, Win, "D2", 1);
[xval_G, yval_G]            = find_peak_loc(xdat, ydat, Win, "G", 1);
[xval_L, yval_L]            = find_peak_loc(xdat, ydat, Win, "L", 1);
[xval_sGL, yval_sGL]        = find_peak_loc(xdat, ydat, Win, "sGL", 1);
[xval_pGL, yval_pGL]        = find_peak_loc(xdat, ydat, Win, "pGL", 1);
[xval_sGLA, yval_sGLA]      = find_peak_loc(xdat, ydat, Win, "sGLA", 1);
[xval_pGLA, yval_pGLA]      = find_peak_loc(xdat, ydat, Win, "pGLA", 1);
% -- Plotting a comparison of all peak finder methods
fig = figure(); fig.Position(3) = 450; fig.Position(4) = 450; hold on; 
xline(0, 'HandleVisibility','off'); yline(1, 'HandleVisibility','off'); 
h = plot(0, 1, 'kx', 'markersize', 10, 'HandleVisibility','off');
cols = lines(20);
plot(xval_R, yval_R, 'kx', 'color', cols(1,:));
plot(xval_S, yval_S, 'rx', 'color', cols(2,:));
plot(xval_G1, yval_G1, 'rx', 'color', cols(3,:));
plot(xval_D1, yval_D1, 'rx', 'color', cols(4,:));
plot(xval_D2, yval_D2, 'rx', 'color', cols(5,:));
plot(xval_G, yval_G, 'kx', 'color', cols(6,:));
plot(xval_L, yval_L, 'rx', 'color', cols(7,:));
plot(xval_sGL, yval_sGL, 'rx', 'color', cols(8,:));
plot(xval_pGL, yval_pGL, 'rx', 'color', cols(9,:));
plot(xval_sGLA, yval_sGLA, 'rx', 'color', cols(10,:));
plot(xval_pGLA, yval_pGLA, 'rx', 'color', cols(11,:));
lgnd = {"Raw", "Spline", "Gaco1", "dydx", "d2ydx2", "G", "L", "sGL", "pGL", "sGLA", "pGLA"};
legend(lgnd, 'location', 'best','fontsize', 7);
axis([-0.015, 0.015, 0.985, 1.015]);
xlabel(' X ', 'fontweight', 'bold');
ylabel(' Y ', 'fontweight', 'bold');
title(sprintf("find_peak_loc()"), 'interpreter', 'none'); 
%% Appendix III - Extracting peak FWHM
close all;
help find_peak_fwhm; 
Win = 0.08 .* [-1, 1];
[fwhm_R, fwhmLocs_R]            = find_peak_fwhm(xdat, ydat, Win, "R", 1);
[fwhm_S, fwhmLocs_S]            = find_peak_fwhm(xdat, ydat, Win, "S", 1);
[fwhm_G1, fwhmLocs_G1]          = find_peak_fwhm(xdat, ydat, Win, "G1", 1);
[fwhm_G, fwhmLocs_G]            = find_peak_fwhm(xdat, ydat, Win, "G", 1);
[fwhm_L, fwhmLocs_L]            = find_peak_fwhm(xdat, ydat, Win, "L", 1);
[fwhm_sGL, fwhmLocs_sGL]        = find_peak_fwhm(xdat, ydat, Win, "sGL", 1);
[fwhm_pGL, fwhmLocs_pGL]        = find_peak_fwhm(xdat, ydat, Win, "pGL", 1);
[fwhm_sGLA, fwhmLocs_sGLA]      = find_peak_fwhm(xdat, ydat, Win, "sGLA", 1);
[fwhm_pGLA, fwhmLocs_pGLA]      = find_peak_fwhm(xdat, ydat, Win, "pGLA", 1);
% -- Plotting a comparison of all peak finder methods
fig = figure(); fig.Position(3) = 450; fig.Position(4) = 450; hold on; 
xline(0.1, 'HandleVisibility','off'); yline(0.1, 'HandleVisibility','off'); 
h = plot(FWHM, FWHM, 'kx', 'markersize', 10, 'HandleVisibility','off');
cols = lines(20);
plot(fwhm_R, fwhm_R, 'kx', 'color', cols(1,:));
plot(fwhm_S, fwhm_S, 'rx', 'color', cols(2,:));
plot(fwhm_G1, fwhm_G1, 'rx', 'color', cols(3,:));
plot(fwhm_G, fwhm_G, 'kx', 'color', cols(4,:));
plot(fwhm_L, fwhm_L, 'rx', 'color', cols(5,:));
plot(fwhm_sGL, fwhm_sGL, 'rx', 'color', cols(6,:));
plot(fwhm_pGL, fwhm_pGL, 'rx', 'color', cols(7,:));
plot(fwhm_sGLA, fwhm_sGLA, 'rx', 'color', cols(8,:));
plot(fwhm_pGLA, fwhm_pGLA, 'rx', 'color', cols(9,:));
lgnd = {"Raw", "Spline", "Gaco1",  "G", "L", "sGL", "pGL", "sGLA", "pGLA"};
legend(lgnd, 'location', 'best','fontsize', 7);
axis([0.08, 0.12, 0.08, 0.12]);
xlabel(' X ', 'fontweight', 'bold');
ylabel(' Y ', 'fontweight', 'bold');
title(sprintf("find_peak_fwhm()"), 'interpreter', 'none'); 
