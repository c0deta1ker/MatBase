close all; clear all;
path_matbase    = what('MatBase'); path_matbase = string(path_matbase.path);
%% 1    :   MATLAB Digitisation of Materials Properties Database (FERTIG!)
run import2matlab_MPD_PCC.m; close all; clear all;
%% 2    :   MATLAB Digitisation of Inelastic Mean Free Path (IMFP) calculators (FERTIG!)
run import2matlab_IMFP_DB_NIST1999.m; close all; clear all;
%% 3    :   MATLAB Digitisation of Photoionisation Binding Energies (FERTIG!)
run import2matlab_BE_DB_1993Moulder.m; close all; clear all;
run import2matlab_BE_DB_2018Trzh.m; close all; clear all;
run import2matlab_BE_DB_2022Cant.m; close all; clear all;
run import2matlab_BE_DB_2023Constantinou.m; close all; clear all;
run run_piefd_meta_analysis.m; close all; clear all;
%% 4    :   MATLAB Digitisation of Photoionisation Cross-Sections & Asymmetry Parameters (FERTIG!)
run import2matlab_XS_DB_Scof1973.m; close all; clear all;
run import2matlab_XS_DB_YehLind1985.m; close all; clear all;
run import2matlab_XS_DB_Trzha2018.m; close all; clear all;
run import2matlab_XS_DB_Cant2022.m; close all; clear all;
run run_pixad_meta_analysis.m; close all; clear all;
%% 5    :   MATLAB Digitisation of Sensitivity Factors (FERTIG!)
run run_sf_vs_hv_meta_analysis.m; close all; clear all;
run run_sf_vs_theta_meta_analysis.m; close all; clear all;
%% 6    :   MATLAB Digitisation of Generic & PES Curve Shapes (FERTIG!)
run Examples_Fitting_Models_Generic.m; close all; clear all;
