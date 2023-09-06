function fig = view_nlayer_pes_model(pes_model, lyr_cols)
% fig = view_nlayer_pes_model(pes_model)
%   This function plots the solutions to the 'nlayer_pes_model()' function;
%   the n-layered sample stack and the corresponding photoelectron
%   contribution from each one of the layers. The user can input the colors
%   for each one of the layers using the 'lyr_cols' argument; this is
%   convenient when you want to color-match the schematic to your own
%   drawings.
%
%   IN:
%   -   pes_model:      data structure that contains all the pes model parameters and variables (from 'nlayer_pes_model()').
%   -   lyr_cols:       Mx1 cell-vector of the [R,G,B] color of each independent layer.
%
%   OUT:
%   -   fig:            MATLAB figure object with the ARPES data plotted.

%% Default parameters
def_cols = flipud(num2cell([0.25,0.25,0.25;lines(length(pes_model.lyr_mat)-1)], 2));
% -- Defining the default parameters
if nargin < 2; lyr_cols = def_cols; end
if isempty(lyr_cols); lyr_cols = def_cols; end

%% 1 - Plotting the the model solutions
% -- Range of photon energies, but a single angle
if length(pes_model.hv) > 1 && length(pes_model.theta) == 1
    fig{1} = view_nlayer_pes_model_vs_hv(pes_model, lyr_cols);
% -- Single photon energy, but a range of angles
elseif length(pes_model.hv) == 1 && length(pes_model.theta) > 1
    fig{1} = view_nlayer_pes_model_vs_theta(pes_model, lyr_cols);
% -- Range of photon energies and angles
elseif length(pes_model.hv) > 1 && length(pes_model.theta) > 1
    fig{1} = view_nlayer_pes_model_vs_hv(pes_model, lyr_cols);
    fig{2} = view_nlayer_pes_model_vs_theta(pes_model, lyr_cols);
% -- Single values only
else; fig{1} = view_nlayer_pes_model_vs_hv(pes_model, lyr_cols);
end