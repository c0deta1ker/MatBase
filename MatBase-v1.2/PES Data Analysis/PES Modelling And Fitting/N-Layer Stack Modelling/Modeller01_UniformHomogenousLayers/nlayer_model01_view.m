function fig = nlayer_model01_view(pes_model)
% fig = nlayer_model01_view(pes_model)
%   This function plots the solutions to the 'nlayer_pes_model()' function;
%   the n-layered sample stack and the corresponding photoelectron
%   contribution from each one of the layers.
%
%   IN:
%   -   pes_model:      data structure that contains all the pes model parameters and variables via 'nlayer_pes_model()'.
%
%   OUT:
%   -   fig:            output MATLAB figure.

%% Plotting the the model solutions
% -- Range of photon energies, but a single angle
if length(pes_model.hv) > 1 && length(pes_model.theta) == 1
    fig{1} = nlayer_model01_view_vs_hv(pes_model);
% -- Single photon energy, but a range of angles
elseif length(pes_model.hv) == 1 && length(pes_model.theta) > 1
    fig{1} = nlayer_model01_view_vs_theta(pes_model);
% -- Range of photon energies and angles
elseif length(pes_model.hv) > 1 && length(pes_model.theta) > 1
    fig{1} = nlayer_model01_view_vs_hv(pes_model);
    fig{2} = nlayer_model01_view_vs_theta(pes_model);
% -- Single values only
else; fig{1} = nlayer_model01_view_vs_hv(pes_model);
end