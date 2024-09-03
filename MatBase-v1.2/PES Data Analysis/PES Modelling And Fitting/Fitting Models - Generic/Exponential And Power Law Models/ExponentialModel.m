function y = ExponentialModel(x, amplitude, decay)
% y = ExponentialModel(x, amplitude, decay)
%   A model based on an exponential decay function. The model has 2 input
%   parameters; amplitude and decay.
%
%   IN:
%   -   x:              N×1 (or 1×N) vector of the input domain
%   -   amplitude:      scalar that defines the amplitude of the exponential decay
%   -   decay:          scalar that defines the decay (or rate) constant
%
%   OUT:
%   -   y:              N×1 (or 1×N) vector of the output range

%% Default parameters
if nargin < 2;          amplitude   = 1; end
if nargin < 3;          decay       = 1; end
if isempty(amplitude);  amplitude   = 1; end
if isempty(decay);      decay       = 1; end
%% Validity checks on the input parameters
% if isrow(x); x = x'; end           % -- Ensure x-data is a column vector
%% - 1 - Determination of the Exponential Model
y = amplitude .* exp(-x / decay);
%% Validity check on the outputs
% if isrow(ydat); ydat = ydat'; end   % -- Ensure y-data is a column vector
y(isnan(y)) = 0;              % -- Ensure all NaN values are zero
end