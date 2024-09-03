function y = PowerLawModel(x, amplitude, exponent)
% ydat = PowerLawModel(x, amplitude, exponent)
%   A model based on a Power Law. The model has 2 input parameters;
%   amplitude and exponent.
%
%   IN:
%   -   x:              N×1 (or 1×N) vector of the input domain
%   -   amplitude:      scalar that defines the amplitude
%   -   exponent:       scalar that defines the exponent
%
%   OUT:
%   -   y:              N×1 (or 1×N) vector of the output range

%% Default parameters
if nargin < 2;          amplitude   = 1; end
if nargin < 3;          exponent       = 1; end
if isempty(amplitude);  amplitude   = 1; end
if isempty(exponent);   exponent       = 1; end
%% Validity checks on the input parameters
% if isrow(x); x = x'; end           % -- Ensure x-data is a column vector
%% - 1 - Determination of the Power Law Model
y = amplitude .* x .^(exponent);
%% Validity check on the outputs
% if isrow(ydat); ydat = ydat'; end   % -- Ensure y-data is a column vector
y(isnan(y)) = 0;              % -- Ensure all NaN values are zero
end