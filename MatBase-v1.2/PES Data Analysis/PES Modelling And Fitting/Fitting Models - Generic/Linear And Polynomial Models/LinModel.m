function y = LinModel(x, slope, intercept)
% y = LinModel(x, slope, intercept)
%   A model based a Linear model, with 2 input parameters; slope and
%   intercept.
%
%   IN:
%   -   x:              N×1 (or 1×N) vector of the input domain
%   -   slope:          scalar that defines the slope (dy/dx)
%   -   intercept:      scalar that defines the intercept value of the y-axis
%
%   OUT:
%   -   y:              N×1 (or 1×N) vector of the output range

%% Default parameters
if nargin < 2;          slope       = 1; end
if nargin < 3;          intercept   = 0; end
if isempty(slope);      slope       = 1; end
if isempty(intercept);  intercept   = 0; end
%% Validity checks on the input parameters
% if isrow(x); x = x'; end           % -- Ensure x-data is a column vector
%% - 1 - Determination of the Linear Model
y = slope.*x + intercept;
%% Validity check on the outputs
% if isrow(y); y = y'; end   % -- Ensure y-data is a column vector
y(isnan(y)) = 0;              % -- Ensure all NaN values are zero
end