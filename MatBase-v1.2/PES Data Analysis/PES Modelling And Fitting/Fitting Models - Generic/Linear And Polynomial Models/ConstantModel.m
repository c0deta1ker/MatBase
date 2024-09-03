function y = ConstantModel(x, constant)
% y = ConstantModel(x, constant)
%   A model based a Constant model, with a single parameter. Note that this 
%   is ‘constant’ in the sense of having no dependence on the independent 
%   variable x, not in the sense of being non-varying. Thus, the single
%   input parameter is expected to be varied in the fitting model.
%
%   IN:
%   -   x:              N×1 (or 1×N) vector of the input domain
%   -   constant:       scalar that defines the constant value
%
%   OUT:
%   -   y:              N×1 (or 1×N) vector of the output range

%% Default parameters
if nargin < 2;          constant   = 0; end
if isempty(constant);   constant   = 0; end
%% Validity checks on the input parameters
% if isrow(x); x = x'; end           % -- Ensure x-data is a column vector
%% - 1 - Determination of the Constant Model
y = ones(size(x)).*constant;
%% Validity check on the outputs
% if isrow(y); y = y'; end   % -- Ensure y-data is a column vector
y(isnan(y)) = 0;              % -- Ensure all NaN values are zero
end