function y = PolynomialModel(x, varargin)
% y = PolynomialModel(x, varargin)
%   A model based a generic Polynomial model. The polynomial is extracted
%   via the MATLAB function y = polyval(p,x), which evaluates the polynomial 
%   p at each point in x. The argument p is a vector of length n+1 whose 
%   elements are the coefficients (in descending powers) of an nth-degree 
%   polynomial.
%
%   IN:
%   -   x:              N×1 (or 1×N) vector of the input domain
%   -   varargin:       scalar values of the polynomial coefficients in descending powers
%
%   OUT:
%   -   y:              N×1 (or 1×N) vector of the output range

%% Default parameters
if nargin < 2; varargin = {}; end
if isempty(varargin); varargin = {}; end
%% Validity checks on the input parameters
% if isrow(x); x = x'; end           % -- Ensure x-data is a column vector
poly_coeffs = cell2mat(varargin);
%% - 1 - Determination of the Polynomial Model
y = polyval(poly_coeffs,x);
%% Validity check on the outputs
% if isrow(y); y = y'; end   % -- Ensure y-data is a column vector
y(isnan(y)) = 0;              % -- Ensure all NaN values are zero
end