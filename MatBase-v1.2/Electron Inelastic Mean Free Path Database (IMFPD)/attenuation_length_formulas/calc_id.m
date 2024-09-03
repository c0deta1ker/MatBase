function id = calc_id(imfp, alpha, P)
% id = calc_id(imfp, alpha, P)
%   Function that calculates the information depth (ID), sometimes referred 
%   to as the sampling depth, is a useful measure of the depth in a sample from
%   which useful information is obtained for a particular material and
%   XPS configuration. An ID can be chosen as the depth from which a specified 
%   percentage, P, of the detected signal originates. Values of P such as 95% or
%   99% are often selected for this purpose. See below literature for more
%   information.
%   [1] C. J. Powell, Practical guide for inelastic mean free paths, effective attenuation lengths, mean escape depths, and information depths in x-ray photoelectron spectroscopy (2020)
%
%   IN:
%   -   imfp:       scalar or N×1 column vector of the electron IMFP values [Angstroms]
%   -   alpha:      scalar or 1×M row vector of the photoelectron take-off angle relative to the surface normal (i.e. normal emission = 0) [degrees]
%   -   P:          scalar of the percentage of the total detected signal between 0% - 100% [Default: 95%]
%
%   OUT:
%   -   id:         scalar or N×M array of the electron IL [Agstroms]

%% Default parameters
if nargin < 3; P = 95; end
if nargin < 2; alpha = 0; end
if isempty(P); P = 95; end
if isempty(alpha); alpha = 0; end
% %% -- Validity check on inputs
% if size(imfp, 2) > 1; imfp = imfp'; end         % imfp should be a column vector
% if size(alpha, 1) > 1; alpha = alpha'; end      % imfp should be a row vector
%% - 1 - Determination of the ID
id = imfp .* cos(deg2rad(alpha)) .* log(1 ./ (1 - (P/100)));
end