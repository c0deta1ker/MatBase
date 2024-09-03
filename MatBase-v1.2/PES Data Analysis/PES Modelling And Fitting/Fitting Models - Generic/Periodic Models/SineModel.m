function y = SineModel(x, amplitude, frequency, phaseshift)
% y = SineModel(x, amplitude, frequency, phaseshift)
%   A model based on a sinusoidal lineshape. The model has 3 input
%   parameters; amplitude, frequency and phaseshift. All are constrained to be 
%   non-negative, and the phaseshift must be between 0 and 2pi.
%
%   IN:
%   -   x:              N×1 (or 1×N) vector of the input domain
%   -   amplitude:      scalar that defines the amplitude
%   -   frequency:      scalar that defines the frequency
%   -   phaseshift:     scalar that defines the phase-shift
%
%   OUT:
%   -   y:              N×1 (or 1×N) vector of the output range

%% Default parameters
if nargin < 2;          amplitude   = 1; end
if nargin < 3;          frequency   = 1; end
if nargin < 4;          phaseshift  = 0; end
if isempty(amplitude);  amplitude   = 1; end
if isempty(frequency);  frequency   = 1; end
if isempty(phaseshift); phaseshift  = 0; end
%% Validity checks on the input parameters
% if isrow(x); x = x'; end           % -- Ensure x-data is a column vector
if amplitude < 0;   amplitude = 0; end          % -- If the amplitude is negative, pad it to zero
if frequency < 0;   frequency = 0; end          % -- If the frequency is negative, pad it to zero
if phaseshift < 0;   phaseshift = 0; end        % -- If the phaseshift is negative, pad it to zero
if phaseshift > 2*pi;   phaseshift = 2*pi; end  % -- If the phaseshift is >2pi, pad it to 2pi
%% - 1 - Determination of the Power Law Model
y = amplitude.*sin(frequency.*x + phaseshift);
%% Validity check on the outputs
% if isrow(y); y = y'; end   % -- Ensure y-data is a column vector
y(isnan(y)) = 0;              % -- Ensure all NaN values are zero
end