function y = TopHatModel(x, center, amplitude, width)
% y = TopHatModel(x, center, amplitude, width)
%   Function that evaluates a Top-Hat curve profile after defining the
%   position of the centre, its amplitude and width. 
%   For all values of x outside of the Top-Hat function, the output is zero, 
%   whereas all values inside the width of the Top-Hat function yield a value
%   equal to the amplitude.
%   Used to simulate a top-hat profile that is ideal & abrupt on each side.
%
%   IN:
%   -   x:              N×1 (or 1×N) vector of the input domain
%   -   center:         scalar of the position of the step along the x-axis
%   -   amplitude:      scalar of the maximum amplitude of the step
%   -   width:          scalar of the total width of the top-hat function and should be a positive number
%
%   OUT:
%   -   y:              N×1 (or 1×N) vector of the output range

%% Default parameters
if nargin < 2; center = 0; end
if nargin < 3; amplitude = 1; end
if nargin < 4; width = 1; end
if isempty(center); center = 0; end
if isempty(amplitude); amplitude = 1; end
if isempty(width); width = 1; end
%% Validity check on the input parameters
% if isrow(x); x = x'; end   % -- Ensure x-data is a column vector
if width < 0; width = 0; end        % -- Ensure width is a positive number
%% - 1 - Determination of the curve intensities
% -- Extract a Heaviside function centered at x0 with a defined width
H1 = heaviside(x - center + 0.5*width);
H2 = heaviside(-x + center + 0.5*width);
% -- Scaling the curve profile to match the desired peak
y = H1 + H2;
y = y - min(y(:));
y = amplitude .* y ./ max(y);
%% Validity check on the outputs
% if isrow(y); y = y'; end   % -- Ensure y-data is a column vector
y(isnan(y)) = 0;              % -- Ensure all NaN values are zero
y(y<amplitude) = 0;                  % -- Ensure all values <I0 are equal to zero
end