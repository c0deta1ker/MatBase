function y = StepModel_LHS(x, center, amplitude)
% y = StepModel_LHS(x, center, amplitude)
%   Function that evaluates a LHS Step curve profile after defining the
%   center position of the step and its amplitude. For all values > x0, 
%   the output is zero, whereas all values < x0 yield a value equal to I0.
%   Used to simulate a step-profile that is ideal & abrupt.
%
%   IN:
%   -   x:              N×1 (or 1×N) vector of the input domain
%   -   center:         scalar of the position of the step along the x-axis
%   -   amplitude:      scalar of the maximum amplitude of the step
%
%   OUT:
%   -   y:              N×1 (or 1×N) vector of the output range

%% Default parameters
if nargin < 2; center = 0; end
if nargin < 3; amplitude = 1; end
if isempty(center); center = 0; end
if isempty(amplitude); amplitude = 1; end
%% Validity check on the input parameters
% if isrow(x); x = x'; end   % -- Ensure x-data is a column vector
%% - 1 - Determination of the curve intensities
% -- Extract a Heaviside function centered at x0
H1      = heaviside(x - center);
% -- Scaling the curve profile to match the desired amplitude of I0
y    = H1;
y    = y - min(y(:));
y    = amplitude .* y ./ max(y);
%% Validity check on the outputs
% if isrow(y); y = y'; end   % -- Ensure y-data is a column vector
y(isnan(y)) = 0;              % -- Ensure all NaN values are zero
y(y<amplitude) = 0;                  % -- Ensure all values <I0 are equal to zero
end