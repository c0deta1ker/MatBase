function y = StepModel_Exp_LHS(x, center, amplitude, cdl)
% y = StepModel_Exp_LHS(x, center, amplitude, cdl)
%   Function that evaluates a LHS Step curve profile whose edge is
%   modulated via an exponential decay function. The center position of the
%   step, its amplitude and characteristic diffusion length (cdl) can 
%   be defined. For all values > x0, the output is the exponetial decay 
%   function, whereas all values < x0 yield a value equal to I0.
%   Used to simulate a step-profile that is exponentially damped.
%
%   IN:
%   -   x:              N×1 (or 1×N) vector of the input domain
%   -   center:         scalar of the position of the step along the x-axis
%   -   amplitude:      scalar of the maximum amplitude of the step
%   -   cdl:            scalar of the characteristic diffusion length (cdl), which defines the decay constant of the exponential and should be a positive number
%
%   OUT:
%   -   y:              N×1 (or 1×N) vector of the output range

%% Default parameters
if nargin < 2; center = 0; end
if nargin < 3; amplitude = 1; end
if nargin < 4; cdl = 1; end
if isempty(center); center = 0; end
if isempty(amplitude); amplitude = 1; end
if isempty(cdl); cdl = 1; end
%% Validity check on the input parameters
% if isrow(x); x = x'; end   % -- Ensure x-data is a column vector
if cdl < 0; cdl = 0; end    % -- Ensure cdl is a positive number
%% - 1 - Determination of the curve intensities
% -- Calculating the functions over the defined domain
F1      = StepModel_LHS(x, center, amplitude);
F3      = amplitude.*exp((x-center)/cdl);
% -- Isolating the function over the given domains
f1      = F1(x>center); if size(f1, 2) > 1; f1 = f1'; end
f3      = F3(x<=center); if size(f3, 2) > 1; f3 = f3'; end
% -- Concatinating the data over the domain
if isrow(f1);   y    = horzcat(f3, f1);
else;           y    = vertcat(f3, f1);
end
%% Validity check on the outputs
% if isrow(y); y = y'; end   % -- Ensure y-data is a column vector
y(isnan(y)) = 0;              % -- Ensure all NaN values are zero
end