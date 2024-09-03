function y = PES_AtmConcCurve(x, Type, varargin)
% y = PES_AtmConcCurve(x, Type, varargin)
%   Function that evaluates a generic atomic concentration curve profile by
%   defining the type, along with the necessary arguments. See below for
%   all available curve types and what input arguments are required for
%   each case. Options for a step-like and tophat-like distribution are
%   available. Additionally, the Error-, Gaussian- and
%   Exponential-Broadened versions of the LHS/RHS or both sides of the
%   distribution are available.
%
%   IN:
%   -   x:          N×1 (or 1×N) vector of the input domain along the x-axis (typically the depth axis)
%   -   Type:       string of the type of atomic concentration curve profile from the following list:
%                       - Step:                                         "StLHS", "StRHS"
%                       - Step Erf Broadened:                           "StErfLHS", "StErfRHS"
%                       - Step Gaussian Broadened:                      "StGaLHS", "StGaRHS"
%                       - Step Gaussian Broadened & Truncated:          "StGaTrLHS", "StGaTrRHS"
%                       - Step Exponential Broadened:                   "StExLHS", "StExRHS"
%                       - Step Exponential Broadened & Truncated:       "StExTrLHS", "StExTrRHS"
%                       - TopHat:                                       "ToHa"
%                       - TopHat Erf Broadened:                         "ToHaErf", "ToHaErfLHS", "ToHaErfRHS"
%                       - TopHat Gaussian Broadened:                    "ToHaGa", "ToHaGaLHS", "ToHaGaRHS"
%                       - TopHat Gaussian Broadened & Truncated:        "ToHaGaTrLHS", "ToHaGaTrRHS"
%                       - TopHat Exponential Broadened:                 "ToHaEx", "ToHaExLHS", "ToHaExRHS"
%                       - TopHat Exponential Broadened & Truncated:     "ToHaExTrLHS", "ToHaExTrRHS"
%   -   varargin:   arguments to be inserted depending on the type of atomic concentration curve profile defined.
%                       - Step:                                         [center, amplitude]
%                       - Step Erf Broadened:                           [center, amplitude, fwhm]
%                       - Step Gaussian Broadened:                      [center, amplitude, fwhm]
%                       - Step Gaussian Broadened & Truncated:          [center, amplitude, fwhm, cutoff]
%                       - Step Exponential Broadened:                   [center, amplitude, cdl]
%                       - Step Exponential Broadened & Truncated:       [center, amplitude, cdl, cutoff]
%                       - TopHat:                                       [center, amplitude, width]
%                       - TopHat Erf Broadened:                         [center, amplitude, width, fwhm]
%                       - TopHat Gaussian Broadened:                    [center, amplitude, width, fwhm]
%                       - TopHat Gaussian Broadened & Truncated:        [center, amplitude, width, fwhm, cutoff]
%                       - TopHat Exponential Broadened:                 [center, amplitude, width, cdl]
%                       - TopHat Exponential Broadened & Truncated:     [center, amplitude, width, cdl, cutoff]
%
%   OUT:
%   -   y:     	    N×1 (or 1×N) vector of the output atomic concentration curve profile

%% Default parameters
if nargin < 3; varargin = {}; end
if isempty(varargin); varargin = {}; end
%% - 1 - Determination of the PES Atomic Concentration Curve
Type = char(lower(Type));
switch Type
    % - Step Functions
    case lower({'StepModel_LHS','StLHS','SLHS'});                           y = StepModel_LHS(x, varargin{:}); label="StLHS";
    case lower({'StepModel_RHS','StRHS','SRHS'});                           y = StepModel_RHS(x, varargin{:}); label="StRHS";
    case lower({'StepErfModel_LHS','StErfLHS','SERFLHS'});                  y = StepModel_Erf_LHS(x, varargin{:}); label="StErfLHS";
    case lower({'StepErfModel_RHS','StErfRHS','SERFRHS'});                  y = StepModel_Erf_RHS(x, varargin{:}); label="StErfRHS";
    case lower({'StepModel_Exp_LHS','StExLHS','SELHS'});                    y = StepModel_Exp_LHS(x, varargin{:}); label="StExLHS";
    case lower({'StepModel_Exp_RHS','StExRHS','SERHS'});                    y = StepModel_Exp_RHS(x, varargin{:}); label="StExRHS";
    case lower({'StepModel_Exp_LHS_trunc','StExTrLHS','SETLHS'});           y = StepModel_Exp_LHS_trunc(x, varargin{:}); label="StExTrLHS";
    case lower({'StepModel_Exp_RHS_trunc','StExTrRHS','SETRHS'});           y = StepModel_Exp_RHS_trunc(x, varargin{:}); label="StExTrRHS";
    case lower({'StepModel_Gauss_LHS','StGaLHS','SGLHS'});                  y = StepModel_Gauss_LHS(x, varargin{:}); label="StGaLHS";
    case lower({'StepModel_Gauss_RHS','StGaRHS','SGLHS'});                  y = StepModel_Gauss_RHS(x, varargin{:}); label="StGaRHS";
    case lower({'StepModel_Gauss_LHS_trunc','StGaTrLHS','SGTLHS'});         y = StepModel_Gauss_LHS_trunc(x, varargin{:}); label="StGaTrLHS";
    case lower({'StepModel_Gauss_RHS_trunc','StGaTrRHS','SGTRHS'});         y = StepModel_Gauss_RHS_trunc(x, varargin{:}); label="StGaTrRHS";
    % - Top-Hat Functions
    case lower({'TopHatModel','ToHa','TH'});                                    y = TopHatModel(x, varargin{:}); label="ToHa";
    case lower({'TopHatModel_Erf','ToHaErf','THERF'});                          y = TopHatModel_Erf(x, varargin{:}); label="ToHaErf";
    case lower({'TopHatModel_Erf_LHS','ToHaErfLHS','THERFLHS'});                y = TopHatModel_Erf_LHS(x, varargin{:}); label="ToHaErfLHS";
    case lower({'TopHatModel_Erf_RHS','ToHaErfRHS','THERFRHS'});                y = TopHatModel_Erf_RHS(x, varargin{:}); label="ToHaErfRHS";
    case lower({'TopHatModel_Exp','ToHaEx','THE'});                             y = TopHatModel_Exp(x, varargin{:}); label="ToHaEx";
    case lower({'TopHatModel_Exp_LHS','ToHaExLHS','THELHS'});                   y = TopHatModel_Exp_LHS(x, varargin{:}); label="ToHaExLHS";
    case lower({'TopHatModel_Exp_RHS','ToHaExRHS','THERHS'});                   y = TopHatModel_Exp_RHS(x, varargin{:}); label="ToHaExRHS";
    case lower({'TopHatModel_Exp_LHS_trunc','ToHaExTrLHS','THETLHS'});          y = TopHatModel_Exp_LHS_trunc(x, varargin{:}); label="ToHaExTrLHS";
    case lower({'TopHatModel_Exp_RHS_trunc','ToHaExTrRHS','THETRHS'});          y = TopHatModel_Exp_RHS_trunc(x, varargin{:}); label="ToHaExTrRHS";
    case lower({'TopHatModel_Gauss','ToHaGa','THG'});                           y = TopHatModel_Gauss(x, varargin{:}); label="ToHaGa";
    case lower({'TopHatModel_Gauss_LHS','ToHaGaLHS','THGLHS'});                 y = TopHatModel_Gauss_LHS(x, varargin{:}); label="ToHaGaLHS";
    case lower({'TopHatModel_Gauss_RHS','ToHaGaRHS','THGRHS'});                 y = TopHatModel_Gauss_RHS(x, varargin{:}); label="ToHaGaRHS";
    case lower({'TopHatModel_Gauss_LHS_trunc','ToHaGaTrLHS','THGTLHS'});        y = TopHatModel_Gauss_LHS_trunc(x, varargin{:}); label="ToHaGaTrLHS";
    case lower({'TopHatModel_Gauss_RHS_trunc','ToHaGaTrRHS','THGTRHS'});        y = TopHatModel_Gauss_RHS_trunc(x, varargin{:}); label="ToHaGaTrRHS";
    otherwise; y = [];
end
%% -- For Debugging
plot_result = 0;
if plot_result == 1
    % - Initialising the figure object
    figure(); hold on;
    plot(x, y, 'b-', 'linewidth', 2);
    title(sprintf("PES_AtmConcCurve(%s)",label), 'interpreter', 'none'); 
    xlabel(' X ', 'fontweight', 'bold');
    ylabel(' Y ', 'fontweight', 'bold');
    % -- Determining the best limits for the plot
    axis([min(x(:)), max(x(:)), min(y(:)), 1.1*max(y(:))]);
end
end