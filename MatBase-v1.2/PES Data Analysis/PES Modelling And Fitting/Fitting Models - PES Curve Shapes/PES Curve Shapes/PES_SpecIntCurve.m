function y = PES_SpecIntCurve(x, Type, BE, INT, FWHM, MR, LSE, LSI, LSW, ASY)
% y = PES_SpecIntCurve(x, Type, BE, INT, FWHM, MR, LSE, LSI, LSW, ASY)
%   Function that evaluates a generic photoelectron spectroscopy (PES)
%   curve, by defining a primary (P) and spin-orbit split (SOS) component.
%
%   IN:
%   -   x:          N×1 (or 1×N) vector of the input domain (binding energy for PES)
%   -   Type:       string of the curve-shape type. Default: ["sGLA"] ("G","L","V","DS","sGL","sGLA","pGL","pGLA")
%   -   BE:      	scalar of the binding energy of PE curve.
%   -   INT:    	scalar of the peak intensity of PE curve.
%   -   FWHM:     	scalar of the FWHM of the PE curve.
%   -   MR:     	scalar of the Mixing Ratio: 0 for pure Gaussian, 1 is for pure Lorentzian.
%   -   LSE:     	scalar of the binding energy of spin-orbit split PE curve.
%   -   LSI:     	scalar of the branching ratio of spin-orbit split PE curve.
%   -   LSW:     	scalar of the additional lorentzian width of spin-orbit split PE curve.
%   -   ASY:     	scalar of the PE curve asymmetry parameter (usually for metallic systems).
%
%   OUT:
%   -   y:          N×1 (or 1×N) vector of the intensity range (spectral intensity for PES)

%% Default parameters
% Default based on inputs
if nargin < 2;  Type	= "sGLA"; end
if nargin < 3;  BE      = 0.00; end
if nargin < 4;  INT     = 1.00; end
if nargin < 5;  FWHM    = 0.25; end
if nargin < 6;  MR      = 0.50; end
if nargin < 7;  LSE     = 0;  end
if nargin < 8;  LSI     = 0;  end
if nargin < 9;  LSW     = 0;  end
if nargin < 10; ASY     = 0;  end
% Default based on empty inputs
if isempty(Type);   Type   = "sGLA"; end
if isempty(BE);     BE      = 0.00; end
if isempty(INT);    INT     = 1.00; end
if isempty(FWHM);   FWHM    = 0.25; end
if isempty(MR);     MR      = 0.50; end
if isempty(LSE);    LSE     = 0; end
if isempty(LSI);    LSI     = 0; end
if isempty(LSW);    LSW     = 0; end
if isempty(ASY);    ASY     = 0; end
%% Validity checks on the input parameters
% if isrow(x); x = x'; end        % -- Ensure x-data is a column vector
if INT < 0; INT = 0; end        % -- If the INT is <0, pad it to 0
if FWHM < 0; FWHM = 0; end      % -- If the FWHM is negative, pad it to zero
if MR < 0; MR = 0; end          % -- If the MR is negative, pad it to zero
if MR > 1; MR = 1; end          % -- If the MR is >1, pad it to 1
if LSI < 0; LSI = 0; end        % -- If the LSI is <0, pad it to 0
if LSW < 0; LSW = 0; end        % -- If the LSW is <0, pad it to 0
if ASY < -1; ASY = -1; end      % -- If the ASY is <-1, pad it to -1
if ASY > 1; ASY = 1; end        % -- If the ASY is >1, pad it to 1
%% - 1 - Determination of the PES Spectral Intensity Curve
switch lower(Type)
    case lower({'Gaussian','Gauss','Gau','G','G - Gaussian'})
        label="Gaussian";
        P_curve     = GaussianModel(x, BE, INT, FWHM);
        SOS_curve 	= GaussianModel(x, BE+LSE, LSI*INT, FWHM+LSW);   
    case lower({'Lorentzian','Lorentz','Lor','L','L - Lorentizan'})
        label="Lorentzian";
        P_curve     = LorentzianModel(x, BE, INT, FWHM);
        SOS_curve 	= LorentzianModel(x, BE+LSE, LSI*INT, FWHM+LSW);
    case lower({'Voigt','V'})
        label="Voigt";
        P_curve     = VoigtModel(x, BE, INT, FWHM, MR);
        SOS_curve 	= VoigtModel(x, BE+LSE, LSI*INT, FWHM+LSW, MR);
    case lower({'DoniachSunjic','Doniach','DonSun','DS','DS - Doniach'})
        label="DoniachSunjic";
        P_curve     = DoniachModel(x, BE, INT, FWHM, ASY);
        SOS_curve 	= DoniachModel(x, BE+LSE, LSI*INT, FWHM+LSW, ASY);
    case lower({'PseudoVoigtSum','PseudoVoigt-sGL','sGL','sGL - Voigt'})
        label="PseudoVoigt-sGL";
        P_curve     = PseudoVoigtModel_sGL(x, BE, INT, FWHM, MR);
        SOS_curve 	= PseudoVoigtModel_sGL(x, BE+LSE, LSI*INT, FWHM+LSW, MR);    
    case lower({'PseudoVoigtSumAsymmetric','PseudoVoigt-sGLA','sGLA','sGLA - Voigt'})
        label="PseudoVoigt-sGLA";
        P_curve     = PseudoVoigtModel_sGLA(x, BE, INT, FWHM, MR, ASY);
        SOS_curve 	= PseudoVoigtModel_sGLA(x, BE+LSE, LSI*INT, FWHM+LSW, MR, ASY);  
    case lower({'PseudoVoigtProduct','PseudoVoigt-pGL','pGL','pGL - Voigt'})
        label="PseudoVoigt-pGL";
        P_curve     = PseudoVoigtModel_pGL(x, BE, INT, FWHM, MR);
        SOS_curve 	= PseudoVoigtModel_pGL(x, BE+LSE, LSI*INT, FWHM+LSW, MR);  
    case lower({'PseudoVoigtProductAsymmetric','PseudoVoigt-pGLA','pGLA','pGLA - Voigt'})
        label="PseudoVoigt-pGLA";
        P_curve     = PseudoVoigtModel_pGLA(x, BE, INT, FWHM, MR, ASY);
        SOS_curve 	= PseudoVoigtModel_pGLA(x, BE+LSE, LSI*INT, FWHM+LSW, MR, ASY);  
    otherwise; error('Curve type not recognised.'); y = []; label="";
end
% - Summing the components together
y = P_curve + SOS_curve;
%% Validity check on the outputs
% if isrow(y); y = y'; end    % -- Ensure y-data is a column vector
y(isnan(y)) = 0;            % -- Ensure all NaN values are zero
%% -- For Debugging
plot_result = 0;
if plot_result == 1
    % - Initialising the figure object
    figure(); hold on;
    % -- Plotting the 1D data
    plot(x, y, 'b-', 'linewidth', 2);
    title(sprintf("PES_SpecIntCurve(%s)",label), 'interpreter', 'none'); 
    xlabel(' X ', 'fontweight', 'bold');
    ylabel(' Y ', 'fontweight', 'bold');
    % -- Determining the best limits for the plot
    axis([min(x(:)), max(x(:)), min(y(:)), 1.1*max(y(:))]);
end
end