function imfp = calc_imfp(ke_dat, formalism, args)
% imfp = calc_imfp(ke_dat, formalism, args)
%   This is a general function that calculates the electron inelastic mean 
%   free path (IMFP) from different sources in the literature. The 
%   Universal ([1] Seah1979), TPP-2M ([2-4] Tanuma1994), S1 & S2 
%   ([5] Seah2011), S3 & S4 ([6] Seah2012) and JTP ([7] JTP2023) formalisms 
%   are available. The user can define a scalar or vector of kinetic energies for the input.
%   If the args is a string of an element / material, it will look it up 
%   the relevant materials parameters in the Material Properties Database 
%   (MPD) ('MPD_PCC.mat'). Otherwise, the args can be manually inserted as
%   defined below. Note, the IMFP is the average distance that an electron 
%   with a given energy travels between successive inelastic collisions.
%
%   IN:
%   -   ke_dat:  	N×1 column vector of the input electron kinetic energy (for PES; KE = BE - HV) [eV]
%   -   formalism:  string for imfp calculator formalism. Default:"S2" ["Universal","TPP2M","TPP2M-avg","Optical","S1","S2","S3","S4","JTP"]
%   -   args:       string or vector of scalars defining the arguments of the chosen formalism;
%                       (1) string of element / material whose parameters are found in materials database; e.g. "Si", "SiO2", "Al2O3"...
%                       (2) vector with manual entry of material parameters:
%                            -> Optical:        no args, database look-up from NIST 1999.
%                            -> Universal:      no args, material independent.
%                            -> TPP2M-avg:      no args, material independent.
%                            -> TPP2M:  4x1     [density(g/cc),atomicweight(amu),egap(eV),valency(valence electrons per atom)]
%                            -> S1:     5x1     [density(g/cc),atomicweight(amu),egap(eV),Z(atomic mass number or average for compounds),stoichiometry]
%                            -> S2:     1x1     [Z(atomic mass number or average for compounds)]
%                            -> S3:     4x1     [density(g/cc),atomicweight(amu),Z(atomic mass number or average for compounds),stoichiometry]
%                            -> S3-organic:     no args, material independent.
%                            -> S4:     1x1     [Z(atomic mass number or average for compounds)]
%                            -> JTP:    4x1     [density(g/cc),atomicweight(amu),egap(eV),valency(valence electrons per atom)]
%
%   OUT:
%   -   imfp:       N×1 column vector of the electron IMFP values [Angstroms]
%
%   SEE REFERENCES:
%       [1] M. P. Seah, Quantitative electron spectroscopy of surfaces A Standard Data Base for Electron Inelastic Mean Free Paths in Solids (1979)
%       [2] S. Tanuma, Calculations of Electron Inelastic Mean Free Paths. V. Data for 14 Organic Compounds over the 50-2000 eV Range (1994)
%       [3] S. Tanuma, Calculation of electron inelastic mean free paths (IMFPs) VII. Reliability of the TPP-2M IMFP predictive equation (2003)
%       [4] S. Tanuma, Calculations of electron inelasticmean free paths. IX. Data for 41 elemental solids over the 50 eV to 30 keV range (2011)
%       [5] M. P. Seah, An accurate and simple universal curve for the energy-dependent electron inelastic mean free path (2011)
%       [6] M. P. Seah, Simple universal curve for the energy‐dependent electron attenuation length (2012)
%       [7] Jablonski A, Tanuma S, Powell CJ. Surf Interface Anal. 2023; 55(8): 609-637. doi:10.1002/sia.7217
%
%   EXAMPLE #1: Calculatng TPP2M IMFP for Silicon at 500 eV elecron energy.
%   Using the Materials Properties Database:            imfp = calc_imfp(500, "tpp2m", "Si")
%   If the material parameters are not known, then:     imfp = calc_imfp(500, "tpp2m", [2.33,28.0850,1.107,4])

%% Default parameters
% - Default formalism
if nargin < 3; args = []; end
if nargin < 2; formalism = "S2"; end
if isempty(formalism); formalism = "S2"; end
if isempty(args); args = []; end
% - Validity check on the inputs
if ischar(formalism); formalism = string(formalism); end
if ischar(args); args = string(args); end
%% - 1 - Determination of the IMFP
% -- Optical NIST1999 Database
if strcmpi(formalism, "Optical") || strcmpi(formalism, "Opt") || strcmpi(formalism, "NIST")
    if isstring(args)
        material = args;
        [imfp, ~] = imfp_optical(ke_dat, material);
    else; msg = 'Material could not be identified. Only use elements 1 - 92; H, He, Li, Be..., Pa, U'; error(msg);
    end
% -- Universal formalism (1979 Seah)
elseif strcmpi(formalism, "Universal") || strcmpi(formalism, "Uni") || strcmpi(formalism, "U") || strcmpi(formalism, "1979Seah") || strcmpi(formalism, "Universal Formalism")
    imfp = imfp_universal(ke_dat);
% -- TPP2M formalism (1994 Tanuma)
elseif strcmpi(formalism, "TPP2M") || strcmpi(formalism, "TPP-2M") || strcmpi(formalism, "1994Tanuma") || strcmpi(formalism, "TPP2M Formalism")
    if isstring(args)
        material = args;
        imfp = imfp_tpp2m_mpd(ke_dat, material);
    else
        rho = args(1); Nv=args(4); M = args(2); Egap = args(3);
        imfp = imfp_tpp2m(ke_dat, rho, Nv, M, Egap);
    end
% -- Average TPP2M formalism (1994 Tanuma)
elseif strcmpi(formalism, "TPP2M-average") || strcmpi(formalism, "TPP2M-avg") || strcmpi(formalism, "TPP2M-mean")
    imfp = imfp_tpp2m_avg(ke_dat);
% -- S1 formalism (2011 Seah)
elseif strcmpi(formalism, "S1") || strcmpi(formalism, "2011Seah-S1")  || strcmpi(formalism, "S1 Formalism")
    if isstring(args)
        material = args;
        imfp = imfp_S1_mpd(ke_dat, material);
    else
        rho = args(1); M = args(2); Egap = args(3); Z=args(4); stoichiometry=args(5);
        imfp = imfp_S1(ke_dat, rho, M, Egap, Z, stoichiometry);
    end
% -- S2 formalism (2011 Seah)
elseif strcmpi(formalism, "S2") || strcmpi(formalism, "2011Seah-S2")  || strcmpi(formalism, "S2 Formalism")
    if isstring(args)
        material = args;
        imfp = imfp_S2_mpd(ke_dat, material);
    else
        Z=args(1);
        imfp = imfp_S2(ke_dat, Z);
    end
% -- S3 formalism (2012 Seah)
elseif strcmpi(formalism, "S3") || strcmpi(formalism, "2012Seah-S3")
    if isstring(args)
        material = args;
        imfp = eal_S3_mpd(ke_dat, material);
    else
        rho = args(1); M = args(2); Z=args(3); stoichiometry=args(4);
        imfp = eal_S3(ke_dat, rho, M, Z, stoichiometry);
    end
% -- S3 formalism for organics (2012 Seah)
elseif strcmpi(formalism, "S3-organic") || strcmpi(formalism, "S3-organics") || strcmpi(formalism, "S3-o") || strcmpi(formalism, "S3O")
    imfp = eal_S3_organic(ke_dat);
% -- S4 formalism (2012 Seah)
elseif strcmpi(formalism, "S4") || strcmpi(formalism, "2012Seah-S4")
    if isstring(args)
        material = args;
        imfp = eal_S4_mpd(ke_dat, material);
    else
        Z=args(1);
        imfp = eal_S4(ke_dat, Z);
    end
% -- JTP formalism (2023 JTP)
elseif strcmpi(formalism, "JTP") || strcmpi(formalism, "2023JTP") || strcmpi(formalism, "JTP2023")
    if isstring(args)
        material = args;
        imfp = imfp_jtp_mpd(ke_dat, material);
    else
        rho = args(1); Nv=args(4); M = args(2); Egap = args(3);
        imfp = imfp_jtp(ke_dat, rho, Nv, M, Egap);
    end
end
%% -- Validity check on outputs
imfp(isnan(imfp)) = 0; if size(imfp, 2) >1; imfp = imfp'; end                  % Ensuring the imfp is a column vector
end