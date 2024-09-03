function [be, cls] = calc_be(element, corelevel, formalism, plot_results)
% [be, cls] = calc_be(element, corelevel, formalism, plot_results)
%   This is a general function that extracts the electron binding energies 
%   from elements with Z from 1 to 98 of the individual subshells. The
%   values are from different sources in the literature. The Moulder (1993) [1],  
%   Trzhaskovskaya (2018) [2] and Cant (2022) [3] formalisms are currently
%   available. The Constantinou (2023) database is a combination of all
%   core-level data from all previous work. The user can define the element and core-level.
%
%   IN:
%   -   element:    	string of the element; e.g. "H", "He", "Si", "In"...
%   -   corelevel:      string or M×1 string vector of the core-levels to be probed; e.g. "5s1", "5p1", "5p3", "5d3", "5d5", "5f5', "5f7"...
%   -   formalism:      string of the photoionization cross-section formalism to use. Default:"Constantinou2023" ["Moulder1993", "Trzhaskovskaya2018", "Cant2022", "Constantinou2023"]
%   -   plot_results:   if 1, will plot figure summary, otherwise it wont.
%
%   OUT:
%   -   be:             M×1 vector of the binding energies of the chosen core-levels [eV]. Returns [] for undefined core-level energies.
%   -   cls:            M×1 vector of the core-level labels. Returns [] for undefined core-level energies.
%
%   SEE REFERENCES:
%   [1] John F. Moulder, Handbook of X-ray Photoelectron Spectroscopy, 1993.
%   [2] M.B. Trzhaskovskaya, V.G. Yarzhemsky. Atomic Data and Nuclear Data Tables 119 (2018) 99–174. Web. doi:10.1016/j.adt.2017.04.003
%   [3] David J. H. Cant, Ben F. Spencer, Wendy R. Flavell, Alexander G. Shard. Surf Interface Anal. 2022; 54(4): 442-454. doi:10.1002/sia.7059

%% Default parameters
if nargin < 2; corelevel = [];  end
if nargin < 3; formalism = "Constantinou2023";  end
if nargin < 4; plot_results = 0;  end
if isempty(corelevel); corelevel = []; end
if isempty(formalism); formalism = "Constantinou2023"; end
if isempty(plot_results); plot_results = 0; end
%% Validity checks on the input parameters
corelevel   = string(corelevel);
% -- Ensuring the inputs are in the form of 1D column vectors
if size(corelevel, 2) > 1;  corelevel = corelevel'; end
%% - 1 - Determination of the photoionization cross section sigma
% -- Moulder1993 formalism
if strcmpi(formalism, "Moulder1993") || strcmpi(formalism, "Mou") || strcmpi(formalism, "M") || strcmpi(formalism, "1993")  || strcmpi(formalism, "M1993")
   [be, cls] = be_moulder1993(element, corelevel, plot_results);
    % --- Search other databases if Eb in NaN
    for i = 1:length(be)
        if isnan(be(i))
            try be(i) = be_trzh2018(element, corelevel, 0); if ~isnan(be(i)); fprintf('Eb was NaN, but a finite value was found from the Trzhaskovskaya2018 database.\n'); end
            catch
                try be(i) = be_cant2022(element, corelevel, 0); if ~isnan(be(i)); fprintf('Eb was NaN, but a finite value was found from the Cant2022 database.\n'); end
                catch; fprintf('Eb not found in any database.\n');
                end
            end
        end
    end
% -- Trzhaskovskaya2018 formalism
elseif strcmpi(formalism, "Trzhaskovskaya2018") || strcmpi(formalism, "Trz") || strcmpi(formalism, "T") || strcmpi(formalism, "2018")  || strcmpi(formalism, "T2018")
   [be, cls] = be_trzh2018(element, corelevel, plot_results);
    % --- Search other databases if Eb in NaN
    for i = 1:length(be)
        if isnan(be(i))
            try be(i) = be_moulder1993(element, corelevel, 0); if ~isnan(be(i)); fprintf('Eb was NaN, but a finite value was found from the Moulder1993 database.\n'); end
            catch
                try be(i) = be_cant2022(element, corelevel, 0); if ~isnan(be(i)); fprintf('Eb was NaN, but a finite value was found from the Cant2022 database.\n'); end
                catch; fprintf('Eb not found in any database.\n');
                end
            end
        end
    end
% -- Cant2020 formalism
elseif strcmpi(formalism, "Cant2022") || strcmpi(formalism, "Cant") || strcmpi(formalism, "C") || strcmpi(formalism, "2022")  || strcmpi(formalism, "C2022")
   [be, cls] = be_cant2022(element, corelevel, plot_results);
    % --- Search other databases if Eb in NaN
    for i = 1:length(be)
        if isnan(be(i))
            try be(i) = be_trzh2018(element, corelevel, 0); if ~isnan(be(i)); fprintf('Eb was NaN, but a finite value was found from the Trzhaskovskaya2018 database.\n'); end
            catch
                try be(i) = be_moulder1993(element, corelevel, 0); if ~isnan(be(i)); fprintf('Eb was NaN, but a finite value was found from the Moulder1993 database.\n'); end
                catch; fprintf('Eb not found in any database.\n');
                end
            end
        end
    end
% -- Constantinou2023 formalism
elseif strcmpi(formalism, "Constantinou2023") || strcmpi(formalism, "Constantinou") || strcmpi(formalism, "2023")  || strcmpi(formalism, "C2023")
    [be, cls] = be_constantinou2023(element, corelevel, plot_results);
else; msg = 'Formalism not found. One of the following must be used: "Moulder1993", "Trzhaskovskaya2018", "Cant2022" or "Constantinou2023".'; error(msg);
end
%% Validity check on the outputs
cls(isnan(be)) = [];
be(isnan(be)) = [];
if isempty(be); be = []; end
end