function ratio = get_ele_ratio_from_mat(material, element)
% ratio = get_ele_ratio_from_mat(material, element)
%   Extracts the elemental ratio of a given element within a material / compund. 
%   For example, Al2O3 has an elemental ratio of 2/5 for Al and 3/5 for O.
%   This is important for elemental analysis.
%
%   IN:
%   -   material:           char/string of the material; e.g. "Si", "SiO2", "Al2O3"...
%   -   element:            char/string of the element of interest in the material
%
%   OUT:
%   -   ratio:              scalar value of the elemental ratio of the element of interest in the material / compound
%
% Examples:     Al = get_ele_ratio_from_mat("Al2O3", "Al"); O = get_ele_ratio_from_mat("Al2O3", "O");
%               Si = get_ele_ratio_from_mat("SiO2", "Si");  O = get_ele_ratio_from_mat("SiO2", "O");
%               In = get_ele_ratio_from_mat("InAs", "In");  As = get_ele_ratio_from_mat("InAs", "As");

%% Validity check on the inputs
material = string(material);
element = string(element);
ratio = [];
%% 1 : Extracting information from the formula
if strcmpi(material, element) == 1; ratio = 1;
else
    str = char(material);
    num=regexp(str,'\d+','match');  % cell array containing the numbers
    D = isstrprop(str,'digit'); %logical array giving location of numbers
    U = isstrprop(str,'upper'); %logical giving location of upper case alphas   
    L = isstrprop(str,'lower'); %logical giving location of lower case alphas
    NumElem = sum(U);   %number of upper case alphas == number of elements in formula
    Formula = struct('element',{},'quantity',{}); %initialize output
    % - If only letters are used in the chemical formula
    if sum(D) == 0
        n = 1;
        for i = 1:NumElem
            Formula(i).element = str(n:n+1);
            Formula(i).quantity = 1;
            n = n+2;
        end
    % - If letters and digits are used in the chemical formula
    else
        % - Loop through formula to extract quantities
        n = 1;
        num_counter = 1;
        for i = 1:NumElem    
            if U(n)
                if U(n) && L(n+1)
                    Formula(i).element = str(n:n+1);
                    n = n+2;
                    if ~D(n)
                        Formula(i).quantity = 1;
                    elseif D(n)
                        Formula(i).quantity = str2num(num{num_counter});
                        n = n+length(num{num_counter});
                        num_counter = num_counter+1;
                    end
                elseif U(n) && ~L(n+1)
                    Formula(i).element = str(n);
                    n = n+1;
                    if D(n)
                        Formula(i).quantity = str2num(num{num_counter});                
                        n = n+length(num{num_counter});
                        num_counter = num_counter+1;      
                    elseif ~D(n)
                        Formula(i).quantity = 1;
                    end
                end       
            else
            n = n+1;
            end
        end
    end
    % - Find the total sum
    total = 0; for i = 1:length(Formula); total = Formula(i).quantity + total; end
    % - Output the ratio of the selected element
    for i = 1:length(Formula); if strcmpi(string(Formula(i).element), string(element)); ratio = Formula(i).quantity ./ total; end; end
end
%% Validity check on output
if isempty(ratio); ratio = 1; end
end