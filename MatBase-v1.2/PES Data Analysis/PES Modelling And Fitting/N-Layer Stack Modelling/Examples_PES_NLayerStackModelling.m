close all; clear all;

%% Modeller 01 - Uniform and Homogenous Layer Stack
%% :   Use the MatBase Application
run App_MatBase_NLayerModeller01;
%% 1    :   Plotting an arbitrary stack
close all; clear all;
% -- Testing with cell array inputs
lyr_mat         = {'SiO2', 'Al', 'Co', 'U', 'InAs'};
lyr_thick       = {1.5, 2.5, 3.5, 4.5, Inf};
view_nlayer_model01_atmconc(lyr_mat, lyr_thick);
view_nlayer_model01_sample(lyr_mat, lyr_thick);
% -- Testing with vector inputs
lyr_mat         = ["SiO2", "Al", "Co", "U", "InAs"];
lyr_thick       = [1.5, 2.5, 3.5, 4.5, Inf];
view_nlayer_model01_atmconc(lyr_mat, lyr_thick);
%% 2    :   Simple 2-layer model
close all; clear all;
lyr_mat         = {'SiO2', 'Si'};               % Defining the material of each layer, from left to right (top-to-bottom on the schematic).
lyr_thick       = {2.0, Inf};                   % Defining the thickness of each layer. For a 2-layer model system, only the upper layer matters, as the last layer is always considered as the bulk.
lyr_ele         = {'Si', 'Si'};                 % Defining the core-level being probed for each layer.
lyr_cls         = {'2p3', '2p3'};               % Defining the core-level being probed for each layer.
% EXTRACTING PHOTON ENERGY DEPENDENCE
hv              = 400:50:5000;
theta           = 0.0;
temp_pes_data    = nlayer_model01_run(lyr_mat, lyr_thick, lyr_ele, lyr_cls, hv, theta);
nlayer_model01_view(temp_pes_data);
% EXTRACTING EMISSION ANGLE DEPENDENCE
hv              = 2500;     % Defining the incident photon energy
theta           = 0:1:90;
temp_pes_data    = nlayer_model01_run(lyr_mat, lyr_thick, lyr_ele, lyr_cls, hv, theta);
nlayer_model01_view(temp_pes_data);
% EXTRACTING A PHOTON ENERGY AND EMISSION ANGLE DEPENDENCE
hv              = 400:400:5000;
theta           = 0:10:90;
temp_pes_data    = nlayer_model01_run(lyr_mat, lyr_thick, lyr_ele, lyr_cls, hv, theta);
nlayer_model01_view(temp_pes_data);
%% 3    :   4-layer model (an ultra-shallow dopant delta-layer in silicon)
close all; clear all;
lyr_mat         = {'SiO2', 'Si', 'As', 'Si'}; 
lyr_thick       = {1.0, 0.5, 1.0, Inf};
lyr_ele         = {'Si', 'Si', 'As', 'Si'};
lyr_cls         = {'2p3', '2p3', '3p3', '2p3'};
% EXTRACTING PHOTON ENERGY DEPENDENCE
hv              = 400:100:5000;
theta           = 0;
pes_model_4L    = nlayer_model01_run(lyr_mat, lyr_thick, lyr_ele, lyr_cls, hv, theta);
nlayer_model01_view(pes_model_4L);

%% Modeller 02 - Uniform and Varying Concentration Layer Stack
%% :   Use the MatBase Application
run App_MatBase_NLayerModeller02;
%% 1    :   Plotting an arbitrary stack
close all; clear all;
lyr_mat         = {{'SiO2'}, {'Si'}};
lyr_conc        = {{100},{100}};
lyr_ele         = {{'Si'}, {'Si'}};
lyr_cls         = {{'2p3'}, {'2p3'}};
lyr_thick       = {2.0, Inf};
view_nlayer_model02_sample(lyr_mat, lyr_thick, lyr_conc);
view_nlayer_model02_atmconc(lyr_mat, lyr_thick, lyr_conc);
%% 2    :   Simple 2-layer model
close all; clear all;
lyr_mat         = {{'SiO2'}, {'Si'}};
lyr_conc        = {{100},{100}};
lyr_ele         = {{'Si'}, {'Si'}};
lyr_cls         = {{'2p3'}, {'2p3'}};
lyr_thick       = {2.0, Inf};
% 2.2    :    EXTRACTING A PHOTON ENERGY DEPENDENCE
hv              = 400:50:5000;
theta           = 0.0;
temp_pes_data    = nlayer_model02_run(lyr_mat, lyr_conc, lyr_ele, lyr_cls, lyr_thick, hv, theta);
nlayer_model02_view(temp_pes_data);
% EXTRACTING EMISSION ANGLE DEPENDENCE
hv              = 2500;
theta           = 0:1:90;
temp_pes_data    = nlayer_model02_run(lyr_mat, lyr_conc, lyr_ele, lyr_cls, lyr_thick, hv, theta);
nlayer_model02_view(temp_pes_data);
%% 3        :    Varying concentration layer
close all;
lyr_mat         = {{'SiO2'}, {'Al'}, {'Al','In'}, {'InAs'}};
lyr_conc        = {{100},{100},{50,50},{100}};
lyr_ele         = {{'Si'},{'Al'},{'Al','In'},{'In'}};
lyr_cls         = {{'2p3'},{'2p3'},{'2p3','3p3'},{'3p3'}};
lyr_thick       = {3, 5, 2, Inf};
view_nlayer_model02_sample(lyr_mat, lyr_thick, lyr_conc);
view_nlayer_model02_atmconc(lyr_mat, lyr_thick, lyr_conc);
% 2.2    :    EXTRACTING A PHOTON ENERGY DEPENDENCE
hv              = 1500:1:8000;
theta           = 0;
temp_pes_data    = nlayer_model02_run(lyr_mat, lyr_conc, lyr_ele, lyr_cls, lyr_thick, hv, theta);
nlayer_model02_view(temp_pes_data);

%% Modeller 03 - Functional Concentration



%% Control Simulations to look for consistency Si all models
hv              = 6000;
tht             = 0:1:90;
%% A    :   Control Simulation of Uniform Stack (Modeller 01)
close all;
% -- Defining the sample stack
lyr_mat         = {"SiO2", "Pb",   "Si"};
lyr_ele         = {"Si",    "Pb",   "Si"};
lyr_cls         = {"2p3",   "4p3",  "2p3"};
lyr_thick       = {3,       5,      Inf};
view_nlayer_model01_atmconc(lyr_mat, lyr_thick);
% -- Running the model simulations
model01 = nlayer_model01_run(lyr_mat, lyr_thick, lyr_ele, lyr_cls, hv, tht, [], [], "Cant2022", "JTP", 0);
% -- Plotting the Si:Pb ratio
In2Pb_model01 = model01.lyr_Inorm{3} ./ model01.lyr_Inorm{2};
figure(); hold on;
plot(tht, In2Pb_model01, 'kx-');
box on;
xlabel('Theta (degrees)', 'FontWeight','bold'); 
ylabel('Si:Pb Ratio (Model 01)', 'FontWeight','bold'); 
axis([0.0, 90., 0, 0.35]); 
% -- Plotting the layer contributions
figure(); hold on;
for i = 1:length(lyr_mat)
    plot(tht, model01.lyr_Inorm{i}, '-', 'linewidth', 2);
end
box on;
xlabel('Theta (degrees)', 'FontWeight','bold'); 
ylabel('Relative Contribution (Model 01)', 'FontWeight','bold'); 
axis([0.0, 90., 0, 1.05]); 
%% B    :   Control Simulation of Uniform Stack (Modeller 02)
close all;
% -- Defining the sample stack
lyr_mat         = {{"SiO2"},   {"Pb"},     {"Si"}};
lyr_ele         = {{"Si"},      {"Pb"},     {"Si"}};
lyr_cls         = {{"2p3"},     {"4p3"},    {"2p3"}};
lyr_conc        = {{100},       {100},      {100}};
lyr_thick       = {3,           5,          Inf};
view_nlayer_model02_atmconc(lyr_mat, lyr_thick, lyr_conc);
% -- Running the model simulations
model02 = nlayer_model02_run(lyr_mat, lyr_conc, lyr_ele, lyr_cls, lyr_thick, hv, tht, [], [], "C2022", "JTP", 0);
% -- Plotting the Si:Pb ratio
In2Pb_model02 = model02.lyr_Inorm{3}{1} ./ model02.lyr_Inorm{2}{1};
figure(); hold on;
plot(tht, In2Pb_model02, 'rx-');
box on;
xlabel('Theta (degrees)', 'FontWeight','bold'); 
ylabel('Si:Pb Ratio (Model 02)', 'FontWeight','bold'); 
axis([0.0, 90., 0, 0.35]); 
% -- Plotting the layer contributions
figure(); hold on;
for i = 1:length(lyr_mat)
    plot(tht, model02.lyr_Inorm{i}{1}, '-', 'linewidth', 2);
end
box on;
xlabel('Theta (degrees)', 'FontWeight','bold'); 
ylabel('Relative Contribution (Model 02)', 'FontWeight','bold'); 
axis([0.0, 90., 0, 1.05]); 
%% C    :   Control Simulation of Uniform Stack (Modeller 03)
close all;
% -- Defining the sample stack
close all;
lyr_mat         = {"SiO2", "Pb",   "Si"};
lyr_ele         = {"Si",    "Pb",   "Si"};
lyr_cls         = {"1s1",   "4p3",  "2p3"};
lyr_type        = {"ToHa",    "ToHaExRHS",   "StLHS"};
lyr_conc        = {100,     100,    100};
lyr_thick       = {3,       5,      Inf};
lyr_cdl         = {0.,      0.,     0.};
view_nlayer_model03_atmconc(lyr_mat, lyr_type, lyr_conc, lyr_thick, lyr_cdl);
% -- Running the model simulations
model03 = nlayer_model03_run(lyr_mat, lyr_ele, lyr_cls, lyr_type, lyr_conc, lyr_thick, lyr_cdl, hv, tht, [], [], "C2022", "JTP", 1);
% -- Plotting the Si:Pb ratio
In2Pb_model03 = model03.lyr_Inorm{3} ./ model03.lyr_Inorm{2};
figure(); hold on;
plot(tht, In2Pb_model03, 'bx-');
box on;
xlabel('Theta (degrees)', 'FontWeight','bold'); 
ylabel('Si:Pb Ratio (Model 03)', 'FontWeight','bold'); 
axis([0.0, 90., 0, 0.35]); 
% -- Plotting the layer contributions
figure(); hold on;
for i = 1:length(lyr_mat)
    plot(tht, model03.lyr_Inorm{i}, '-', 'linewidth', 2);
end
box on;
xlabel('Theta (degrees)', 'FontWeight','bold'); 
ylabel('Relative Contribution (Model 03)', 'FontWeight','bold'); 
axis([0.0, 90., 0, 1.05]); 
%% D    :   Comparing all simulations
% 3.1 Comparing Si:Pb Ratio
figure(); hold on;
plot(tht, In2Pb_model01, 'kd-', 'MarkerSize', 5);
plot(tht, In2Pb_model02, 'rs-', 'MarkerSize', 5);
plot(tht, In2Pb_model03, 'bo-', 'MarkerSize', 5);
box on;
xlabel('Theta (degrees)', 'FontWeight','bold'); 
ylabel('Si:Pb Ratio', 'FontWeight','bold'); 
axis([0.0, 90., 0, 0.025]); 
legend({'Modeller01', 'Modeller02','Modeller03'}, 'location', 'best');
% 3.3 Comparing all PES intensities
figure(); hold on;
cols = {'k', 'b', 'r'};
for i = 1:length(lyr_mat); plot(tht, model01.lyr_Inorm{i}, '-', 'linewidth', 1, 'color', cols{i}); end
for i = 1:length(lyr_mat); plot(tht, model02.lyr_Inorm{i}{1}, 'x:', 'linewidth', 1, 'color', cols{i}); end
for i = 1:length(lyr_mat); plot(tht, model03.lyr_Inorm{i}, '--', 'linewidth', 1, 'color', cols{i}); end
box on;
xlabel('Theta (degrees)', 'FontWeight','bold'); 
ylabel('Relative Intensity', 'FontWeight','bold'); 
axis([0.0, 90., 0, 1.05]); 
