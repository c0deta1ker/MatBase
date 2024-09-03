%% MATLAB Digitisation of Materials Properties Database (MPD)
% The Material Properties Database (MPD) local MATLAB reference file 
% 'MPD_PCC.mat' is created and saved here.
%
% The data that forms the materials property database is taken from a 
% range of sources, where the 'average' values are used for parameters 
% that had more than 1 unique value;
%   [1] https://www.wolframalpha.com/ (electronegativity, electron affinity, ionisation energies, temperatures and crystal structures)
%   [2] https://www.schoolmykids.com/learn/periodic-table-of-elements/ (temperatures, crystal structures, unit cell parameters, electronic and magnetic properties).
%   [3] 10.1063/1.3253115 (Compilation of Energy Band Gaps in Elemental and Binary Compound Semiconductors and Insulators)
%   [4] https://en.wikipedia.org/ (Band-gap estimates, temperatures, crystal structures, unit cell parameters, electronic and magnetic properties)
%   [5] https://pubchem.ncbi.nlm.nih.gov/ (Atomic weights of elements and compounds)
%
% In order to connect a Microsoft Access database (*.mdb, *.accdb) to MATLAB, 
% you must ensure that you have the correct ODBC connected for Windows and 
% a step-by-step guide can be found here: 
%   https://ch.mathworks.com/help/releases/R2017a/database/ug/microsoft-access-odbc-windows.html
%
% Just in case you don't have the correct ODBC drivers, I found that 
% installing these three drivers from Microsoft made the MATLAB Database Explorer work fine;
%   https://docs.microsoft.com/en-us/sql/ssms/download-sql-server-management-studio-ssms?redirectedfrom=MSDN&view=sql-server-ver15
%   https://docs.microsoft.com/en-us/sql/connect/odbc/download-odbc-driver-for-sql-server?view=sql-server-ver15
%   https://www.microsoft.com/en-us/download/details.aspx?id=13255
%
% Step-by-step guide to connecting the database to MATLAB: 
%   https://ch.mathworks.com/help/database/ug/mysql-native-interface-for-windows.html
close all; clear all;
path_matbase    = what('MatBase'); path_matbase = string(path_matbase.path);
path_data       = path_matbase + "\MatBase-v1.2\Material Properties Database (MPD)\data\";
%% (1)     :    Connecting to the MPD
conn            = database('MPD_PCC','admin','admin');
MPD_version     = "MPD_v1";
MPD_PCC         = sqlread(conn,MPD_version); close(conn);
MPD_PCC
%% (2)     :    Saving the latest version of the MPD
file_name = "MPD_PCC"; save_fullfile = path_data + file_name;
save(char(save_fullfile + ".mat"), 'MPD_PCC');
%% (3)     :    Testing that the MPD works
material_props = get_mpd_props('U'); 
material_props
%% (4)     :    Saving the MPD as a table
all_fig = findall(0, 'type', 'figure'); close(all_fig);
fig = uifigure();
fig.Position(3) = 1000;
fig.Position(4) = 500;
fig.Name = sprintf("Materials Properties Database - %s", datetime);
uit = uitable(fig,"Data", MPD_PCC, "Position",[20 20 950 450]);
saveas(fig,path_data+"MPD_PCC.fig");
