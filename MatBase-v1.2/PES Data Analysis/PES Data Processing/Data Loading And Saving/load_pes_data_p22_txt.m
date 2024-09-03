function dataStr = load_pes_data_p22_txt(PathName, FileName)
% dataStr = load_pes_data_p22_txt(PathName, FileName)
%   Loads a text file from the P22 beamline that consists of many columns. 
%   The first column is the binding energy, all intermediate columns are 
%   spectral intensities from each consecutive sweep and the final column 
%   is the average spectral intensity over all sweeps.
%
%   IN:
%   -   PathName:       string of the full directory path to the data file
%   -   FileName:       string of the filename of the data file to be loaded
%
%   OUT:
%   -   dataStr:        MATLAB data structure for PES data

%% Default parameters
if nargin < 1; PathName = ''; end
if nargin < 2; FileName = ''; end
if isempty(FileName);   FileName = '';  end
if isempty(PathName);   PathName = '';  end
disp('Loading in P22 text-file data...')
%% 1 - Loading and reading in the .txt file
data_table          = readtable(string(PathName) + string(FileName));
dataStr.TimeStamp    = datetime;
dataStr.PathName     = PathName;
dataStr.FileName     = FileName;
dataStr.Type         = "PES";
dataStr.hv           = [];
dataStr.thtM         = [];
dataStr.sweeps       = size(data_table(:,2:end-1), 2);
dataStr.ydat_sweeps  = table2array(data_table(:,2:end-1));
dataStr.xdat         = -1.*round(table2array(data_table(:,1)),3);
dataStr.ydat         = round(mean(table2array(data_table(:,end)),2),3);
dataStr.xdat_lims    = round([min(dataStr.xdat(:)), max(dataStr.xdat(:))],3);
dataStr.ydat_lims    = round([min(dataStr.ydat(:)), max(dataStr.ydat(:))],3);
end