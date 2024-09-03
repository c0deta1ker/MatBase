function pesStr = create_pes_data()
% pesStr = create_pes_data()
%   This function generates a template data structure used for
%   analysing PES data. The fields can be filled in, after which the user
%   can then use the data structure within all the processing / fitting
%   tools.
pesStr              = struct();
pesStr.TimeStamp    = datetime;
pesStr.PathName     = '';
pesStr.FileName     = '';
pesStr.Type         = "PES";
pesStr.hv           = [];
pesStr.thtM         = [];
pesStr.xdat         = [];
pesStr.ydat         = [];
pesStr.xdat_lims    = [];
pesStr.ydat_lims    = [];
end