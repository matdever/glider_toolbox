function [oxy] = load_oxy_wphase(csv_root)
%==========================================================================
% DESCRIPTION:
%       [oxy] = load_owy_wphase(csv_root)
%       extracts the oxygen data from the optode csv file on the glider
%
% INPUTS:
%       'csv_root' is the core filename (or path) of the csv mission file
%
% OUTPUTS:
%       'oxy' is a structure that contains:
%           - Latitude (°)
%           - Longitude (°)
%           - Depth (m)
%           - Time (matlab datenum format)
%           - Temperature fomr optode
%           - Oxygen concentration (µM)
%           - Oxygen saturation (%)
%
% DEPENDENCIES:
%       N/A
%
% AUTHOR:   Mathieu Dever 05-11-2014
% UPDATES:
%==========================================================================

% CTD
% Concatenate the CTD suffix to the core filename
fn = strcat(csv_root,'_oxy3835.csv');
clear cvs_root

% Because some of the fiesld are text, the textscan function must be used
% for this dataset
fid = fopen(fn);
% Get the header
header_tmp = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s',1,'delimiter',',');
for ii = 1:length(header_tmp)
    header{ii} = header_tmp{ii}{1};
end; clear ii header_tmp ans

% Get the data
M = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','HeaderLines',1,'delimiter',',');
fclose(fid); clear fid ans

% Longitude
temp = strfind(header,'lon');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        oxy.lon = cellfun(@str2num,M{:,ii});
    end
end; clear ii temp

% Latitude
temp = strfind(header,'lat');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        oxy.lat = cellfun(@str2num,M{:,ii});
    end
end; clear ii temp

% Time
temp = strfind(header,'time');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        oxy.time = datenum(1970,1,1) + cellfun(@str2num,M{:,ii}) ./ (24*60*60);
    end
end; clear ii temp

% Depth
temp = strfind(header,'depth');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        oxy.depth = cellfun(@str2num,M{:,ii});
    end
end; clear ii temp

% Temperature from optode - use *_temp field if wphase_temp is empty
% Assign NaN is both variable are empty
temp = strfind(header,'_temp');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0 && strcmp(M{1,ii}{1},'None')==0
        oxy.temp = cellfun(@str2num,M{:,ii});
        break
    elseif isempty(temp{ii}) == 0 && strcmp(M{1,ii}{1},'None')==1 && ii == length(header)
        warning('Temperature from the optode has not been measured')
        oxy.temp = NaN;
    end
end; clear ii temp

% Oxygen concentration - use *_oxygen field if wphase_oxygen is empty
% Assign NaN is both variable are empty
temp = strfind(header,'_oxygen');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0 && strcmp(M{1,ii}{1},'None')==0
        oxy.wphase_o2 = cellfun(@str2num,M{:,ii});
        break
    elseif isempty(temp{ii}) == 0 && strcmp(M{1,ii}{1},'None')==1 && ii == length(header)
        warning('Oxygen concentration has not been measured')
        oxy.wphase_o2 = NaN;
    end
end; clear ii temp

% Oxygen saturation - use *_saturation field if wphase_saturation is empty
% Assign NaN is both variable are empty
temp = strfind(header,'_saturation');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0 && strcmp(M{1,ii}{1},'None')==0
        oxy.wphase_sat = cellfun(@str2num,M{:,ii});
        break
    elseif isempty(temp{ii}) == 0 && strcmp(M{1,ii}{1},'None')==1 && ii == length(header)
        warning('Oxygen saturation has not been measured')
        oxy.wphase_sat = NaN;
    end
end; clear ii temp

% Issue a warning message if a variable is not created
if isfield(oxy,'lat')~=1
    warning('The latitude has not been extracted')
end
if isfield(oxy,'lon')~=1
    warning('The longitude has not been extracted')
end
if isfield(oxy,'depth')~=1
    warning('The depth has not been extracted')
end
if isfield(oxy,'time')~=1
    warning('The timestamp has not been extracted')
end
if isfield(oxy,'temp')~=1
    warning('The temperature from the optode has not been extracted')
end
if isfield(oxy,'wphase_o2')~=1
    warning('The oxygen concentration has not been extracted')
end
if isfield(oxy,'wphase_sat')~=1
    warning('The oxygen saturation has not been extracted')
end