function [flsv2] = load_flsv2(csv_root)
%==========================================================================
% DESCRIPTION:
%       [flsv2] = load_flsv2(csv_root)
%       extracts the fluorescence (v2) data from the glider csv file
%
% INPUTS:
%       'csv_root' is the core filename (or path) of the csv mission file
%
% OUTPUTS:
%       'flsv2' is a structure that contains:
%           - Latitude (°)
%           - Longitude (°)
%           - Depth (m)
%           - Time (matlab datenum format)
%           - Chlorophyll (mg.m-3)
%           - Backscatter at 470 µm (m-1.sr-1)
%           - Backscatter at 532 µm (m-1.sr-1)
%
% DEPENDENCIES:
%       N/A
%
% AUTHOR:   Mathieu Dever 05-11-2014
% UPDATES:
%==========================================================================

% CTD
% Concatenate the CTD suffix to the core filename
fn = strcat(csv_root,'_bb2flsv2.csv');
clear cvs_root

% Gets the header from the file to make sure the right column is used for
% each variables
fid = fopen(fn);

% Takes the change in format of the source files into consideration
if datenum(csv_root(1:10),'yyyy-mm-dd') < datenum('2012-05-24','yyyy-mm-dd')
    header_tmp = textscan(fid,'%s %s %s %s %s %s %s',1,'delimiter',',');
else
    header_tmp = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s',1,'delimiter',',');
end
fclose(fid); clear fid

for ii = 1:length(header_tmp)
    header{ii} = header_tmp{ii}{1};
end; clear ii header_tmp ans

% Read the csv file, skipping the header
M = csvread(fn,1,0);

% Checks if the header and the datasets have the same dimension
if length(header)~=size(M,2)
    error('The dimensions of the dataset and the corresponding header do not match')
end

% Longitude
temp = strfind(header,'lon');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        flsv2.lon = M(:,ii);
    end
end; clear ii temp

% Latitude
temp = strfind(header,'lat');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        flsv2.lat = M(:,ii);
    end
end; clear ii temp

% Time
temp = strfind(header,'time');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        flsv2.time = datenum(1970,1,1) + M(:,ii) ./ (24*60*60);
    end
end; clear ii temp

% Depth
temp = strfind(header,'depth');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        flsv2.depth = M(:,ii);
    end
end; clear ii temp

% Chlorophyll
temp = strfind(header,'chl');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        flsv2.chl = M(:,ii);
    end
end; clear ii temp

% Backscatter 470 µm
temp = strfind(header,'b470');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        flsv2.b470 = M(:,ii);
    end
end; clear ii temp

% Backscatter 532 µm
temp = strfind(header,'b532');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        flsv2.b532 = M(:,ii);
    end
end; clear ii temp

% Issue a warning message if a variable is not created
if isfield(flsv2,'lat')~=1
    warning('The latitude has not been extracted')
end
if isfield(flsv2,'lon')~=1
    warning('The longitude has not been extracted')
end
if isfield(flsv2,'depth')~=1
    warning('The depth has not been extracted')
end
if isfield(flsv2,'time')~=1
    warning('The timestamp has not been extracted')
end
if isfield(flsv2,'chl')~=1
    warning('The Chlorophyll has not been extracted')
end
if isfield(flsv2,'b470')~=1
    warning('The Backscatter at 470 µm has not been extracted')
end
if isfield(flsv2,'b532')~=1
    warning('The Backscatter at 532 µm has not been extracted')
end