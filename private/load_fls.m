function [fls] = load_fls(csv_root)
%==========================================================================
% DESCRIPTION:
%       [fls] = load_fls(csv_root)
%       extracts the fluorescence data from the glider csv file
%
% INPUTS:
%       'csv_root' is the core filename (or path) of the csv mission file
%
% OUTPUTS:
%       'fls' is a structure that contains:
%           - Latitude (°)
%           - Longitude (°)
%           - Depth (m)
%           - Time (matlab datenum format)
%           - CDOM (ppb)
%           - Backscatter at 660 µm (m-1.sr-1)
%           - Backscatter at 880 µm (m-1.sr-1)
%
% DEPENDENCIES:
%       N/A
%
% AUTHOR:   Mathieu Dever 05-11-2014
% UPDATES:
%==========================================================================

% FLS
% Concatenate the FLS suffix to the core filename
fn = strcat(csv_root,'_bb2fls.csv');

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

% Longitude
temp = strfind(header,'lon');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        fls.lon = M(:,ii);
    end
end; clear ii temp

% Latitude
temp = strfind(header,'lat');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        fls.lat = M(:,ii);
    end
end; clear ii temp

% Time
temp = strfind(header,'time');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        fls.time = datenum(1970,1,1) + M(:,ii) ./ (24*60*60);
    end
end; clear ii temp

% Depth
temp = strfind(header,'depth');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        fls.depth = M(:,ii);
    end
end; clear ii temp

% CDOM
temp = strfind(header,'cdom');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        fls.cdom = M(:,ii);
    end
end; clear ii temp

% Backscatter 660 µm
temp = strfind(header,'b660');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        fls.b660 = M(:,ii);
    end
end; clear ii temp

% Backscatter 880 µm
temp = strfind(header,'b880');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        fls.b880 = M(:,ii);
    end
end; clear ii temp

% Issue a warning message if a variable is not created
if isfield(fls,'lat')~=1
    warning('The latitude has not been extracted')
end
if isfield(fls,'lon')~=1
    warning('The longitude has not been extracted')
end
if isfield(fls,'depth')~=1
    warning('The depth has not been extracted')
end
if isfield(fls,'time')~=1
    warning('The timestamp has not been extracted')
end
if isfield(fls,'cdom')~=1
    warning('The CDOM has not been extracted')
end
if isfield(fls,'b660')~=1
    warning('The Backscatter at 660 µm has not been extracted')
end
if isfield(fls,'b880')~=1
    warning('The Backscatter at 880 µm has not been extracted')
end