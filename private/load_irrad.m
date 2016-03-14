function [irrad] = load_irrad(csv_root)
%==========================================================================
% DESCRIPTION:
%       [irrad] = load_irrad(csv_root)
%       extracts the irradiance data from the glider csv file 
%
% INPUTS:
%       'csv_root' is the core filename (or path) of the csv mission file
%
% OUTPUTS:
%       'irrad' is a structure that contains:
%           - Latitude (°)
%           - Longitude (°)
%           - Depth (m)
%           - Time (matlab datenum format)
%           - Irradiance #1 (µW cm-2 nm-1)
%           - Irradiance #2 (µW cm-2 nm-1)
%           - Irradiance #3 (µW cm-2 nm-1)
%           - Irradiance #4 (µW cm-2 nm-1)
%
% DEPENDENCIES:
%       N/A
%
% AUTHOR:   Mathieu Dever 05-11-2014
% UPDATES:  
%==========================================================================

% CTD
% Concatenate the CTD suffix to the core filename
fn = strcat(csv_root,'_ocr504i.csv');
clear cvs_root

% Gets the header from the file to make sure the right column is used for
% each variables
fid = fopen(fn);
header_tmp = textscan(fid,'%s %s %s %s %s %s %s %s',1,'delimiter',',');
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
        irrad.lon = M(:,ii);
    end
end; clear ii temp

% Latitude
temp = strfind(header,'lat');
for ii = 1:length(header)
	if isempty(temp{ii}) == 0
        irrad.lat = M(:,ii);
    end
end; clear ii temp

% Time
temp = strfind(header,'time');
for ii = 1:length(header)
	if isempty(temp{ii}) == 0
        irrad.time = datenum(1970,1,1) + M(:,ii) ./ (24*60*60);
    end
end; clear ii temp
    
% Depth
temp = strfind(header,'depth');
for ii = 1:length(header)
	if isempty(temp{ii}) == 0
        irrad.depth = M(:,ii);
    end
end; clear ii temp

% Irradiance #1
temp = strfind(header,'irrad1');
for ii = 1:length(header)
	if isempty(temp{ii}) == 0
        irrad.irrad1 = M(:,ii);
    end
end; clear ii temp

% Irradiance #2
temp = strfind(header,'irrad2');
for ii = 1:length(header)
	if isempty(temp{ii}) == 0
        irrad.irrad2 = M(:,ii);
    end
end; clear ii temp

% Irradiance #3
temp = strfind(header,'irrad3');
for ii = 1:length(header)
	if isempty(temp{ii}) == 0
        irrad.irrad3 = M(:,ii);
    end
end; clear ii temp

% Irradiance #4
temp = strfind(header,'irrad4');
for ii = 1:length(header)
	if isempty(temp{ii}) == 0
        irrad.irrad4 = M(:,ii);
    end
end; clear ii temp


% Issue a warning message if a variable is not created
if isfield(irrad,'lat')~=1
    warning('The latitude has not been extracted')
end
if isfield(irrad,'lon')~=1
    warning('The longitude has not been extracted')
end
if isfield(irrad,'depth')~=1
    warning('The depth has not been extracted')
end
if isfield(irrad,'time')~=1
    warning('The timestamp has not been extracted')
end
if isfield(irrad,'irrad1')~=1
    warning('The Irradiance #1 has not been extracted')
end
if isfield(irrad,'irrad2')~=1
    warning('The Irradiance #2 has not been extracted')
end
if isfield(irrad,'irrad3')~=1
    warning('The Irradiance #3 has not been extracted')
end
if isfield(irrad,'irrad4')~=1
    warning('The Irradiance #4 has not been extracted')
end