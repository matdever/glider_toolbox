function [ctd] = load_ctd(csv_root)
%==========================================================================
% DESCRIPTION:
%       [ctd] = load_ctd(csv_root)
%       extracts the CTD data from the glider csv file and calculates:
%           - The Practical Salinity (SP)
%           - The potential Density
%       using the GSW TEOS-10 routine package
%
% INPUTS:
%       'csv_root' is the core filename (or path) of the csv mission file
%
% OUTPUTS:
%       'ctd' is a structure that contains:
%           - Conductivity (S/m)
%           - Latitude (°)
%           - Longitude (°)
%           - Pressure (db)
%           - Depth (m)
%           - Time (matlab datenum format)
%           - In-Situ Temperature (°C ITS-90)
%           - Practical Salinity (SP)
%           - Potential Density
%
% DEPENDENCIES:
%       the GSW TEOS-10 mfiles
%
% AUTHOR:   Mathieu Dever 05-11-2014
% UPDATES:
%==========================================================================

% CTD
% Concatenate the CTD suffix to the core filename
fn = strcat(csv_root,'_water.csv');
clear csv_root

% Gets the header from the file to make sure the right column is used for
% each variables
fid = fopen(fn);
header_tmp = textscan(fid,'%s %s %s %s %s %s %s',1,'delimiter',',');
fclose(fid); clear fid
for ii = 1:length(header_tmp)
    header{ii} = header_tmp{ii}{1};
end; clear ii header_tmp ans

% Read the csv file, skipping the header
M = csvread(fn,1,0);


% Conductivity
temp = strfind(header,'cond');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        ind = find(M(:,ii)>1);
        ctd.cond = M(ind,ii);
    end
end; clear ii temp

% Latitude
temp = strfind(header,'lat');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        ctd.lat = M(ind,ii);
    end
end; clear ii temp

% Longitude
temp = strfind(header,'lon');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        ctd.lon = M(ind,ii);
    end
end; clear ii temp

% Pressure
temp = strfind(header,'pressure');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        ctd.press = M(ind,ii)*10;
    end
end; clear ii temp

% Depth
temp = strfind(header,'depth');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        ctd.depth = M(ind,ii);
    end
end; clear ii temp

% Time
temp = strfind(header,'time');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        ctd.time = datenum(1970,1,1)+M(ind,ii)/(24*60*60);
    end
end; clear ii temp

% In-Situ Temperature
temp = strfind(header,'temp');
for ii = 1:length(header)
    if isempty(temp{ii}) == 0
        ctd.tw = M(ind,ii);
    end
end; clear ii temp

% Practical Salinity
ctd.s = gsw_SP_from_C(ctd.cond*10,ctd.tw,ctd.press); % factor 10 for S/m --> mS/cm

% Potential Density
SA = gsw_SA_from_SP(ctd.s,ctd.press,ctd.lon,ctd.lat); % Absolute Salinity
CT = gsw_CT_from_t(SA,ctd.tw,ctd.press); % Conservative Temperature

ctd.pden = gsw_rho(SA,CT,0) - 1000;

% Issue a warning message if a variable is not created
if isfield(ctd,'cond')~=1
    warning('The conductivity has not been extracted')
end
if isfield(ctd,'lat')~=1
    warning('The latitude has not been extracted')
end
if isfield(ctd,'lon')~=1
    warning('The longitude has not been extracted')
end
if isfield(ctd,'press')~=1
    warning('The pressure has not been extracted')
end
if isfield(ctd,'depth')~=1
    warning('The depth has not been extracted')
end
if isfield(ctd,'time')~=1
    warning('The timestamp has not been extracted')
end
if isfield(ctd,'tw')~=1
    warning('The temperature has not been extracted')
end
if isfield(ctd,'s')~=1
    warning('The practical salinity has not been calculated')
end
if isfield(ctd,'pden')~=1
    warning('The potential density has not been calculated')
end



