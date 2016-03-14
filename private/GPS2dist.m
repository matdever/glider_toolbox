function [distance]=GPS2dist(ref_lon,ref_lat,lon,lat)

% GPS2dist  Distance between two lat,lon coordinates
%===================================================================
%
% USAGE:  [distance]=GPS2dist(lon1,lat1,lon2,lat2)
%
% DESCRIPTION:
%   Calculate distance between a position and a reference position (e.g. 
%   starting point) on globe using the "haversine formula". 
%
% INPUT:
%    ref_lon      = reference longitude in decimal degrees (+ve E, -ve W) [- 180.. +180]
%    ref_lat      = reference latitude in decimal degrees (+ve N, -ve S) [- 90.. +90]
%   !! If they areleft empty, Devil's Island coordinates are assigned to
%   the reference point
%
%    lat          = latitude in decimal degrees (+ve N, -ve S) [-90..+90]
%    lon          = longitude in decimal degrees (+ve E, -ve W) [-180..+180]
%
% OUTPUT:
%    distance       = distance between positions in m
%
% AUTHOR:   Mathieu Dever 17-11-2011
%
% REFERENCE:
%    http://en.wikipedia.org/wiki/Haversine_formula
%==================================================================

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------

% check the number of inputs
if nargin ~=4
  error('sw_dist.m: wrong number of arguments, 2 pairs of coordinates (lon,lat) are needed')
end

% assign Devil's Island coordinate is reference point left empty
if isempty(ref_lon)==1 && isempty(ref_lat)==1
    ref_lon=-63.4572; ref_lat=44.5824;
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Devils Island Coordinate are used as reference')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end


% check is the reference positions are scalars
if length(ref_lat)~=1 || length(ref_lon)~=1
   error('GPS2dist.m: the reference (lon,lat) must ba a pair of scalar, not vectors or matrices')
end

% check if the positions are vectors, not matrices
[mlat,nlat] = size(lat);
if mlat~=1 && nlat~=1
   error('GPS2dist.m: lat, lon must be vectors.  No matrices allowed')
end

% convert coordinates to radians
ref_lon=deg2rad(ref_lon);
ref_lat=deg2rad(ref_lat);
lon=deg2rad(lon);
lat=deg2rad(lat);

% compute the latitude and longitude differences
dlon=lon-ref_lon;
dlat=lat-ref_lat;

% constants
Rt=6371000; %Earth Radius in m

% compute the distance
distance=2.*Rt.*asin(sqrt((sin(dlat)/2).^2+cos(ref_lat).*cos(lat).*(sin(dlon/2)).^2));