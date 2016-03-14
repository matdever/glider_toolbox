function [output,X,Z,counter] = glider_gridder(input,lon,lat,depth,varargin)
%
% Moving block technique to grid the data from glider's transect.
%==========================================================================
%
% USAGE:    [output] = glider_grider_v2(input,lon,lat,depth,varargin)
%
% DESCRIPTION:
%   It applies the "moving block" technique to grid the gliders data. There
%   are no weigth to the points within a block, because there are a lot of
%   data per blocks.
%
% INPUTS:
%   input    = Data to grid (Temperature, Salinity, Density, ...)
%   lon      = Longitude of each data point, in decimal degrees
%   lat      = Latitude of each data point, in decimal degrees
%   depth    = Depth of each data point, in meters
%
% OPTIONAL INPUTS (varargin):
%
%       'horiz_grid' - 2 or 3 element vector containing the lower and upper
%       limits of the grid and the resolution to use, in m (default = 1000 m)
%       [min_grid horizontal_rez max_grid]. Resolution could be ommited.
%       default: [0 1000 200000]
%
%       'vert_grid' - 2 or 3 element vector containing the lower and upper
%       limits of the grid and the resolution to use, in m (default = 0.5 m)
%       [min_grid vertical_rez max_grid]. Resolution could be ommited.
%       default: [0 0.5 205]
%
%       'origin' - 2 element vector containing the coordinates of the point
%       of origin used to calculate the distance, in decimal degrees.
%       Default values are the coordinates of Devil's Island in Halifax
%       Harbour ([-63.4572, 44.5824])
%
%       'counter' - 'on' or 'off', this returns a matrix the same size as
%       input giving the number of datapoints averaged from input to obtain
%       the output matrix (Default = 'off')
%
%
% OUTPUTS:
%       output  = Grided input
%       X       = Horizontal vector (from origin)
%       Z       = Vertical vector, aka depth (in m)
%       counter = number of datapoints averaged per cell, if requested in
%       optional inputs
%
% AUTHOR:    Mathieu Dever 24-10-2012
%
% DEPENDENCIES:
%       - GPS2dist.m
%
% REFERENCE:
%
% UPDATES:
%=========================================================================

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------

if nargin < 4
    error('myApp:argChk', 'Inputs must include the data and the (lat,lon,depth)')
end
if nargout ~= 3 && nargout ~= 4
    error('myApp:argChk', 'Wrong number of outputs')
end

vin = varargin;
for ii = 1:2:length(vin)
    % Define the horizontal grid vector
    if isequal(vin{ii},'horiz_grid')
        min_hgrid = vin{ii+1}(1);
        rez_hgrid = vin{ii+1}(2);
        max_hgrid = vin{ii+1}(3);
        
        % QC tests
        if max_hgrid<min_hgrid || length(vin{ii+1}(:))<2
            error('horiz_grid: Please check your definition of the desired grid')
        elseif length(vin{ii+1}(:))==2
            rez_hgrid = 1000; % Default value
        end
        
        % Define the vertical grid vector
    elseif isequal(vin{ii},'vert_grid')
        min_vgrid = vin{ii+1}(1);
        rez_vgrid = vin{ii+1}(2);
        max_vgrid = vin{ii+1}(3);
        
        % QC tests
        if max_vgrid<min_vgrid || length(vin{ii+1}(:))<2
            error('vert_grid:Please check your definition of the desired grid')
        elseif length(vin{ii+1}(:))==2
            rez_vgrid = 0.5; % Default value
        end
        
        % Define the origin point for horizontal axis
    elseif isequal(vin{ii},'origin')
        lon0 = vin{ii+1}(1);
        lat0 = vin{ii+1}(2);
        
        % Define the origin point for horizontal axis
    elseif isequal(vin{ii},'counter')
        if isequal(vin{ii+1},'on');
            counting = 'yes';
        end
        
    else
        error([vin{ii},': Unknown optional input'])
    end
end

%----------------------
% DEFINE DEFAULT VALUES
%----------------------
if exist('min_hgrid','var')==0
    min_hgrid = 0;
    rez_hgrid = 1000;
    max_hgrid = 200000;
end
if exist('min_vgrid','var')==0
    min_vgrid = 0;
    rez_vgrid = 0.5;
    max_vgrid = 205;
end
if exist('lon0','var')==0
    % Devil's Island coordinates
    lon0 = -63.4572;
    lat0 = 44.5824;
end
if exist('counting','var')==0
    counting = 'no';
end

% Distance vector
dist = GPS2dist(lon0,lat0,lon,lat);

% The grid
X = min_hgrid:rez_hgrid:max_hgrid;
Z = min_vgrid:rez_vgrid:max_vgrid;

output = NaN*ones(length(X),length(Z));

for xx=1:size(X,2)
    for zz=1:size(Z,2)
        
        % b.c. corner 1
        if xx==1 && zz==1
            % Find the data points in the moving box
            ind = find(dist<X(xx+1) & dist>X(xx) & depth<Z(zz+1) & depth>Z(zz));
            if isempty(ind)==1
                output(xx,zz) = NaN;
            else
                output(xx,zz)=nanmean(input(ind));
            end
            
            counter(xx,zz) = length(ind);
            
            % b.c. corner 2
        elseif xx==1 && zz==size(Z,2)
            ind=find(dist<X(xx+1) & dist>X(xx) & depth<Z(zz-1) & depth>Z(zz));
            if isempty(ind)==1
                output(xx,zz)=NaN;
            else
                output(xx,zz)=nanmean(input(ind));
            end
            counter(xx,zz) = length(ind);
            
            % b.c. corner 3
        elseif xx==size(X,2) && zz==size(Z,2)
            % Find the data points in the moving box
            ind=find(dist<X(xx-1) & dist>X(xx) & depth<Z(zz-1) & depth>Z(zz));
            if isempty(ind)==1
                output(xx,zz)=NaN;
            else
                output(xx,zz)=nanmean(input(ind));
            end
            counter(xx,zz) = length(ind);
            
            % b.c. corner 4
        elseif xx==size(X,2) && zz==1
            ind=find(dist<X(xx-1) & dist>X(xx) & depth<Z(zz+1) & depth>Z(zz));
            if isempty(ind)==1
                output(xx,zz)=NaN;
            else
                output(xx,zz)=nanmean(input(ind));
            end
            counter(xx,zz) = length(ind);
            
            % b.c. side 1
        elseif zz==1
            % Find the data points in the moving box
            ind=find(dist<X(xx+1) & dist>X(xx-1) & depth<Z(zz+1) & depth>Z(zz));
            if isempty(ind)==1
                output(xx,zz)=NaN;
            else
                output(xx,zz)=nanmean(input(ind));
            end
            counter(xx,zz) = length(ind);
            
            % b.c. side 2
        elseif xx==1
            % Find the data points in the moving box
            ind=find(dist<X(xx+1) & dist>X(xx) & depth<Z(zz+1) & depth>Z(zz-1));
            if isempty(ind)==1
                output(xx,zz)=NaN;
            else
                output(xx,zz)=nanmean(input(ind));
            end
            counter(xx,zz) = length(ind);
            
            % b.c. side 3
        elseif zz==size(Z,2)
            % Find the data points in the moving box
            ind=find(dist<X(xx+1) & dist>X(xx-1) & depth<Z(zz-1) & depth>Z(zz));
            if isempty(ind)==1
                output(xx,zz)=NaN;
            else
                output(xx,zz)=nanmean(input(ind));
            end
            counter(xx,zz) = length(ind);
            
            % b.c. side 4
        elseif xx==size(X,2)
            % Find the data points in the moving box
            ind=find(dist<X(xx-1) & dist>X(xx) & depth<Z(zz+1) & depth>Z(zz-1));
            if isempty(ind)==1
                output(xx,zz)=NaN;
            else
                output(xx,zz)=nanmean(input(ind));
            end
            counter(xx,zz) = length(ind);
            
        else
            % Find the data points in the moving box
            ind=find(dist<X(xx+1) & dist>X(xx-1) & depth<Z(zz+1) & depth>Z(zz-1));
            if isempty(ind)==1
                output(xx,zz)=NaN;
            else
                output(xx,zz)=nanmean(input(ind));
            end
            counter(xx,zz) = length(ind);
            
        end
    end
end

% % % Saves variables in structure
% % temp = output; clear output
% % output.output = temp; clear temp
% % if isequal(counting,'yes')
% % output.counter = counter;
% % end
% % output.X = X;
% % output.Z = Z;