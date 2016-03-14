function sections = glider_sectionner(lon,lat)

% glider_sectionner returns the 1st and last indices for each selected
% sections of the glider mission
%===================================================================
%
% USAGE:  glider_sectionner(lon,lat);
%
% DESCRIPTION: This routine extracts the first and last indices
%               corresponding of each section of the glider mission
%               selected by the user
%
% INPUT:
%       - lon is a vector and is the longitude (in decimal degrees)
%       - lat is a vector and is the latitude (in decimal degrees)
%
% OUTPUT:
%       - sections lists the indices starting (1st column) and finishing
%           (2nd column) the section. Dimensions are thus [Nx2], where N is
%           the number of sections.
%
% AUTHOR:   Mathieu Dever 03-03-2016
%
% DEPENDENCIES:
%
% REFERENCE:
%
% UPDATES:
%==================================================================

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------

% Test the number of I/O
if nargin ~= 2
    error('myApp:argChk', 'Wrong number of inputs')
end
if nargout ~= 1
    error('myApp:argChk', 'Wrong number of outputs')
end

% Verify that inputs are vectors
if isvector(lon)==0 || isvector(lat)==0
    error('glider_sectionner: inputs must be vectors')
end

% Verify that inputs have identical lengths
if length(lon)~=length(lat)
    error('glider_sectionner: both inputs must have the same length')
end

%----------------------
% CORE CODE
%----------------------


% This code separates glider missions into transect - manually. It
% uses the minimums in latitude to detect turning points.

% Distance vector (use default location)
dist = GPS2dist(lon(1),lat(1),lon,lat);

% Calculate distance travelled
for xx=1:length(dist)
    if xx==1
        trav(xx)=dist(xx);
    else
        trav(xx)=abs(dist(xx)-dist(xx-1))+trav(xx-1);
    end
end
clear xx

satisfied = 'No';
while strcmp(satisfied,'Yes')~=1
    
    FF = figure;
    plotyy(1:length(dist),dist/1000,1:length(dist),trav/1000);
    xlabel('data point #');ylabel('distance from shore (km)')
    hold on
    
    N = inputdlg('How many sections to extract from this mission?',...
        'Number of transect from mission',...
        1,{'2'});
    
    [x,~] = ginput (str2double(N{1})*2); clear y
    
    x(x<0)=1;
    x(x>length(lat)) = length(lat);
    scatter(round(x),dist(round(x))/1000,35,'or')
    
    
    % Ask the user to select the end-point and start point of the HL transect
    satisfied = questdlg('Are you satisfied?');
end
close(FF)

sections(:,1) = round(x(1:2:end));
sections(:,2) = round(x(2:2:end));

end

