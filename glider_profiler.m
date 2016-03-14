function cast_struct = glider_profiler(var,depth,varargin)

% glider_profiler  Separates upcasts from downcasts in a glider mission
%===================================================================
%
% USAGE:  cast_struct = glider_profiler(var,depth,varargin);
%
% DESCRIPTION: This routine was developed to automatically separate
%       downcasts from upcasts, based on the pressure gradient. It returns
%       a structure that includes a matrix of each profiles.
%
% INPUT:
%       - var is the variable that wants to be splitted in upcasts and
%           downcasts (e.g. Temperature, Salinity, ...). var can be a
%           vector OR a matrix, but one of its dimensions must match the
%           length of the depth vector.
%       - depth is the pressure or depth field used to find inflexion
%           points. Must be a vector in meters, or decibars.
%       - varargin are the optional arguments available:
%               - (...,'plot','yes') will produce a plot with the inflexion
%                   points
%
% OUTPUT:
%       - cast_struct is a structure that includes 2D matrices where 1
%       dimensions is the number of measurements and second dimension is
%       the cast #.
%
%       It includes a matrix with the pressure associated to each
%       measurements, as well as indivual profiles of var
%
% AUTHOR:   Mathieu Dever 08-04-2015
%
% REFERENCE:
%
% UPDATES:
%       - Mathieu Dever on 30-04-2015: Add an abs() on lines 80 and 98. It
%          makes sure that even if the gradient is negative, it finds
%          datapoints separated by +/- 10 m.
%       - Mathieu Dever & Tetjane Ross on 01-05-2015: Introduced a smoothed
%           depth vector to filter out short term variations in depth
%           gradient. Degree of smoothing is linked to the threshold value
%           used throughout the code (default = 10m)
%       - Mathieu Dever & Tetjana Ross on 12-05-2015: Corrected a typo on
%           line 181, and added an IF-LOOP that produces the depth field
%           only once, as it is the same for every variables
%==================================================================

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------

% Check the number of inputs
if nargin <2
    error(['glider_caster.m: wrong number of arguments, at least 2 inputs',...
        '(variable and pressure or depth field) are required'])
end


% Optional inputs.
vin = varargin;
for ii = 1:2:length(vin)
    % If plot is required
    if isequal(vin{ii},'plot')
        plot_request = vin{ii+1};
    else
        error('-plot is the only optional argument available... for now')
    end
end; clear ii

% If optional argument is not specified, assumes 'no'
if nargin == 2
    plot_request = 'no';
end

%----------------------
% CONSTANTS DEFINITION AND CHECKS
%----------------------

% Check that "depth" is a vector
if isvector(depth)==0
    error('glider_caster.m: "depth" must be a vector')
elseif iscolumn(depth)==0
    depth = depth';
end


[M,N] = size(var);
% Check if "var" is a vector or a matrix
if ismatrix(var)==1
    
    % Time as 1st dimension
    if N == length(depth)
        var = var';
        
        % Returns an error if none of "var"'s dimensions match the length of
        % "depth"
    elseif M ~= length(depth) && N ~= length(depth)
        error(['glider_caster.m: one of the dimensions of "var" must match',...
            ' the length of the "depth" vector'])
    end
    
else
    % Error if "var" has more than 2 dimensions
    error('glider_caster.m: "var" must be a vector or 2D matrix')
end


% Threshold used throughout the code (in meters)
threshold = 10; %10


%----------------------
% CORE CODE
%----------------------

% Gets rid of all measurements made during surfacings (shallower than 1.5m)
var(abs(depth) < 1.5,:) = [];
depth(abs(depth) < 1.5) = [];

% Stores the input depth
depthvar = depth;
% Smooth the depth
depth = smooth(depthvar,threshold/nanmean(abs(diff(depthvar)))/2);

% Separates downcast from upcast
ind_up = find(diff(depth)<0);
ind_down = find(diff(depth)>0);

% Considers boundary conditions
% If starts on an upcast
if abs(depth(1))>threshold
    ind = find(abs(diff(abs(depth(ind_up))))>threshold);
    start_up = [1;ind_up(ind+1)];
    finish_up = ind_up(ind);
    start_down = finish_up+1;
    finish_down = start_up(2:end)-1;
    
    if abs(depth(end))<threshold % If finishes on an upcast
        finish_up = [finish_up;length(depth)];
        
    elseif abs(depth(end))>threshold % If finishes on a downcast
        ind = find(diff(abs(depth(ind_down)))>threshold);
        finish_up = [finish_up;ind_down(ind(end)+1)];
        start_down = [start_down;ind_down(ind(end)+1)+1];
        finish_down = [finish_down;length(depth)];
    end
    
    % If starts on a downcast
elseif abs(depth(1))<threshold
    ind = find(abs(diff(abs(depth(ind_down))))>threshold);
    start_up = ind_down(ind)+1;
    finish_up = ind_down(ind+1)-1;
    start_down = [1;finish_up];
    finish_down = start_up;
    
    if abs(depth(end))<threshold % If finishes on an upcast
        ind = find(diff(abs(depth(ind_up)))>threshold);
        finish_down = [finish_down;ind_up(ind(end)+1)-1];
        start_up = [start_up;ind_up(ind(end)+1)];
        finish_up = [finish_up;length(depth)];
        
    elseif abs(depth(end))>threshold % If finishes on a downcast
        finish_down = [finish_down;length(depth)];
    end
    
end

% NOTE:
% Threshold in meters can be converted to a threshold in the number of
% observations, if the vertical resolution is known:
% threshold in meters / average vertical resolution = threshold in number
% of observation.

for numvar = 1:N % For each variable
    
    % If any cast is longer than the threshold (i.e. if upcasts exist).
    if max(finish_up-start_up) > threshold/nanmean(abs(diff(depthvar)))
        
        % Creates the empty matrices
        cast_struct.upcast{1,numvar} = NaN*ones(max(finish_up - start_up)+1,length(start_up));
        cast_struct.upcast_depth = NaN*ones(max(finish_up - start_up)+1,length(start_up));
        
        for ii = 1:length(start_up)
            cast_struct.upcast{1,numvar}(1:finish_up(ii)-start_up(ii)+1,ii) = var(start_up(ii):finish_up(ii),numvar);%%TARR substituted numvar for N (May 11, 2015)
            if numvar==N %%TARR added if statement (May 11, 2015)
                cast_struct.upcast_depth(1:finish_up(ii)-start_up(ii)+1,ii) = depthvar(start_up(ii):finish_up(ii));
            end
        end
        
    end
    
    % If any cast is longer than the threshold (i.e. if downcasts exist).
    if max(finish_down-start_down) > threshold/nanmean(abs(diff(depthvar)))
        
        cast_struct.downcast{1,numvar} = NaN*ones(max(finish_down - start_down)+1,length(start_down));
        cast_struct.downcast_depth = NaN*ones(max(finish_down - start_down)+1,length(start_down));
        
        for ii = 1:length(start_down)
            
            cast_struct.downcast{1,numvar}(1:finish_down(ii)-start_down(ii)+1,ii) = var(start_down(ii):finish_down(ii),numvar);%%TARR substituted numvar for N (May 11, 2015)
            if numvar==N %%TARR added if statement (May 11, 2015)
                cast_struct.upcast_depth(1:finish_up(ii)-start_up(ii)+1,ii) = depthvar(start_up(ii):finish_up(ii));
            end
        end
    end
end


% PLOTS

% Still in development for better visualization
if isequal(plot_request,'yes')
    figure
    % Plot upcasts
    subplot(2,1,1)
    if isfield(cast_struct,'upcast')
        h0 = plot(ind_up,depth(ind_up),'ok'); axis ij
        hold on
        title(['Upcasts (',num2str(size(cast_struct.upcast{1},2)),' profiles)'])
        % Upcasts end members
        h1 = plot(start_up,depth(start_up),'sr','markerfacecolor','r');
        h2 = plot(finish_up,depth(finish_up),'og','markerfacecolor','g');
        legend([h0 h1 h2],'upcast','start upcast','end upcast','location','best')
    else
        text(0.25,0.5,'No upcast profiles','fontsize',14)
    end
    
    % Plot downcasts
    subplot(2,1,2)
    if isfield(cast_struct,'downcast')
        h0 = plot(ind_down,depth(ind_down),'ok'); axis ij
        hold on
        title(['Downcasts (',num2str(size(cast_struct.downcast{1},2)),' profiles)'])
        % Downcasts end members
        h1 = plot(start_down,depth(start_down),'sr','markerfacecolor','r');
        h2 = plot(finish_down,depth(finish_down),'og','markerfacecolor','g');
        legend([h0 h1 h2],'downcast','start downcast','end downcast','location','best')
    else
        text(0.25,0.5,'No downcast profiles','fontsize',14)
    end
end