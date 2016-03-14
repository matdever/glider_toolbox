% This is a small piece of code provided as an example as how the function
% can be used.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTS THE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Imports the data from glider csv files. Including a map of the glider
% track for the imported mission
glider_data_csv2mat('private/2012-05-24_view_sci','plot','yes');


%%%%%%%%%%%%%%%%%%%% Extracts upcast from section %%%%%%%%%%%%%%%%%%%%%%%%%

% Extracts the upcasts and downscasts of the temperature and salinity field
clear
load test_case/test_dataset.mat
UP_AND_DOWNS = glider_profiler(cat(2,ctd.tw,ctd.s),...
    ctd.depth,...
    'plot','yes');


%%%%%%%%%%%%%%% Extracts sections from the mission %%%%%%%%%%%%%%%%%%%%%%

S = glider_sectionner(ctd.lon,ctd.lat);

%%%%%%%%%%%%%%%%%%%% Grids the data from section %%%%%%%%%%%%%%%%%%%%%%%%%

% Extracts the upcasts and downscasts of the temperature and salinity field
[gridded.T,gridded.X,gridded.Z,gridded.counter] =...
    glider_gridder(ctd.tw,ctd.lon,ctd.lat,ctd.depth,'counter','on');
figure
pcolor(gridded.X,-gridded.Z,gridded.T'); shading flat