function glider_data_csv2mat(csv_root,varargin)

% glider_data_csv2mat extracts csv-file into a mat-file
%===================================================================
%
% USAGE:  glider_data_csv2mat(csv_root,varargin);
%
% DESCRIPTION: This routine import the glider data from the source
%               csv-files and saves a mat-file. The mat-file includes
%               a matlab structure for each instrument mounted on the
%               glider:
%               - ctd (self-explanatory)
%               - fls for backscatter sensor
%               - flsv2 for fluorometer sensor
%               - irrad for irrandiance sensor
%               - oxy for the oxygen sensor
%
% INPUT:
%       - csv_root is the prefix used for each csv-file corresponding to a
%           single mission [e.g. '2012-05-24_view_sci']
%       - varargin are the optional arguments available:
%               - (...,'plot','yes') will produce a plot of the glider
%               track.
%
% OUTPUT:
%       No outputs for this function, as it directly saves the variables
%       as a mat-file.
%
% AUTHOR:   Mathieu Dever 03-03-2016
%
% DEPENDENCIES:
%       - load_ctd, import CTD data
%       - load_fls, import backscatter data
%       - load_flsv2, import fluorescence data
%       - load_irrad, import irradiance data
%       - load_oxy_wphase, import oxygen data
%       - cal_oxy_wphase, calibrate oxygen data using temperature and
%       salinity from CTD
%       - TEOS-10 GSW package (http://www.teos-10.org/software.htm#1)
%
% REFERENCE:
%       - TEOS-10 GSW package (http://www.teos-10.org/software.htm#1)
%
% UPDATES:
%==================================================================

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------

% Check the number of inputs
if iseven(nargin)==1
    error(['glider_caster.m: wrong number of arguments. Optional argument',...
        'must be specified as pairs.'])
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
end; clear ii vin

% If optional argument is not specified, assumes 'no'
if nargin == 1
    plot_request = 'no';
end
clear varargin


%----------------------
% CORE CODE
%----------------------

%% Load each science component
% CTD
if exist([csv_root,'_water.csv'],'file')==2
    tic
    ctd = load_ctd(csv_root);
    disp('%%%%%%%%%%%%%%%%%%%%')
    disp(['CTD extracted in ',num2str(toc),' seconds !'])
    disp('%%%%%%%%%%%%%%%%%%%%')
else
    warning('No CTD file found!')
end

% fluorescence
if exist([csv_root,'_bb2fls.csv'],'file')==2
    tic
    fls = load_fls(csv_root);
    disp('%%%%%%%%%%%%%%%%%%%%')
    disp(['FLS extracted in ',num2str(toc),' seconds !'])
    disp('%%%%%%%%%%%%%%%%%%%%')
else
    warning('No Backscatter file found!')
end

% fluorescence 2
if exist([csv_root,'_bb2flsv2.csv'],'file')==2
    tic
    flsv2 = load_flsv2(csv_root);
    disp('%%%%%%%%%%%%%%%%%%%%')
    disp(['FLSv2 extracted in ',num2str(toc),' seconds !'])
    disp('%%%%%%%%%%%%%%%%%%%%')
else
    warning('No Fluorescence file found!')
end

% Light sensor
if exist([csv_root,'_ocr504i.csv'],'file')==2
    tic
    irrad = load_irrad(csv_root);
    disp('%%%%%%%%%%%%%%%%%%%%')
    disp(['IRRAD extracted in ',num2str(toc),' seconds !'])
    disp('%%%%%%%%%%%%%%%%%%%%')
else
    warning('No Irradiance file found!')
end

% Oxygen sensor
if exist([csv_root,'_oxy3835.csv'],'file')==2
    tic
    oxy = load_oxy_wphase(csv_root);
    disp('%%%%%%%%%%%%%%%%%%%%')
    disp(['OXY extracted in ',num2str(toc),' seconds !'])
    disp('%%%%%%%%%%%%%%%%%%%%')
    % re-compute O2 based on the ctd-measured S and T
    if exist('ctd','var')==1
        oxy = cal_oxy_wphase(ctd,oxy);
        disp('%%%%%%%%%%%%%%%%%%%%')
        disp(['OXY calibrated in ',num2str(toc),' seconds !'])
        disp('%%%%%%%%%%%%%%%%%%%%')
    else
        warning(['Oxygen could not be re-calibrated using temperature and',...
            'salinity fomr the CTD, as no CTD data could be found'])
    end
else
    warning('No Oxygen file found!')
end

%% PLOT the track

if strcmp(plot_request,'yes')==1
    figure
    if exist('ctd','var')==1
        plot(ctd.lon,ctd.lat,'.k');
    elseif exist('fls','var')==1
        plot(fls.lon,fls.lat,'.k');
    elseif exist('flsv2','var')==1
        plot(flsv2.lon,flsv2.lat,'.k');
    elseif exist('irrad','var')==1
        plot(irrad.lon,irrad.lat,'.k');
    elseif exist('oxy','var')==1
        plot(oxy.lon,oxy.lat,'.k');
    end
    grid on; xlabel('Longitude (°)'); ylabel('latitude(°)')
    set(gca,'fontsize',16);
end
 
clear plot_request

%% Save the mat file

C = who;
[FileName,PathName,~] = uiputfile('*.mat','Save mat-file',[csv_root,'.mat']);
save([PathName,FileName],C{:,1});clear C

disp('%%%%%%%%%%%%%%%%%%%%')
disp('Mat-file saved!')
disp('%%%%%%%%%%%%%%%%%%%%')