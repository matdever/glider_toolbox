function [PA,PB,PC] = glider_maxwells_triangle(TA,SA,TB,SB,TC,SC,temperature,salinity,varargin)
%
% Project maxwell's triangle onto a mixing triangle to assign a
% 2-dimensional colorbar. Returns the % of each endmember for a given
% datapoint.
%==========================================================================
%
% USAGE:
%   [PA,PB,PC] = glider_maxwells_triangle(TA,SA,TB,SB,TC,SC,temperature,salinity)
%
% DESCRIPTION:
%   This function plots the T-S distribution, with a unique color
%   associated to each datapoint, obtained by projecting maxwell's triangle
%   onto the mixing triangle. It returns the percentage of each water mass
%   forming a given datapoint.
%
% INPUTS:
%   (TA,SA,TB,SB,TC,SC)     = Temperature and Salinity characterizing each
%                           endmember.
%   temperature             = temperature distribution
%   salinity                = salinity distribution
%
% OPTIONAL INPUTS (varargin):
%
%       'colorbar' - plots an additional figure showing the 2-dimensional
%       colorbar projected onto the mixing triangle defined by
%       (TA,SA,TB,SB,TC,SC)
%
%
% OUTPUTS:
%       PA,PB,PC    = Percentage of watermass A,B and C for a given
%                       datapoint. Has the same dimensions as the inputs
%                       (temperature,salinity)
%
% AUTHOR:    Mathieu Dever 20-04-2016
%
% DEPENDENCIES:
%       - gsw_rho.m: GIBBS SEAWATER PACKAGE
%
% REFERENCE:
%       Judd, D. B. (1935). A maxwell triangle yielding uniform
%           chromaticity scales. Journal of the Optical Society of America,
%           25(1), 24?35.
%
% UPDATES:
%=========================================================================

%% ----------------------
% CHECK INPUT ARGUMENTS
%------------------------

if nargin < 8
    error('myApp:argChk', 'Minimum 8 inputs required: please refer to the function help.')
end
if nargout ~= 3
    error('myApp:argChk', 'Wrong number of outputs. Please refer to the function help.')
end

vin = varargin;
if length(vin)>1
    error('Only 1 possible optional input (for now!). Please refer to the function help.')
end

% No colorbar plot by default
colorbar_plot = 0;
for ii = 1:length(vin)
    % Request the 2-dimensional colorbar plot
    if strcmp(vin{ii},'colorbar')==1
        colorbar_plot = 1;
    end
end

%% ----------------------
% CORE CODE
%------------------------

%% 2-dimensional colorbar plot - if requested 
if colorbar_plot==1
    
    % Plot the triangular colorbar (Maxwell's triangle)
    % Create base vectors
    grid.s = 30:0.01:36;
    grid.t = 0:0.01:12;
    [grid.s,grid.t] = meshgrid(grid.s,grid.t);
    in = inpolygon(grid.s,grid.t,[SA SB SC],[TA TB TC]);
    
    SX = grid.s(:);
    TX = grid.t(:);
    
    % Distance between mixing triangle's vertices
    AB = sqrt((SB-SA).^2+(TB-TA).^2);
    AC = sqrt((SA-SC).^2+(TA-TC).^2);
    BC = sqrt((SB-SC).^2+(TB-TC).^2);
    
    % Distance between gridpoints and mixing triangle's vertices
    AX = sqrt((SX-SA).^2+(TX-TA).^2);
    BX = sqrt((SB-SX).^2+(TB-TX).^2);
    CX = sqrt((SC-SX).^2+(TC-TX).^2);
    
    % Amount of C
    % distance between C and AB
    h0 = sqrt(BC.^2 - ((AB.^2 - AC.^2 + BC.^2)/(2*AB)).^2);
    % distance between X and AB
    h = sqrt(BX.^2 - ((AB.^2 - AX.^2 + BX.^2)/(2*AB)).^2);
    C = h/h0;
    clear h h0
    
    % Amount of A
    % distance between A and BC
    h0 = sqrt(AC.^2 - ((BC.^2 - AB.^2 + AC.^2)/(2*BC)).^2);
    % distance between X and BC
    h = sqrt(CX.^2 - ((BC.^2 - BX.^2 + CX.^2)/(2*BC)).^2);
    A = h/h0;
    clear h h0
    
    % Amount of B
    % distance between B and AC
    h0 = sqrt(AB.^2 - ((AC.^2 - BC.^2 + AB.^2)/(2*AC)).^2);
    % distance between X and AC
    h = sqrt(AX.^2 - ((AC.^2 - CX.^2 + AX.^2)/(2*AC)).^2);
    B = h/h0;
    clear h h0
    
    
    %% Create the plot
    
    % Plot the 2-dimensional colorbar
    figure
    % Isopycnals in the background
    [XS,XT] = meshgrid(29:37,-5:20);
    RHO = gsw_rho(XS,XT,10.1325);
    [CC,h]=contour(XS,XT,RHO,'color',[0.7 0.7 0.7]);
    clabel(CC,h,'labelspacing',3000,'fontsize',50)
    hold on
    
    % Plot the colored gridpoints
    scatter(grid.s(in),grid.t(in),35,...
        [A(in)./mean([A(in) B(in) C(in)],2) B(in)./mean([A(in) B(in) C(in)],2) C(in)./mean([A(in) B(in) C(in)],2)],...
        'filled')
    hold on
    
    % Identify the endmembers
    scatter([SA SB SC],[TA TB TC],35,'ok','filled')
    plot([SA SB SC SA],[TA TB TC TA],'--','color','k')
    text(SA-.3,TA,'A','fontsize',16)
    text(SB+.3,TB,'B','fontsize',16)
    text(SC,TC-.8,'C','fontsize',16)
    
    % misc
    axis([29 37 -5 20])
    set(gca,'fontsize',24)
    ylabel('Temperature (°C)','fontsize',24)
    xlabel('Salinity','fontsize',24)
    axis square
    
    % Clean-up
    clearvars -except TA TB TC SA SB SC temperature salinity
    
end

%% Figure 2 - TS distribution with 2D colorbar

temperature2 = temperature(:);
salinity2 = salinity(:);

% Identify datapoints outside of the mixing triangle
IN = inpolygon(salinity,temperature,[SA SB SC],[TA TB TC]);

% Projects datapoints outside of the mixing triangle onto the mixing
% line for color purposes

% Project on line AB datapoints above the triangle
ind2 = find(IN == 0 & (interp1([SA SB],[TA TB],salinity,'Linear','extrap')-temperature)<0);
temperature2(ind2) = (TB - TA)/(SB - SA)*salinity(ind2) + (TA - (TB - TA)/(SB - SA)*SA);
clear ind2
% Project on line AC datapoints below and left of the triangle
ind2 = find(IN == 0 & (interp1([SA SC],[TA TC],salinity,'Linear','extrap')-reshape(temperature2,size(temperature,1),size(temperature,2)))>0);
temperature2(ind2) = (TC - TA)/(SC - SA)*salinity(ind2) + (TA - (TC - TA)/(SC - SA)*SA);
clear ind2
% Project on line BC datapoints below and right of the triangle
ind2 = find(IN == 0 & (interp1([SB SC],[TB TC],salinity,'Linear','extrap')-temperature)>0);
temperature2(ind2) = (TC - TB)/(SC - SB)*salinity(ind2) + (TB - (TC - TB)/(SC - SB)*SB);
clear ind2
ind2 = find(salinity2<SA);
temperature2(ind2) = TA;
salinity2(ind2) = SA;

temperature2 = reshape(temperature2,size(temperature,1),size(temperature,2));
salinity2 = reshape(salinity2,size(salinity,1),size(salinity,2));

% Calculate the percentage of all three sort of water masses
PC = (temperature2*SB - temperature2*SA - TA*(SB - salinity2) -...
    TB*(salinity2-SA))/(TA*(SC-SB) + TB*(SA-SC) + TC*(SB-SA));
PB = (salinity2 - SA - PC*(SC-SA))/(SB-SA);
PA = 1 - PB - PC;
PC(PC<0) = 0;
PB(PB<0) = 0;
PA(PA<0) = 0;
PA = PA*100; PB = PB*100; PC = PC*100;

% Computes color associated with measurements
SX = salinity2(:);
TX = temperature2(:);
AB = sqrt((SB-SA).^2+(TB-TA).^2);
AC = sqrt((SA-SC).^2+(TA-TC).^2);
BC = sqrt((SB-SC).^2+(TB-TC).^2);
AX = sqrt((SX-SA).^2+(TX-TA).^2);
BX = sqrt((SB-SX).^2+(TB-TX).^2);
CX = sqrt((SC-SX).^2+(TC-TX).^2);

% Amount of C
% distance between C and AB
h0 = sqrt(BC.^2 - ((AB.^2 - AC.^2 + BC.^2)/(2*AB)).^2);
% distance between X and AB
h = sqrt(BX.^2 - ((AB.^2 - AX.^2 + BX.^2)/(2*AB)).^2);
C = abs(h)/h0;
clear h h0

% Amount of A
% distance between A and BC
h0 = sqrt(AC.^2 - ((BC.^2 - AB.^2 + AC.^2)/(2*BC)).^2);
% distance between X and BC
h = sqrt(CX.^2 - ((BC.^2 - BX.^2 + CX.^2)/(2*BC)).^2);
A = abs(h)/h0;
clear h h0

% Amount of B
% distance between B and AC
h0 = sqrt(AB.^2 - ((AC.^2 - BC.^2 + AB.^2)/(2*AC)).^2);
% distance between X and AC
h = sqrt(AX.^2 - ((AC.^2 - CX.^2 + AX.^2)/(2*AC)).^2);
B = abs(h)/h0;
clear h h0

% Plot the T-S distribution
figure
% Isopycnals in the background
[XS,XT] = meshgrid(29:37,-5:20);
RHO = gsw_rho(XS,XT,10.1325);
[CC,h] = contour(XS,XT,RHO,'color',[0.7 0.7 0.7]);
clabel(CC,h,'labelspacing',3000,'fontsize',50)
hold on

% Plot the T-S distribution
scatter(salinity(:),temperature(:),35,[A./mean([A B C],2) B./mean([A B C],2) C./mean([A B C],2)],'filled')
hold on

% Identify the endmembers
scatter([SA SB SC],[TA TB TC],35,'ok','filled')
plot([SA SB SC SA],[TA TB TC TA],'--','color','k')
text(SA-.3,TA,'A','fontsize',16)
text(SB+.3,TB,'B','fontsize',16)
text(SC,TC-.8,'C','fontsize',16)

% misc
axis([29 37 -5 20])
set(gca,'fontsize',24)
ylabel('Temperature (°C)','fontsize',24)
xlabel('Salinity','fontsize',24)
axis square