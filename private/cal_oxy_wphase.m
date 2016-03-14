function [oxy] = cal_oxy_wphase(ctd,oxy)
%==========================================================================
% DESCRIPTION:
%       [oxy] = cal_oxy_wphase(csv_root)
%       Re-calculates the oxygen concentration and oxygen saturation
%       according for salinity/temperature from CTD
%
% INPUTS:
%       'csv_root' is the core filename (or path) of the csv mission file
%
% OUTPUTS:
%       'oxy' is a structure that contains:
%           - Latitude (°)
%           - Longitude (°)
%           - Depth (m)
%           - Time (matlab datenum format)
%           - Temperature fomr optode
%           - Oxygen concentration (µM)
%           - Re-calculated Oxygen concentration (µmol/kg)
%           - Oxygen saturation (%)
%           - Recalculated Oxygen saturation (%)
%           - Oxygen saturation level (µmol/kg)
%
% DEPENDENCIES:
%       the GSW TEOS-10 mfiles
%
% AUTHOR:   Mathieu Dever 05-11-2014
% UPDATES:
%==========================================================================


S0=35.;
B0=-6.24097e-3; B1=-6.93498e-3; B2=-6.90358e-3; B3=-4.29155e-3; C0=-3.1168e-7;

% Find the nearest (in time) temperature and salinity measurements from CTD
% for each of the optopde measures

if length(oxy.wphase_o2)==1 && isnan(oxy.wphase_o2)==1
    oxy.sat = NaN;
    oxy.o2percentsat = NaN;
    oxy.o2 = NaN;
    warning('Oxygen measurements are not available, NaN have been attributed to o2, o2sat and o2percentsat')
    
else
    for ii = 1:length(oxy.time)
        % Finds the closest (S,T)
        [~, ind] = min(abs(ctd.time-oxy.time(ii)));
        % records the averag S and T
        S(ii) = nanmean(ctd.s(ind)); 
        T(ii) = nanmean(ctd.tw(ind));
        P(ii) = nanmean(ctd.press(ind));
    end
    
    % Re-compute the oxygen concentration
    Ts = log((298.15-T)./(273.15+T)); % ??
    oxy.o2(:,1) = oxy.wphase_o2'.*exp((S-S0).*(B0+B1.*Ts+B2.*Ts.^2+B3.*Ts.^3)+C0.*(S.^2-S0.^2));
    
    % Calculates the saturation level
    SA = gsw_SA_from_SP(S,P,oxy.lon,oxy.lat); % Absolute Salinity
    CT = gsw_CT_from_t(SA,T,P); % Conservative Temperature
    oxy.sat(:,1) = gsw_O2sol(SA,CT,P,oxy.lon,oxy.lat);
    %o2sat=sw_satO2(S,T)*44.6589; % factor to convert ml/l to uM
    
    oxy.o2percentsat = oxy.o2./oxy.sat.*100;
end


