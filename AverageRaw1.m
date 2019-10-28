function [ vTime, vS_in, vU, vT_a, vP_a, vRH, vP_d, vVP_a, vCloudiness] = AverageRaw1( kTS_Days, Data, Humidity, Precip, Pressure, SolarRad, Temp, Time, WindSpeed, Cloudiness)
% AverageAWSData:
% Average every n elements in each vector of raw data input

% Number of timesteps in one day
kN_per_Day = 60 * 24 / Data.Resolution;
% Number of timesteps to average
kN_Ave = kN_per_Day * kTS_Days;

% Time vector
vTime = arrayfun(@(i) nanmean(Time(i:i + kN_Ave - 1)), ...
    1:kN_Ave:length(Time) - kN_Ave + 1)';                                                               
% Incoming shortwave radiation (W m^2)
vS_in = arrayfun(@(i) nanmean(SolarRad(i:i + kN_Ave - 1)), ...
    1:kN_Ave:length(SolarRad) - kN_Ave + 1)';                                                            
% Wind speed at 2m above the surface (m/s)
vU = arrayfun(@(i) nanmean(WindSpeed(i:i + kN_Ave - 1)), ...
    1:kN_Ave:length(WindSpeed) - kN_Ave + 1)';                                                            
% Temperature of the atmosphere (C)
vT_a = arrayfun(@(i) nanmean(Temp(i:i + kN_Ave - 1)), ...
    1:kN_Ave:length(Temp) - kN_Ave + 1)';                                                        
% Actual air pressure (hPa)
vP_a = arrayfun(@(i) nanmean(Pressure(i:i + kN_Ave - 1)), ...
    1:kN_Ave:length(Pressure) - kN_Ave + 1)';                                                  
% Relative humidity (%)
vRH = arrayfun(@(i) nanmean(Humidity(i:i + kN_Ave - 1)), ...
    1:kN_Ave:length(Humidity) - kN_Ave + 1)';                                                         
% Precipitation (depth) rate (mm TS^-1)
vP_d = arrayfun(@(i) nansum(Precip(i:i + kN_Ave - 1)), ...
    1:kN_Ave:length(Precip) - kN_Ave + 1)';                                                                            
% Cloudiness (unitless fraction)
vCloudiness = arrayfun(@(i) nanmean(Cloudiness(i:i + kN_Ave - 1)), ...
    1:kN_Ave:length(Cloudiness) - kN_Ave + 1)';
% % Incoming longwave radiation (W m^-2)
% vLW_in_air = arrayfun(@(i) nanmean(Longwave_in_air(i:i + kN_Ave - 1)), ...
%     1:kN_Ave:length(Longwave_in_air) - kN_Ave + 1)';
% % Measured incoming longwave radiation (W m^-2)
% vLW_in_Meas = arrayfun(@(i) nanmean(LongwaveIn(i:i + kN_Ave - 1)), ...
%     1:kN_Ave:length(LongwaveIn) - kN_Ave + 1)'; 
% % Measured outgoing longwave radiation (W m^-2)
% vLW_out_Meas = arrayfun(@(i) nanmean(LongwaveOut(i:i + kN_Ave - 1)), ...
%     1:kN_Ave:length(LongwaveOut) - kN_Ave + 1)'; 
% % Measured net radiation (W m^-2)
% vNetRad_Meas = arrayfun(@(i) nanmean(NetRadiation(i:i + kN_Ave - 1)), ...
%     1:kN_Ave:length(NetRadiation) - kN_Ave + 1)'; 
% % Measured outgoing shortwave radiation (W m^-2)
% vSW_out_Meas = arrayfun(@(i) nanmean(ShortwaveOut(i:i + kN_Ave - 1)), ...
%     1:kN_Ave:length(ShortwaveOut) - kN_Ave + 1)'; 
% % Snow depth (cm)
% vSnowDepth_Meas = arrayfun(@(i) nanmean(SnowDepth(i:i + kN_Ave - 1)), ...
%     1:kN_Ave:length(SnowDepth) - kN_Ave + 1)'; 
% % Specific humidity (hPa)
% vSpecHum = arrayfun(@(i) nanmean(SpecificHumidity(i:i + kN_Ave - 1)), ...
%     1:kN_Ave:length(SpecificHumidity) - kN_Ave + 1)'; 


vVP_a = (6.11 * exp((2.5e6 / 461) * (1 / 273 - 1 ./ (273 + vT_a)))) ...
    .* vRH / 100;

























end

