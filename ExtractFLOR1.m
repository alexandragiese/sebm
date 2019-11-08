function [ vLW_in, vPrecip, vPressure, vSpecificHumidity, vRH, vT_a, vSW_in, vWindSpeed, stTime, Data] = ExtractFLOR1(kLat, kLong, sFileDir, iEnsNum) 

%% Extract FLOR data

% addpath(genpath('/uufs/chpc.utah.edu/common/home/steenburgh-group6/Eric/GFDL_FLOR_Data'))
addpath(genpath(sFileDir))
% cd(sFileDir)


vLat = double(ncread('Model_atmos_static_elev.nc', 'lat'));
% mLat = flipud(mLat');
vLong = double(ncread('Model_atmos_static_elev.nc','lon'));
% mLong = flipud(mLong');

mLat = repmat(flipud(vLat),1,length(vLong));
mLong = repmat(vLong',length(vLat),1);

vLat = mLat(:);
vLong = mLong(:);

[~, iLatLongIdx] = min((vLat - kLat).^2 + (vLong - kLong).^2);
[m_Idx, n_Idx] = ind2sub([28, 55], iLatLongIdx);

%% Grid Elevation

mAWS_Alt = double(ncread('Model_atmos_static_elev.nc', 'zsurf'));
mAWS_Alt = flipud(mAWS_Alt');
Data.AWS_Alt = mAWS_Alt(m_Idx, n_Idx);

%% Extract Variables

% iEnsNum = 1;

m3LW_in1 = ncread('lwdn_sfc_ens_01-10_1961_2020.nc','lwdn_sfc',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3LW_in2 = ncread('lwdn_sfc_ens_11-20_1961_2020.nc','lwdn_sfc',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3LW_in3 = ncread('lwdn_sfc_ens_21-30_1961_2020.nc','lwdn_sfc',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3LW_in = cat(4,m3LW_in1,m3LW_in2,m3LW_in3);
% vLW_in = squeeze(nanmean(m3LW_in,4));
vLW_in = squeeze(m3LW_in(:,:,:,iEnsNum));

% mm s^-1 to m day^-1
m3Precip1 = ncread('precip_ens_01-10_1961_2020.nc','precip',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3Precip2 = ncread('precip_ens_11-20_1961_2020.nc','precip',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3Precip3 = ncread('precip_ens_21-30_1961_2020.nc','precip',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3Precip = cat(4,m3Precip1,m3Precip2,m3Precip3);
% vPrecip = squeeze(nanmean(m3Precip,4))*24*60*60/1000;
vPrecip = squeeze(m3Precip(:,:,:,iEnsNum))*24*60*60/1000;

% hPa
m3Pressure1 = ncread('ps_ens_01-10_1961_2020.nc','ps',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3Pressure2 = ncread('ps_ens_11-20_1961_2020.nc','ps',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3Pressure3 = ncread('ps_ens_21-30_1961_2020.nc','ps',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3Pressure = cat(4,m3Pressure1,m3Pressure2,m3Pressure3);
% vPressure = squeeze(nanmean(m3Pressure,4))/100;
vPressure = squeeze(m3Pressure(:,:,:,iEnsNum))/100;

m3SpecificHumidity1 = ncread('q_ref_ens_01-10_1961_2020.nc','q_ref',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3SpecificHumidity2 = ncread('q_ref_ens_11-20_1961_2020.nc','q_ref',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3SpecificHumidity3 = ncread('q_ref_ens_21-30_1961_2020.nc','q_ref',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3SpecificHumidity = cat(4,m3SpecificHumidity1,m3SpecificHumidity2,m3SpecificHumidity3);
% vSpecificHumidity = squeeze(nanmean(m3SpecificHumidity,4));
vSpecificHumidity = squeeze(m3SpecificHumidity(:,:,:,iEnsNum));

% %
m3RH1 = ncread('rh_ref_ens_01-10_1961_2020.nc','rh_ref',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3RH2 = ncread('rh_ref_ens_11-20_1961_2020.nc','rh_ref',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3RH3 = ncread('rh_ref_ens_21-30_1961_2020.nc','rh_ref',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3RH = cat(4,m3RH1,m3RH2,m3RH3);
% vRH = squeeze(nanmean(m3RH,4));
vRH = squeeze(m3RH(:,:,:,iEnsNum));

% deg C
m3T_a1 = ncread('t_ref_ens_01-10_1961_2020.nc','t_ref',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3T_a2 = ncread('t_ref_ens_11-20_1961_2020.nc','t_ref',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3T_a3 = ncread('t_ref_ens_21-30_1961_2020.nc','t_ref',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3T_a = cat(4,m3T_a1,m3T_a2,m3T_a3);
% vT_a = squeeze(nanmean(m3T_a,4)) - 273;
vT_a = squeeze(m3T_a(:,:,:,iEnsNum)) - 273;

m3SW_in1 = ncread('swdn_sfc_ens_01-10_1961_2020.nc','swdn_sfc',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3SW_in2 = ncread('swdn_sfc_ens_11-20_1961_2020.nc','swdn_sfc',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3SW_in3 = ncread('swdn_sfc_ens_21-30_1961_2020.nc','swdn_sfc',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3SW_in = cat(4, m3SW_in1, m3SW_in2, m3SW_in3);
% vSW_in = squeeze(nanmean(m3SW_in,4));
vSW_in = squeeze(m3SW_in(:,:,:,iEnsNum));

m3WindSpeed1 = ncread('wind_ens_01-10_1961_2020.nc','wind',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3WindSpeed2 = ncread('wind_ens_11-20_1961_2020.nc','wind',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3WindSpeed3 = ncread('wind_ens_21-30_1961_2020.nc','wind',[n_Idx,28-m_Idx+1,1,1],[1,1,21915,10],[1,1,1,1]);
m3WindSpeed = cat(4,m3WindSpeed1,m3WindSpeed2,m3WindSpeed3);
% vWindSpeed = squeeze(nanmean(m3WindSpeed,4));
vWindSpeed = squeeze(m3WindSpeed(:,:,:,iEnsNum));

%% Other Data

% Temporal resolution of the raw data (minutes)
Data.Resolution_mins = 60*24;
Data.AWS_Lat = kLat;
Data.AWS_Long = kLong;
Data.TimeZone = 0;
Data.FirstDay = 0;

%% Calculate time vectors

stTime.year = zeros(length(vT_a),1);
stTime.month = zeros(length(vT_a),1);
stTime.day = zeros(length(vT_a),1);

kYearLength = 365;

for iYear = 1961:2020

    if rem(iYear,4) > 0
        iLeapYear = 0;
    elseif rem(iYear,4) == 0
        iLeapYear = 1;
    end
    
    iYearLength = kYearLength + iLeapYear;
    iLast_idx = find(stTime.year,1,'last');
    
    if isempty(iLast_idx)
        iLast_idx = 0;
    end
    
    stTime.year(iLast_idx+1:iLast_idx+iYearLength) = iYear;
    vMonth_temp = [repmat(1,[31,1]);repmat(2,[28+iLeapYear,1]);repmat(3,[31,1]);repmat(4,[30,1]);...
        repmat(5,[31,1]);repmat(6,[30,1]);repmat(7,[31,1]);repmat(8,[31,1]);...
        repmat(9,[30,1]);repmat(10,[31,1]);repmat(11,[30,1]);repmat(12,[31,1])];
    stTime.month(iLast_idx+1:iLast_idx+iYearLength) = vMonth_temp;
    vDay_temp = [1:31, 1:28+iLeapYear, 1:31, 1:30, 1:31, 1:30, 1:31, 1:31, 1:30, 1:31, 1:30, 1:31]';
    stTime.day(iLast_idx+1:iLast_idx+iYearLength) = vDay_temp;

end

stTime.hour = zeros(size(stTime.year));
stTime.min = zeros(size(stTime.year));
stTime.sec = zeros(size(stTime.year));
stTime.UTC = 5;


















end















































% m3LW_in1 = ncread('lwdn_sfc_ens_01-10_1961_2020.nc','lwdn_sfc');
% m3LW_in2 = ncread('lwdn_sfc_ens_11-20_1961_2020.nc','lwdn_sfc');
% m3LW_in3 = ncread('lwdn_sfc_ens_21-30_1961_2020.nc','lwdn_sfc');
% m3LW_in = cat(4,m3LW_in1,m3LW_in2,m3LW_in3);
% m3LW_in = nanmean(m3LW_in,4);
% m3LW_in = permute(m3LW_in,[2,1,3]);
% m3LW_in = m3LW_in(fliplr(1:size(m3LW_in,1)),:,:);
% vLW_in = squeeze(m3LW_in(m_Idx, n_Idx, :));
% 
% m3Precip1 = ncread('precip_ens_01-10_1961_2020.nc','precip');
% m3Precip2 = ncread('precip_ens_11-20_1961_2020.nc','precip');
% m3Precip3 = ncread('precip_ens_21-30_1961_2020.nc','precip');
% m3Precip = cat(4,m3Precip1,m3Precip2,m3Precip3);
% m3Precip = nanmean(m3Precip,4);
% m3Precip = permute(m3Precip,[2,1,3]);
% m3Precip = m3Precip(fliplr(1:size(m3Precip,1)),:,:);
% vPrecip = squeeze(m3Precip(m_Idx, n_Idx, :));
% 
% m3Pressure1 = ncread('ps_ens_01-10_1961_2020.nc','ps');
% m3Pressure2 = ncread('ps_ens_11-20_1961_2020.nc','ps');
% m3Pressure3 = ncread('ps_ens_21-30_1961_2020.nc','ps');
% m3Pressure = cat(4,m3Pressure1,m3Pressure2,m3Pressure3);
% m3Pressure = nanmean(m3Pressure,4);
% m3Pressure = permute(m3Pressure,[2,1,3]);
% m3Pressure = m3Pressure(fliplr(1:size(m3Pressure,1)),:,:);
% vPressure = squeeze(m3Pressure(m_Idx, n_Idx, :));
% 
% m3SpecificHumidity1 = ncread('q_ref_ens_01-10_1961_2020.nc','q_ref');
% m3SpecificHumidity2 = ncread('q_ref_ens_11-20_1961_2020.nc','q_ref');
% m3SpecificHumidity3 = ncread('q_ref_ens_21-30_1961_2020.nc','q_ref');
% m3SpecificHumidity = cat(4,m3SpecificHumidity1,m3SpecificHumidity2,m3SpecificHumidity3);
% m3SpecificHumidity = nanmean(m3SpecificHumidity,4);
% m3SpecificHumidity = permute(m3SpecificHumidity,[2,1,3]);
% m3SpecificHumidity = m3SpecificHumidity(fliplr(1:size(m3SpecificHumidity,1)),:,:);
% vSpecificHumidity = squeeze(m3SpecificHumidity(m_Idx, n_Idx, :));
% 
% m3RH1 = ncread('rh_ref_ens_01-10_1961_2020.nc','rh_ref');
% m3RH2 = ncread('rh_ref_ens_11-20_1961_2020.nc','rh_ref');
% m3RH3 = ncread('rh_ref_ens_21-30_1961_2020.nc','rh_ref');
% m3RH = cat(4,m3RH1,m3RH2,m3RH3);
% m3RH = nanmean(m3RH,4);
% m3RH = permute(m3RH,[2,1,3]);
% m3RH = m3RH(fliplr(1:size(m3RH,1)),:,:);
% vRH = squeeze(m3RH(m_Idx, n_Idx, :));
% 
% m3T_a1 = ncread('t_ref_ens_01-10_1961_2020.nc','t_ref');
% m3T_a2 = ncread('t_ref_ens_11-20_1961_2020.nc','t_ref');
% m3T_a3 = ncread('t_ref_ens_21-30_1961_2020.nc','t_ref');
% m3T_a = cat(4,m3T_a1,m3T_a2,m3T_a3);
% m3T_a = nanmean(m3T_a,4);
% m3T_a = permute(m3T_a,[2,1,3]);
% m3T_a = m3T_a(fliplr(1:size(m3T_a,1)),:,:);
% vT_a = squeeze(m3T_a(m_Idx, n_Idx, :));
% 
% m3SW_in1 = ncread('swdn_sfc_ens_01-10_1961_2020.nc','swdn_sfc');
% m3SW_in2 = ncread('swdn_sfc_ens_11-20_1961_2020.nc','swdn_sfc');
% m3SW_in3 = ncread('swdn_sfc_ens_21-30_1961_2020.nc','swdn_sfc');
% m3SW_in = cat(4, m3SW_in1, m3SW_in2, m3SW_in3);
% m3SW_in = nanmean(m3SW_in,4);
% m3SW_in = permute(m3SW_in,[2,1,3]);
% m3SW_in = m3SW_in(fliplr(1:size(m3SW_in,1)),:,:);
% vSW_in = squeeze(m3SW_in(m_Idx, n_Idx, :));
% 
% m3WindSpeed1 = ncread('wind_ens_01-10_1961_2020.nc','wind');
% m3WindSpeed2 = ncread('wind_ens_11-20_1961_2020.nc','wind');
% m3WindSpeed3 = ncread('wind_ens_21-30_1961_2020.nc','wind');
% m3WindSpeed = cat(4,m3WindSpeed1,m3WindSpeed2,m3WindSpeed3);
% m3WindSpeed = nanmean(m3WindSpeed,4);
% m3WindSpeed = permute(m3WindSpeed,[2,1,3]);
% m3WindSpeed = m3WindSpeed(fliplr(1:size(m3WindSpeed,1)),:,:);
% vWindSpeed = squeeze(m3WindSpeed(m_Idx, n_Idx, :));
% 


% %% Calculate time vectors
% 
% stTime.year = zeros(length(vT_a),1);
% stTime.month = zeros(length(vT_a),1);
% stTime.day = zeros(length(vT_a),1);
% 
% kYearLength = 8760;
% 
% for iYear = 1961:2020
% 
%     if rem(iYear,4) > 0
%         iLeapYear = 0;
%     elseif rem(iYear,4) == 0
%         iLeapYear = 1;
%     end
%     
%     iYearLength = kYearLength + iLeapYear*24;
%     iLast_idx = find(stTime.year,1,'last');
%     
%     if isempty(iLast_idx)
%         iLast_idx = 0;
%     end
%     
%     stTime.year(iLast_idx+1:iLast_idx+iYearLength) = iYear;
%     vMonth_temp = [repmat(1,[31,1]);repmat(2,[28+iLeapYear,1]);repmat(3,[31,1]);repmat(4,[30,1]);...
%         repmat(5,[31,1]);repmat(6,[30,1]);repmat(7,[31,1]);repmat(8,[31,1]);...
%         repmat(9,[30,1]);repmat(10,[31,1]);repmat(11,[30,1]);repmat(12,[31,1])];
%     stTime.month(iLast_idx+1:iLast_idx+iYearLength) = repelem(vMonth_temp,24);
%     vDay_temp = [1:31, 1:28+iLeapYear, 1:31, 1:30, 1:31, 1:30, 1:31, 1:31, 1:30, 1:31, 1:30, 1:31]';
%     vDay_temp = repelem(vDay_temp,24);
%     stTime.day(iLast_idx+1:iLast_idx+iYearLength) = vDay_temp;
% 
% end
% 
% stTime.hour = repelem((1:24)',length(stTime.year)/24);
% stTime.min = zeros(size(stTime.year));
% stTime.sec = zeros(size(stTime.year));
% stTime.UTC = 5;


































