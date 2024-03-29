function [ vLapseRates ] = GCM_Surface_Lapse_Rates1( sGCM, sTempRes, kLat, kLong )
%% Extract surface temperature lapse rates from FLOR, HAR, or WRF



%% File locations and specifics for various GCM data

if strcmp(sGCM,'HAR')
    sFileDirectory = '/uufs/chpc.utah.edu/common/home/steenburgh-group6/Eric/HAR_Data/Temp';
elseif strcmp(sGCM,'WRF')
    sFileDirectory = '';
elseif strcmp(sGCM,'FLOR')
    sFileDirectory = '/uufs/chpc.utah.edu/common/home/steenburgh-group6/Eric/GFDL_FLOR_Data';
else
    error('Error: Incorrect GCM input')
end

% Add data directory and all subfolders to path
addpath(genpath(sFileDirectory))

%% Load temperature data

% HAR
if strcmp(sGCM,'HAR')
    
    mLat = double(ncread('har_d10km_h_2d_t2_2001.nc', 'lat'));
    mLong = double(ncread('har_d10km_h_2d_t2_2001.nc','lon'));

    vLat = mLat(:);
    vLong = mLong(:);
    [~, iLatLongIdx] = min((vLat - kLat).^2 + (vLong - kLong).^2);
    [m_Idx, n_Idx] = ind2sub([270, 180], iLatLongIdx);
    
    m3T_a_3x3 = [];

    if strcmp(sTempRes,'hourly')
        
        for i = 0:14

            if i < 10
                sFileName = ['har_d10km_h_2d_t2_200', num2str(i), '_rechunked.nc'];
            else
                sFileName = ['har_d10km_h_2d_t2_20', num2str(i), '_rechunked.nc'];
            end

            if rem(i,4) == 0
                Leap = 1;
            else
                Leap = 0;
            end
if m_Idx == 1
    m_Idx = 2;
    disp('changed m index to 2 in lapse rates script')
end
            vTemp_Data_3x3 = double(squeeze(ncread(sFileName, 't2',[m_Idx-1,n_Idx-1,1], [3,3,8760+24*Leap], [1,1,1])));
            m3T_a_3x3 = cat(3, m3T_a_3x3, vTemp_Data_3x3);

        end
     
        % Convert Kelvin to Celsius
        m3T_a_3x3 = m3T_a_3x3 - 273.15;
        
        % Remove leap day data
        m3T_a_3x3(:,:,[24*58+1:24*59,24*1518+1:24*1519,24*2978+1:24*2979,24*4438+1:24*4439]) = [];
        % Average hourly data for each day
        m3HA_Daily_Temp = zeros(size(m3T_a_3x3,1),size(m3T_a_3x3,2),size(m3T_a_3x3,3)/24);
        
        for t = 1:size(m3T_a_3x3,3)/24
            m3HA_Daily_Temp(:,:,t) = mean(m3T_a_3x3(:,:,(t-1)*24+1:t*24),3); 
        end
        
    elseif strcmp(sTempRes,'daily')
        
        for i = 0:14

            if i < 10
                sFileName = ['har_d10km_d_2d_t2_200', num2str(i), '_rechunked.nc'];
            else
                sFileName = ['har_d10km_d_2d_t2_20', num2str(i), '_rechunked.nc'];
            end

            if rem(i,4) == 0
                Leap = 1;
            else
                Leap = 0;
            end

            vTemp_Data_3x3 = double(squeeze(ncread(sFileName, 't2',[m_Idx-1,n_Idx-1,1], [3,3,365+Leap], [1,1,1])));
            m3T_a_3x3 = cat(3, m3T_a_3x3, vTemp_Data_3x3);

        end
        
        % Convert Kelvin to Celsius
        m3T_a_3x3 = m3T_a_3x3 - 273.15;
        
        % Remove leap day data
        m3T_a_3x3(:,:,[59,1519,2979,4439]) = [];
        % Average hourly data for each day
        m3HA_Daily_Temp = m3T_a_3x3;
        
    end
    
    % Average daily data for each year
    m3DA_Annual_Temp = zeros(size(m3T_a_3x3,1),size(m3T_a_3x3,2),365);
    
    for t = 1:365
    
        m3DA_Annual_Temp(:,:,t) = nanmean(m3HA_Daily_Temp(:,:,t:365:end),3);
    
    end
    
    % Get elevation data for HAR
    addpath('/uufs/chpc.utah.edu/common/home/steenburgh-group6/Eric/HAR_Data')
    mElevations =  ncread('/uufs/chpc.utah.edu/common/home/steenburgh-group6/Eric/HAR_Data/har_d10km_static_hgt.nc','hgt');
    mElevations_3x3 = mElevations(m_Idx-1:m_Idx+1, n_Idx-1:n_Idx+1);
    
% WRF

% FLOR
elseif strcmp(sGCM,'FLOR')
    
    addpath(genpath(sFileDirectory))
    iEnsNum = 8;

    vLat = double(ncread('Model_atmos_static_elev.nc', 'lat'));
    vLong = double(ncread('Model_atmos_static_elev.nc','lon'));
    mLat = repmat(flipud(vLat),1,length(vLong));
    mLong = repmat(vLong',length(vLat),1);

    vLat = mLat(:);
    vLong = mLong(:);

    [~, iLatLongIdx] = min((vLat - kLat).^2 + (vLong - kLong).^2);
    [m_Idx, n_Idx] = ind2sub([28, 55], iLatLongIdx);
    
    m3T_a1 = ncread('t_ref_ens_01-10_1961_2020.nc','t_ref',[n_Idx-1,28-m_Idx,1,1],[3,3,21915,10],[1,1,1,1]);
    m3T_a2 = ncread('t_ref_ens_11-20_1961_2020.nc','t_ref',[n_Idx-1,28-m_Idx,1,1],[3,3,21915,10],[1,1,1,1]);
    m3T_a3 = ncread('t_ref_ens_21-30_1961_2020.nc','t_ref',[n_Idx-1,28-m_Idx,1,1],[3,3,21915,10],[1,1,1,1]);
%     m3T_a1 = ncread('t_ref_ens_01-10_1961_2020.nc','t_ref',[n_Idx-1:n_Idx+1,(28-m_Idx):(28-m_Idx+2),1,1,1],[1,1,21915,10],[1,1,1,1]);
%     m3T_a2 = ncread('t_ref_ens_11-20_1961_2020.nc','t_ref',[n_Idx-1:n_Idx+1,(28-m_Idx):(28-m_Idx+2),1,1,1],[1,1,21915,10],[1,1,1,1]);
%     m3T_a3 = ncread('t_ref_ens_21-30_1961_2020.nc','t_ref',[n_Idx-1:n_Idx+1,(28-m_Idx):(28-m_Idx+2),1,1,1],[1,1,21915,10],[1,1,1,1]);
    m3T_a = cat(4,m3T_a1,m3T_a2,m3T_a3);
    m3T_a_3x3 = squeeze(m3T_a(:,:,:,iEnsNum)) - 273;
    
    % Remove leap day data
    m3T_a_3x3(:,:,[365*2+31+29,365*6+31+29+1,365*10+31+29+2,365*14+31+29+3,365*18+31+29+4,365*22+31+29+5,365*26+31+29+6,...
        365*30+31+29+7,365*34+31+29+8,365*38+31+29+9,365*42+31+29+10,365*46+31+29+11]) = [];
    
    % Average daily data for each year
    m3DA_Annual_Temp = zeros(size(m3T_a_3x3,1),size(m3T_a_3x3,2),365);
    
    for t = 1:365
    
        m3DA_Annual_Temp(:,:,t) = nanmean(m3T_a_3x3(:,:,t:365:end),3);
    
    end

    mElevations = double(ncread('Model_atmos_static_elev.nc', 'zsurf'));
    mElevations = flipud(mElevations');
    mElevations_3x3 = mElevations(m_Idx-1:m_Idx+1, n_Idx-1:n_Idx+1);
    
else
    error('Error: You havent added functionality for that GCM yet. Way to go pinhead.')
end

%% Fit regression line to temperature data

vLapseRates = zeros(365,1);

warning('off','all')

for t = 1:365

    % Extract temp and elevation data for gridcell + neighbors
    vAltNeighbors_temp = mElevations_3x3(:); %9x1
    mTempNeighbors_temp = squeeze(m3DA_Annual_Temp(:,:,t)); %3x3
    % Calculate lapse rate (deg C km^-1)
%     scatter(vAltNeighbors_temp,mTempNeighbors_temp(:),'b')
    kLapseRate = polyfit(vAltNeighbors_temp,mTempNeighbors_temp(:),1); 
    kLapseRate = kLapseRate(1) * 1000;
    vLapseRates(t) = kLapseRate;   
end
warning('on','all')
