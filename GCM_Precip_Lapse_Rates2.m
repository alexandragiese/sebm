function [ vPrcpLapseRates, mPolyfit_prcp ] = GCM_Precip_Lapse_Rates2( sGCM, sTempRes, kLat, kLong, stTime )
%% Extract precip lapse rates from HAR

%% File locations and specifics for various GCM data

if strcmp(sGCM,'HAR')
    sFileDirectory = '/uufs/chpc.utah.edu/common/home/steenburgh-group6/Eric/HAR_Data/Temp';
else
    error('Error: Incorrect GCM input')
end

% Add data directory and all subfolders to path
addpath(genpath(sFileDirectory))

%% Load precipitation data

% HAR
if strcmp(sGCM,'HAR')
    
    mLat = double(ncread('har_d10km_h_2d_t2_2001.nc', 'lat'));
    mLong = double(ncread('har_d10km_h_2d_t2_2001.nc','lon'));

    vLat = mLat(:);
    vLong = mLong(:);
% kLat = 40.4885;
% kLong = 101.7554;
    [~, iLatLongIdx] = min((vLat - kLat).^2 + (vLong - kLong).^2); %suppress actual min, get index in vectors vLat,vLong
    [m_Idx, n_Idx] = ind2sub([270, 180], iLatLongIdx); %270, 180 is hardcoded for size(mLat) = size(mLong)
    %m_Idx: 1 - 270 longitude
    %n_Idx: 1 - 180 latitude
    
    m3P_3x3 = [];
    m3YN_3x3 = [];

    if strcmp(sTempRes,'hourly')
        
        for i = 0:14

        if i < 10
            sFileName = ['har_d10km_h_2d_prcp_200', num2str(i), '_rechunked.nc'];
        else
            sFileName = ['har_d10km_h_2d_prcp_20', num2str(i), '_rechunked.nc'];
        end
        
            if rem(i,4) == 0
                Leap = 1;
            else
                Leap = 0;
            end
            
% %problem line below!!!
% if (180-m_Idx) >= 178
%     m_Idx = 2; % If 178 doesn't work, try 177 in the line above, and 3 in this line
% end

            vPrecip_Data_3x3 = double(squeeze(ncread(sFileName, 'prcp',[m_Idx-1,n_Idx-1,1], [3,3,8760+24*Leap], [1,1,1]))); %3x3 all hr in yr (full is 270x180x8784)
% keyboard
            vPrecip_YN_3x3 = zeros(size(vPrecip_Data_3x3));
                vPrecip_YN_3x3(1,:,:) = spones(squeeze(vPrecip_Data_3x3(1,:,:)));
                vPrecip_YN_3x3(2,:,:) = spones(squeeze(vPrecip_Data_3x3(2,:,:)));
                vPrecip_YN_3x3(3,:,:) = spones(squeeze(vPrecip_Data_3x3(3,:,:)));

            m3P_3x3 = cat(3, m3P_3x3, vPrecip_Data_3x3); %gridded precip data around point of interest (3x3 in all hr for whole stTime)
            m3YN_3x3 = cat(3, m3YN_3x3, vPrecip_YN_3x3);


        end
 
    % Sum Precipitation (monthly total & monthly count of +P hours) for each month
    m3MT_Precip    = nan(3,3,12); 
    m3MC_posPrecip = nan(3,3,12); 


    % fill 3D-matrices <-- checked one by hand
    t=1;
    for y = 2001:2013
        for m = 1:12
            m3MT_Precip(:,:,t) = nansum(m3P_3x3(:,:,stTime.month==m & stTime.year==y),3); %mm each month over all time / # yr
            m3MC_posPrecip(:,:,t) = nansum(m3YN_3x3(:,:,stTime.month==m & stTime.year==y),3);
            t = t+1;
        end
    end

%     elseif strcmp(sTempRes,'daily') %CURRENTLY NO PRECIP LAPSE RATE CALC FOR DT = 1DAY
%         
%         for i = 0:14
% 
%             if i < 10
%                 sFileName = ['har_d10km_d_2d_prcp_200', num2str(i), '_rechunked.nc'];
%             else
%                 sFileName = ['har_d10km_d_2d_prcp_20', num2str(i), '_rechunked.nc'];
%             end
% 
%             if rem(i,4) == 0
%                 Leap = 1;
%             else
%                 Leap = 0;
%             end
% 
%             vPrcp_Data_3x3 = double(squeeze(ncread(sFileName, 'prcp',[n_Idx-1,180-m_Idx,1], [3,3,365+Leap], [1,1,1])));
%             m3P_3x3 = cat(3, m3P_3x3, vTemp_Data_3x3);
% 
%         end
%         
%         % Remove leap day data
%         m3P_3x3(:,:,[59,1519,2979,4439]) = [];
%         % Average hourly data for each day
%         m3HA_Daily_Prcp = m3P_3x3;       
    end

%Calculate monthly data, averaged over all 13 years
mPrcpNeighbors_ini  = m3MT_Precip./m3MC_posPrecip; %3x3x12*13
mPrcpNeighbors_ini2 = zeros(size(m3P_3x3,1),size(m3P_3x3,2),12);
    
for m = 1:12
    mPrcpNeighbors_ini2(:,:,m) = nanmean(mPrcpNeighbors_ini(:,:,m:12:156),3); %3x3x12
end

    % Get elevation data for HAR
    addpath('/uufs/chpc.utah.edu/common/home/steenburgh-group6/Eric/HAR_Data')
    mElevations =  ncread('/uufs/chpc.utah.edu/common/home/steenburgh-group6/Eric/HAR_Data/har_d10km_static_hgt.nc','hgt'); %270 x 180
%   keyboard
%     mElevations = flipud(mElevations'); %why? orientation for a map but
%     not matching what was extraxted from precip
    mElevations_3x3 = mElevations(m_Idx-1:m_Idx+1, n_Idx-1:n_Idx+1);
else
    error('Error: You havent added functionality for that GCM yet. Way to go pinhead.')
end    
%% Fit regression line to temperature data
vPrcpLapseRates = zeros(12,1);
mPolyfit_prcp = zeros(12,2);
warning('off','all')

for t = 1:12
    % Extract prcp and elevation data for gridcell + neighbors
    vAltNeighbors_temp = mElevations_3x3(:); %9x1
    mPrcpNeighbors_temp = rot90(squeeze(mPrcpNeighbors_ini2(:,:,t))); %3x3
%     mPrcpNeighbors_temp = (squeeze(mPrcpNeighbors_ini2(:,:,t))); %3x3

    % Calculate lapse rate (mm P km^-1)
%       figure(t+20); scatter(vAltNeighbors_temp,mPrcpNeighbors_temp(:),'b'); 
        mPolyfit_prcp(t,:) = polyfit(vAltNeighbors_temp,mPrcpNeighbors_temp(:),1); %sort x to plot as line
        f = polyval(mPolyfit_prcp(t,:),vAltNeighbors_temp); 
   figure(t+40); plot(vAltNeighbors_temp,mPrcpNeighbors_temp(:),'o',vAltNeighbors_temp,f,'-'), legend('data','linear fit')
    kLapseRate = mPolyfit_prcp(t,1) * 1000;
    vPrcpLapseRates(t) = kLapseRate; 
end

warning('on','all')