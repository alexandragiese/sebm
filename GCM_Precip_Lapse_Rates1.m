function [ vPrcpLapseRates, mPolyfit_prcp ] = GCM_Precip_Lapse_Rates1( sGCM, sTempRes, kLat, kLong, stTime )
%% Extract precip lapse rates from HAR
pl=0;


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

    [~, iLatLongIdx] = min((vLat - kLat).^2 + (vLong - kLong).^2);
    [m_Idx, n_Idx] = ind2sub([270, 180], iLatLongIdx); %270, 180???
    
    m3P_3x3 = [];
    m3YN_3x3 = [];
vP = nan; %single-point testing
vYN = nan; %single-point testing
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
            
%problem line below!!!
if (180-m_Idx) >= 178
    m_Idx = 2; % If 178 doesn't work, try 177 in the line above, and 3 in this line
end

            vPrecip_Data_3x3 = double(squeeze(ncread(sFileName, 'prcp',[n_Idx-1,180-m_Idx,1], [3,3,8760+24*Leap], [1,1,1]))); %3x3 all hr in yr
   vPrecip_Data_1x1 = double(squeeze(ncread(sFileName, 'prcp',[n_Idx-1,180-m_Idx,1], [1,1,8760+24*Leap], [1,1,1]))); %1x1 all hr in yr (single-point testing)

            vPrecip_YN_3x3 = zeros(size(vPrecip_Data_3x3));
                vPrecip_YN_3x3(1,:,:) = spones(squeeze(vPrecip_Data_3x3(1,:,:)));
                vPrecip_YN_3x3(2,:,:) = spones(squeeze(vPrecip_Data_3x3(2,:,:)));
                vPrecip_YN_3x3(3,:,:) = spones(squeeze(vPrecip_Data_3x3(3,:,:)));
   vPrecip_YN_1x1 = spones(vPrecip_Data_1x1); %single-point testing

            m3P_3x3 = cat(3, m3P_3x3, vPrecip_Data_3x3); %gridded precip data around point of interest (3x3 in all hr for whole stTime)
            m3YN_3x3 = cat(3, m3YN_3x3, vPrecip_YN_3x3);

        vP  = [vP; vPrecip_Data_1x1];
        vYN = [vYN; vPrecip_YN_1x1];

        end
 vP = vP(2:end); %single-point testing
 vYN = vYN(2:end); %single-point testing
 
    % Sum Precipitation (monthly total & monthly count of +P hours) for each month
    m3MT_Precip    = nan(3,3,12); 
    m3MC_posPrecip = nan(3,3,12); 

% keyboard
% load take_P_local.mat
t=1;
    % fill 3D-matrices <-- check one by hand
    for y = 2001:2013
        for m = 1:12
            m3MT_Precip(:,:,t) = nansum(m3P_3x3(:,:,stTime.month==m & stTime.year==y),3); %mm each month over all time / # yr
            m3MC_posPrecip(:,:,t) = nansum(m3YN_3x3(:,:,stTime.month==m & stTime.year==y),3);
% vMT(t) = nansum(vP(stTime.month==t & stTime.year>=2001 & stTime.year<=2013))/13;
% vMC(t) = nansum(vYN(stTime.month==t & stTime.year>=2001 & stTime.year<=2013))/13;
            t = t+1;
        end
    end
% Plots based on single point (to show how prcp lapse rate calc works)
k=1;    
for y = 2001:2013
    P_Jan(k) = nansum(vP(stTime.month==1 & stTime.year==y) );
    P_Feb(k) = nansum(vP(stTime.month==2 & stTime.year==y) );
    P_Mar(k) = nansum(vP(stTime.month==3 & stTime.year==y) );
    P_Apr(k) = nansum(vP(stTime.month==4 & stTime.year==y) );   
    P_May(k) = nansum(vP(stTime.month==5 & stTime.year==y) );
    P_Jun(k) = nansum(vP(stTime.month==6 & stTime.year==y) );    
    P_Jul(k) = nansum(vP(stTime.month==7 & stTime.year==y) );
    P_Aug(k) = nansum(vP(stTime.month==8 & stTime.year==y) );    
    P_Sep(k) = nansum(vP(stTime.month==9 & stTime.year==y) );
    P_Oct(k) = nansum(vP(stTime.month==10 & stTime.year==y) );    
    P_Nov(k) = nansum(vP(stTime.month==11 & stTime.year==y) );
    P_Dec(k) = nansum(vP(stTime.month==12 & stTime.year==y) );  
    
    C_Jan(k) = nansum(vYN(stTime.month==1 & stTime.year==y) );
    C_Feb(k) = nansum(vYN(stTime.month==2 & stTime.year==y) );
    C_Mar(k) = nansum(vYN(stTime.month==3 & stTime.year==y) );
    C_Apr(k) = nansum(vYN(stTime.month==4 & stTime.year==y) );   
    C_May(k) = nansum(vYN(stTime.month==5 & stTime.year==y) );
    C_Jun(k) = nansum(vYN(stTime.month==6 & stTime.year==y) );    
    C_Jul(k) = nansum(vYN(stTime.month==7 & stTime.year==y) );
    C_Aug(k) = nansum(vYN(stTime.month==8 & stTime.year==y) );    
    C_Sep(k) = nansum(vYN(stTime.month==9 & stTime.year==y) );
    C_Oct(k) = nansum(vYN(stTime.month==10 & stTime.year==y) );    
    C_Nov(k) = nansum(vYN(stTime.month==11 & stTime.year==y) );
    C_Dec(k) = nansum(vYN(stTime.month==12 & stTime.year==y) );  
    k=k+1;
end   
if pl==1
close all
figure;
subplot(341); plot(P_Jan,'.'); title('Jan'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(342); plot(P_Feb,'.'); title('Feb'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(343); plot(P_Mar,'.'); title('Mar'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(344); plot(P_Apr,'.'); title('Apr'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(345); plot(P_May,'.'); title('May'); ylabel('Cumulative mm Precip'); 
                                             xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(346); plot(P_Jun,'.'); title('Jun'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(347); plot(P_Jul,'.'); title('Jul'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(348); plot(P_Aug,'.'); title('Aug'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(349); plot(P_Sep,'.'); title('Sep'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(3,4,10); plot(P_Oct,'.'); title('Oct'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(3,4,11); plot(P_Nov,'.'); title('Nov'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(3,4,12); plot(P_Dec,'.'); title('Dec'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})

figure;
subplot(341); plot(C_Jan,'.'); title('Jan'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(342); plot(C_Feb,'.'); title('Feb'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(343); plot(C_Mar,'.'); title('Mar'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(344); plot(C_Apr,'.'); title('Apr'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(345); plot(C_May,'.'); title('May'); ylabel('Hours of Precip (24*30=720)'); 
                                             xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(346); plot(C_Jun,'.'); title('Jun'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(347); plot(C_Jul,'.'); title('Jul'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(348); plot(C_Aug,'.'); title('Aug'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(349); plot(C_Sep,'.'); title('Sep'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(3,4,10); plot(C_Oct,'.'); title('Oct'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(3,4,11); plot(C_Nov,'.'); title('Nov'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
subplot(3,4,12); plot(C_Dec,'.'); title('Dec'); xticks([1 2 3 4 5 6 7 8 9 10 11 12 13]); xticklabels({'2001','','','','','06','','','','','11','',''})
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
    mPrcpNeighbors_ini2(:,:,m) = mean(mPrcpNeighbors_ini(:,:,m:12:156),3); %3x3x12
end

    % Get elevation data for HAR
    addpath('/uufs/chpc.utah.edu/common/home/steenburgh-group6/Eric/HAR_Data')
    mElevations =  ncread('/uufs/chpc.utah.edu/common/home/steenburgh-group6/Eric/HAR_Data/har_d10km_static_hgt.nc','hgt');
    mElevations = flipud(mElevations');
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
%     figure(t); scatter(vAltNeighbors_temp,mPrcpNeighbors_temp(:),'b'); 
        mPolyfit_prcp(t,:) = polyfit(vAltNeighbors_temp,mPrcpNeighbors_temp(:),1); %sort x to plot as line
        f = polyval(mPolyfit_prcp(t,:),vAltNeighbors_temp); 
%         plot(vAltNeighbors_temp,mPrcpNeighbors_temp(:),'o',vAltNeighbors_temp,f,'-'), legend('data','linear fit')
    kLapseRate = mPolyfit_prcp(t,1) * 1000;
    vPrcpLapseRates(t) = kLapseRate; 
end
warning('on','all')
