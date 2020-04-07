%   Extract climate data for a single 10km x 10km grid point for use with
%   surface and energy mass balance model
clear
close all

sDataDirectory = '/uufs/chpc.utah.edu/common/home/steenburgh-group6/Eric/HAR_Data';
addpath(genpath(sDataDirectory))

load UIB_triangle.mat
    polyin = polyshape(lon_full,lat_full); 
    
cd(sDataDirectory)
mLat = double(ncread('har_d10km_h_2d_t2_2001.nc', 'lat'));
mLong = double(ncread('har_d10km_h_2d_t2_2001.nc','lon'));

vLat = mLat(:);
vLong = mLong(:);

% [~, iLatLongIdx] = min((vLat - Lat).^2 + (vLong - Long).^2);
% [m_Idx, n_Idx] = ind2sub([270, 180], iLatLongIdx);

% Coordinates
in_uib = inpolygon(vLong,vLat,lon_full,lat_full);

plot(polyin); hold on
    plot(vLong(in_uib==1),vLat(in_uib==1),'r.')
    plot(vLong(in_uib==0),vLat(in_uib==0),'b.')
    
Idx = find(in_uib==1);

% Elev_HAR = nan(length(vLat),1);
% Precip_HAR = nan(length(vLat),131496); 

% Grid Elevation
ncid = netcdf.open('har_d10km_static_hgt.nc');
    mAWS_Alt = double(ncread('har_d10km_static_hgt.nc', 'hgt'));
    Data.AWS_Alt = mAWS_Alt(:); 
netcdf.close(ncid)

Elev_HAR = Data.AWS_Alt(Idx);

%% Precipitation (mm/HOUR!)
for j = 1:length(Idx)
    iLatLongIdx = Idx(j); %Long,Lat combos in UIB
    [m_Idx, n_Idx] = ind2sub([270, 180], iLatLongIdx); 

    Precip = [];
    for i = 0:14

        if i < 10
            sFileName = ['har_d10km_h_2d_prcp_200', num2str(i), '_rechunked.nc'];
        else
            sFileName = ['har_d10km_h_2d_prcp_20', num2str(i), '_rechunked.nc'];
        end
        
        ncid = netcdf.open(sFileName);
            if rem(i,4) == 0
                Leap = 1;
            else
                Leap = 0;
            end
            vPrecip_Data = double(squeeze(ncread(sFileName, 'prcp',[m_Idx,n_Idx,1], [1,1,8760+24*Leap], [1,1,1])));
            Precip = [Precip;vPrecip_Data];

        netcdf.close(ncid)

    end

    % Convert to m / hr
    Precip = Precip / 1000; %1/1/00-1/1/15 for each location of 270 x 180 (48600)
    Precip_HAR(j,:) = Precip'; %131496x1
    j
end

%% Create time vectors
t1 = datetime(2000,1,1,0,0,0);
t2 = datetime(2014,12,31,23,0,0);

tt = t1:hours(1):t2;  %131496

cd ../../../u6027899/sebm
varinfo=whos('Precip_HAR');
saveopt='';
if varinfo.bytes >= 2^31
  saveopt='-v7.3';
end

save('HAR_ElevPrecip_exact.mat',saveopt)
% load ../../steenburgh-group6/Ali2/HAR_ElevPrecip.mat
