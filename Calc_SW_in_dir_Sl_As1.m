function [ mSW_in_dir_Sl_As ] = Calc_SW_in_dir_Sl_As1( mSW_in_dir, mLat, mLong, mGlacAlt, mSlope, mAspect, sTime, mGlacMask, t )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nansum(mSW_in_dir(:)) > 0

%% Calculate date and time for solar radiation calculation

time.year = sTime.year(t);
time.month = sTime.month(t);
time.day = sTime.day(t);
time.hour = sTime.hour(t);
time.min = 0;
time.sec = 0;
time.UTC = 5;

% Preallocate for 'for loop'
mSunZenith = zeros(size(mGlacMask));
mSunAzimuth = zeros(size(mGlacMask));

% parfor m = 1:size(mGlacMask,1)
%     
%     vSunZenith_temp = zeros(1,size(mGlacMask,2));
%     vSunAzimuth_temp = zeros(1,size(mGlacMask,2));
%     
%     for n = 1:size(mGlacMask,2)
% 
%         if mGlacMask(m,n) == 1
% 
%             longitude = mLong(m,n);
%             latitude = mLat(m,n);
%             altitude = mGlacAlt(m,n);
%             
%             [zenith, azimuth] = sun_position_parfor(time,latitude,longitude,altitude);
%             
%             vSunZenith_temp(n) = zenith;
%             vSunAzimuth_temp(n) = azimuth;
% 
%         end
% 
%     end
%     
%     mSunZenith(m,:) = vSunZenith_temp;
%     mSunAzimuth(m,:) = vSunAzimuth_temp;
% 
% end
[mSunZenith, mSunAzimuth] = sun_position_parfor_matrix(time,mLat,mLong,mGlacAlt);

mCos_Theta = abs(cosd(mSlope) .* cosd(mSunZenith) + sind(mSlope) .* sind(mSunZenith) .* cosd(mSunAzimuth - mAspect));

else

mCos_Theta = 1;
    
end

mSW_in_dir_Sl_As = mCos_Theta .* mSW_in_dir;




















































end

