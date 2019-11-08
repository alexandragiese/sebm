function [ mSW_in_dir_Sl_As ] = Calc_SW_in_dir_Sl_As1( mSW_in_dir, mLat, mLong, mGlacAlt, mSlope, mAspect, sTime, mGlacMask, kTimeStep_days, sGCM, t )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



if nansum(mSW_in_dir(:)) > 0
    
    % Date and time for solar radiation calculation
    time.year = sTime.year(t);
    time.month = sTime.month(t);
    time.day = sTime.day(t);
    time.hour = sTime.hour(t);
    time.min = 0;
    time.sec = 0;
    
    % Set time zone for HAR vs FLOR
    if strcmp(sGCM,'HAR')
        time.UTC = 5; 
    elseif strcmp(sGCM,'FLOR')
        time.UTC = 0;
    end

    if kTimeStep_days < 1
    
        %% Hourly calculations

        [mSunZenith, mSunAzimuth] = sun_position_parfor_matrix(time,mLat,mLong,mGlacAlt);

        mCos_Theta = abs(cosd(mSlope) .* cosd(mSunZenith) + sind(mSlope) .* sind(mSunZenith) .* cosd(mSunAzimuth - mAspect));

    elseif kTimeStep_days >= 1

       %% Daily calculations
       
%        m3SW_in = zeros([size(mGlacMask),24]);
%        
%         for g = 1:24
%             
%             time.hour = sTime.hour(g);
%             
%             [mSunZenith, mSunAzimuth] = sun_position_parfor_matrix(time,mLat,mLong,mGlacAlt);
% 
%             mCos_Theta = abs(cosd(mSlope) .* cosd(mSunZenith) + sind(mSlope) .* sind(mSunZenith) .* cosd(mSunAzimuth - mAspect));
%             
%         end
        
    mCos_Theta = 1; % Not currently downscaling SW_in for daily runs
        
        
    end


else

mCos_Theta = 1;
    
end

mSW_in_dir_Sl_As = mCos_Theta .* mSW_in_dir;




















































end

