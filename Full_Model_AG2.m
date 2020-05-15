function [R1, R2, RAve] = Full_Model_AG2( GUI_Input )
threeDsave = GUI_Input.threeDsave;
close all
addpath('/uufs/chpc.utah.edu/common/home/u6027899/sebm');
% addpath('/uufs/chpc.utah.edu/common/home/u6027899/sebm/SkyViewFolder');

%% Load files
% mGlacNum    = load('ALOS_90m_Glacier_Numbers_Ind_Gang_Brahm1.mat','mGlacNum'); %***
mGlacNum    = load('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_30m_Numbers.mat','mGlacNum');
mGlacNum    = mGlacNum.mGlacNum;
% mDebrisMask = load('ALOS_90m_Debris_Mask_Ind_Gang_Brahm1_forAG.mat','mDebrisMask_90m'); % Added 11/15/19 -ESJ ***
% mDebrisMask = mDebrisMask.mDebrisMask_90m; % Added 11/15/19 -ESJ
mDebrisMask = load('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_30m_Debris_Mask.mat','mDebrisMask'); 
mDebrisMask = mDebrisMask.mDebrisMask; 

%% Specifics for larger region of interest

kColLn = size(mGlacNum,1); %ACTUALLY NUMBER OF ROWS, EJ: 15600; 
kRowLn = size(mGlacNum,2); %ACTYALLY NUMBER OF COL,  EJ: 45600; 

%% Glacier to extract
iGlacierNumber  = GUI_Input.glacier_number; 

% Border (in matrix indices) necessary for sky view factor calculation
kBorder         = 180;

% Find indices of glacier edges
[row, col]  = find(mGlacNum==iGlacierNumber); 
iM_i        = min(row) - kBorder;
iM_f        = max(row) + kBorder;
iN_i        = min(col) - kBorder;
iN_f        = max(col) + kBorder;

mGlacNum_Extract    = mGlacNum(iM_i:iM_f,iN_i:iN_f); %+1 removed
mDebrisMask_Extract = mDebrisMask(min(row):max(row),min(col):max(col)); 

% Extract DEMs and Masks for given glacier
mDEM    = ncread('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_30m_DEM.nc','mGlacAlt',[iM_i, iN_i],[iM_f-iM_i+1, iN_f-iN_i+1]); %+1 needed for ncread 'count'
mMask   = ncread('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_30m_Glacier_Mask.nc','mGlacMask',[iM_i, iN_i],[iM_f-iM_i+1, iN_f-iN_i+1]); %+1 needed for ncread 'count'

% Remove other glaciers from the glacier mask
mMask(mGlacNum_Extract ~= iGlacierNumber) = 0;

mGlacAlt    = double(mDEM);
mGlacMask   = mMask;
mDebMask    = mDebrisMask_Extract;

% Make lat/lon matrices
vLat    = (linspace(28.999861, 37.117361,kColLn))'; %28.999861 37.117361
mLat    = flipud(repmat(vLat,1,kRowLn));
mLat    = mLat(iM_i:iM_f,iN_i:iN_f); %+1 removed
vLong   = linspace(67.618750, 82.500694,kRowLn); %67.6187499 82.50069444
mLong   = repmat(vLong,kColLn,1);
mLong   = mLong(iM_i:iM_f,iN_i:iN_f); %+1 removed
% keyboard
kLat    = nanmean(mLat(mMask==1)); % Find average latitude of glacier
kLong   = nanmean(mLong(mMask==1)); % Find average longitude of glacier

clearvars mGlacNum mTileNumber row col iM_i iM_f iN_i iN_f ...
    vTileNumbers mGlacNum_Extract mDEM mMask t iTileNumber mM_Temp mN_Temp

%% Glacier slope and view factor

iLat_min    = min(mLat(:));
iLat_max    = max(mLat(:));
iLong_min   = min(mLong(:));
iLong_max   = max(mLong(:));

% Calculate geographic reference matrix for glacier mask
R_Geo = makerefmat('RasterSize',size(mGlacMask),'Latlim',[iLat_min, iLat_max],...
    'Lonlim',[iLong_min, iLong_max]);

% Slope of the glacier surface (degrees)
[mAspect, mSlope, ~, ~] = gradientm(mGlacAlt,R_Geo);
mAspect(isnan(mAspect)) = nanmean(mAspect(mGlacMask==1));
% Sky view factor (unitless fraction)
kRadius                 = GUI_Input.kRadius; %<-- what is this?
[mSVF]                  = SkyViewFactor1(mGlacMask,mGlacAlt,mLat,mLong,kRadius); % Gridcells (1 gridcell = 30m ALG)
mSVF(isnan(mSVF) == 1)  = 1;

%% Remove added borders (--> mGlacMask to final size)
mGlacMask(1:kBorder,:)  = []; mGlacMask(end-kBorder+1:end,:) = [];
mGlacMask(:,1:kBorder)  = []; mGlacMask(:,end-kBorder+1:end) = [];
mGlacAlt(1:kBorder,:)   = []; mGlacAlt(end-kBorder+1:end,:) = [];
mGlacAlt(:,1:kBorder)   = []; mGlacAlt(:,end-kBorder+1:end) = [];
mSlope(1:kBorder,:)     = []; mSlope(end-kBorder+1:end,:) = [];
mSlope(:,1:kBorder)     = []; mSlope(:,end-kBorder+1:end) = [];
mAspect(1:kBorder,:)    = []; mAspect(end-kBorder+1:end,:) = [];
mAspect(:,1:kBorder)    = []; mAspect(:,end-kBorder+1:end) = [];
mSVF(1:kBorder,:)       = []; mSVF(end-kBorder+1:end,:) = [];
mSVF(:,1:kBorder)       = []; mSVF(:,end-kBorder+1:end) = [];
mLat(1:kBorder,:)       = []; mLat(end-kBorder+1:end,:) = [];
mLat(:,1:kBorder)       = []; mLat(:,end-kBorder+1:end) = [];
mLong(1:kBorder,:)      = []; mLong(end-kBorder+1:end,:) = [];
mLong(:,1:kBorder)      = []; mLong(:,end-kBorder+1:end) = [];
  
%% Extract GCM Data    

sGCM = GUI_Input.sGCM;

if strcmp(sGCM,'HAR')
    sDataDirectory = '/uufs/chpc.utah.edu/common/home/steenburgh-group6/Eric/HAR_Data';
elseif strcmp(sGCM,'FLOR')
    sDataDirectory = '/uufs/chpc.utah.edu/common/home/steenburgh-group6/Eric/GFDL_FLOR_Data';
elseif strcmp(sGCM,'WRF')
    
end

sTempRes = GUI_Input.HAR_temp_res;
if strcmp(sGCM,'FLOR')
    iEnsNum         = GUI_Input.ensemble_number; % Only matters if sGCM = FLOR
    [ vLW_in, vP_d, vP_a, ~, vRH, vT_a, vSW_in, vU, stTime, Data ] = ExtractFLOR1(kLat, kLong, sDataDirectory, iEnsNum);
elseif strcmp(sGCM,'HAR')
    if strcmp(sTempRes,'Hourly') % I can never remember whether or not to capitalize things
        sTempRes    = 'hourly';
    elseif strcmp(sTempRes,'Daily')
        sTempRes    = 'daily';
    end
    [ vSW_in, vU, vT_a, vP_d, vP_a, vRH, ~, vLW_in, stTime, Data  ] = ExtractHAR1( kLat, kLong, sTempRes, sDataDirectory );
elseif strcmp(sGCM,'WRF')
    
end
kAWS_Alt = Data.AWS_Alt;

%% Start and End Dates

% Start Date
StartMonth  = str2double(GUI_Input.start_date(1:2));
StartDay    = str2double(GUI_Input.start_date(4:5));
StartYear   = str2double(GUI_Input.start_date(7:10));

kStartDay   = find(stTime.day==StartDay & stTime.month==StartMonth & stTime.year==StartYear,1,'first');

% End Date
EndMonth    = str2double(GUI_Input.end_date(1:2));
EndDay      = str2double(GUI_Input.end_date(4:5));
EndYear     = str2double(GUI_Input.end_date(7:10));

kEndDay     = find(stTime.day==EndDay & stTime.month==EndMonth & stTime.year==EndYear,1,'first');

% % Date to start recording MB (1975-2000)
% MBMonth1975 = str2double(GUI_Input.mb_date_1975(1:2));
% MBDay1975   = str2double(GUI_Input.mb_date_1975(4:5));
% MBYear1975  = str2double(GUI_Input.mb_date_1975(7:10));
% 
% % kMBDay1975 = find(stTime.day==MBDay1975 & stTime.month==MBMonth1975 & stTime.year==MBYear1975,1,'first');
% 
% % Date to start recording MB (2000-2016)
% MBMonth2000 = str2double(GUI_Input.mb_date_2000(1:2));
% MBDay2000   = str2double(GUI_Input.mb_date_2000(4:5));
% MBYear2000  = str2double(GUI_Input.mb_date_2000(7:10));
% 
% % kMBDay2000  = find(stTime.day==MBDay2000 & stTime.month==MBMonth2000 & stTime.year==MBYear2000,1,'first');

%% Time Step

% Data resolution (minutes)
kData_Resolution = Data.Resolution_mins;
% Time step (days)
kTimeStep_days  = kData_Resolution / (24 * 60);   
% Time step (seconds)
kTimeStep_sec   = kTimeStep_days * 24 * 60 * 60;
% Number of timesteps in one day
kN_per_Day      = 1 ./ kTimeStep_days;

%% Inputs

vTime           = 1:length(vT_a);     

% Density of the snow (kg m^-3)
kRho_snow       = GUI_Input.snow_density;

% Change in temperature, Forcing (C)
kDeltaT_a       = GUI_Input.temp_forcing;

% Roughness length
kZ_0m_snow  = GUI_Input.iZ_0m_snow;
kZ_0m_ice   = GUI_Input.iZ_0m_ice;
kZ_0T_snow  = GUI_Input.iZ_0T_snow;
kZ_0T_ice   = GUI_Input.iZ_0T_ice;
kZ_0q_snow  = GUI_Input.iZ_0q_snow;
kZ_0q_ice   = GUI_Input.iZ_0q_ice;

% Precipitation event threshold (m day^-1)
kPrecipThresh = GUI_Input.precip_threshold;

%% Lapse rate (C km^-1)
vLapseRates                      = GCM_Surface_Lapse_Rates1( sGCM, sTempRes, kLat, kLong );
[vPrcpLapseRates, mPolyfit_prcp] = GCM_Precip_Lapse_Rates1( sGCM, sTempRes, kLat, kLong, stTime );

[ vLapseRates_smooth ]  = SmoothLapseRates1( vLapseRates  );
vPrcpLapseRates_extended = [vPrcpLapseRates(:); vPrcpLapseRates(1)]; %Include Dec - Jan transition
vPrcpRates_smooth_extended = smooth(vPrcpLapseRates_extended); %to be smoothed
vPrcpLapseRates_smooth = GUI_Input.cf*vPrcpRates_smooth_extended(1:12); %
% figure; plot(vPrcpLapseRates,'.'); hold on; plot(vPrcpLapseRates_smooth)
%% Scale Wind Speed from 10m to 2m (for HAR)

if strcmp(sGCM,'HAR')
    kZ_0    = nanmean([kZ_0m_snow,kZ_0m_ice]);
    vU      = vU .* (log(2/kZ_0) ./ log(10/kZ_0));
end

%% Average weather station data inputs- n/a

% Assume no precip where data is missing
vP_d(isnan(vP_d)) = 0;
% Convert precip from mm hr^-1 to m hr^-1
% vP_d = vP_d / 1000;
% Only consider precip events greater than kPrecipThresh
% vP_d(vP_d < kPrecipThresh) = 0;

vT_a = vT_a + kDeltaT_a; %uniform T ch

%% Time Step

% Where to start and stop model (in elements in averaged vectors)
if strcmp(sGCM,'FLOR')
    [~, kStart] = min(abs(vTime-kStartDay*kN_per_Day-(kN_per_Day-1)));
    [~, kEnd]   = min(abs(vTime-kEndDay*kN_per_Day));
else
    kStart  = kStartDay;
    kEnd    = kEndDay;
end
                                        
%% Constants

% Latent heat of sublimation (J kg^-1)
kL_s        = 2.848e6;                                                          % From Rupper and Roe, 2008; Molg and Hardy, 2004
% Latent heat of vaporization for water (J kg^-1)
kL_v        = 2.514e6;                                                        	% From Rupper and Roe, 2008; Molg and Hardy, 2004
% Latent heat of fusion (J kg^-1)
kL_f        = 3.34e5;                                                           % From Rupper and Roe, 2008
% Albedo of fresh snow (unitless)
kAlpha_fs   = 0.85;                                                           	% 0.9 from Molg and Hardy, 2004; 0.85 from Physics of Glaciers
% Albedo of firn (unitless)
kAlpha_fi   = 0.40; %AG. was 0.55;                                              % 0.53 from Molg and Hardy, 2004; 0.55 from Physics of Glaciers
% Albedo of glacier ice (unitless)
kAlpha_i    = 0.30; %AG. was 0.35;                                              % 0.45 from Molg and Hardy, 2004; 0.35 from Physics of Glaciers
% e-folding constant for aging snow (days)
kE_t        = 21.9;                                                          	% From Molg and Hardy, 2004
% e-folding constant for snow depth (m)
kE_d        = 0.032;                                                            % From Molg and Hardy, 2004                                             
% Density of glacier ice (kg m^-3)
kRho_ice    = 900;

%% Make glacier matrices

mTotalMelt  = zeros(size(mGlacMask));

%% Initial snow depth across glacier
      
% Initial snow depth across the glacier is not constant (m)
mSnowDepth  = zeros(size(mGlacAlt)) + GUI_Input.snow_depth;

% Ice surface height (m)
mIceSurface = zeros(size(mGlacAlt));

%% Snow clock

mSnowClock  = zeros(size(mGlacMask)) - 1;
mTotalAccum = zeros(size(mGlacMask));
mTotalSub   = zeros(size(mGlacMask));
mTotalEvap  = zeros(size(mGlacMask));

%% Tracking Variables

% Precipitation fraction (snow/total)
RAve.PrecipFrac = NaN(length(vTime),1);

% 
% mNaNGlacMask = mGlacMask;
% mNaNGlacMask(mNaNGlacMask == 0) = nan;

%% Surface Temperature Constants

% Emissivity
kE_i        = 1; %AG. was 0.98;                                                               
% Stefan-Boltzmann constant (W m^-2 K^-4)
kSigma      = 5.6703737e-8;                                                   	% Google
% Air Density at standard sea level (kg m^-3)
kRho_a_sl   = 1.29;                                                          	% From Rupper and Roe, 2008
% Air pressure at sea level (hPa)
kP_a_sl     = 1013;                                                             % From Rupper and Roe, 2008
% Specific heat capacity of air at constant pressure (J kg^-1 K^-1)
kC_p        = 1010;                                                           	% From Rupper and Roe, 2008
% AWS measurement height (m)
kZ          = 2;                                                                % From Molg and Hardy, 2004
% Von Karman's constant (unitless)
kK_0        = 0.4;                                                          	% From Anderson and Others, 2010
% Gravity constant (m s^-2)
kG          = 9.8076;                                                        	% From Wikipedia
% Specific heat capacity of water (J kg^-1 K^-1)
kC_w        = 4181.3;                                                       	% Wikipedia, http://en.wikipedia.org/wiki/Heat_capacity
% Specific heat capacity of ice (J kg^-1 K^-1)
kC_i        = 2097;                                                          	% From Arnold et al., 2006
% Thermal conductivity of ice (W m^-1 K^-1)
kK          = 2.10;                                                         	% From Paterson, 1994; The Physics of Glaciers, 3rd ed., pg. 205   
% Precipitation phase transition threshold (C) [AG added this to lines
% where hard-coded as 0: approx. 427, 522, 566, 692, 696
kRainThreshold = 2;

%% Preallocate 2 layer transient model

mT_s_Act_s = zeros(size(mGlacMask));
mT_s_Act_2 = zeros(size(mGlacMask));
% mT_s_Act_3 = zeros(size(mGlacMask));

% Thickness of surface layers (m)
kLayerThick_s = 0.22;
kLayerThick_2 = 2.78;

% % Calculate temperature of the main body (DEEP) of the glacier (mean annual air temp) (C)
mT_s_Act_3 = (nanmean(vT_a) + (kAWS_Alt - mGlacAlt) .* nanmean(vLapseRates_smooth) / 1000) .* mGlacMask;
mT_s_Act_3(mT_s_Act_3 > 0) = 0;

% Thermal diffusivity of ice (m^2 s^-1)
kK_ice  = 1.16E-06;
% Thermal diffusivity of snow (m^2 s^-1)
kK_snow = 0.4E-06;

% Preallocate for 'for' loop
RAve.Q_L    = zeros(length(vTime),1);
RAve.Q_S    = zeros(length(vTime),1);
RAve.S_net  = zeros(length(vTime),1);
RAve.S_in   = zeros(length(vTime),1);
RAve.L_in   = zeros(length(vTime),1);
RAve.L_out  = zeros(length(vTime),1);
RAve.Q_P    = zeros(length(vTime),1);
RAve.Q_G    = zeros(length(vTime),1);
RAve.T_s_1  = zeros(length(vTime),1);
RAve.T_s_2  = zeros(length(vTime),1);
RAve.Q_m    = zeros(length(vTime),1);
RAve.Q_net  = zeros(length(vTime),1);
RAve.Sub    = zeros(length(vTime),1);
RAve.Q_T_s  = zeros(length(vTime),1); %AG
RAve.Q_T_2  = zeros(length(vTime),1); %AG
RAve.Lapse_Rate = zeros(length(vTime),1);
RAve.SnowDepth  = zeros(length(vTime),1);
RAve.SnowMelt   = zeros(length(vTime),1);
RAve.IceSurface = zeros(length(vTime),1);
RAve.IceMelt    = zeros(length(vTime),1);
RAve.Emissivity = zeros(length(vTime),1);
RAve.TotalMelt  = zeros(length(vTime),1);
RAve.DeltaM_net = zeros(length(vTime),1);
RAve.Alpha      = zeros(length(vTime),1);
RAve.T_a        = zeros(length(vTime),1);
RAve.Precip     = zeros(length(vTime),1);
RAve.P_a        = zeros(length(vTime),1);
RAve.Accumulation       = zeros(length(vTime),1);
RAve.Rain_Freeze_Energy = zeros(length(vTime),1);
RAve.Rain_Freeze_Amount = zeros(length(vTime),1);

if threeDsave == 1
% FOR SAVING FULLY DISTRIBUTED VARIABLES - PREALLOCATE SPACE
m3TotalMelt = zeros([size(mGlacMask),kEnd-kStart+1]); %Total Melt (mm/time)
m3TotalAccum= zeros([size(mGlacMask),kEnd-kStart+1]); %Total Acc (mm/time)
% m3Alpha     = zeros([size(mGlacMask),kEnd-kStart+1]); %Albedo (-)
% m3S_net     = zeros([size(mGlacMask),kEnd-kStart+1]); %New SW (W/m^2)
% m3L_out     = zeros([size(mGlacMask),kEnd-kStart+1]); %Emited LW (W/m^2)
% m3L_in      = zeros([size(mGlacMask),kEnd-kStart+1]); %Incoming LW (W/m^2)
% m3Q_S       = zeros([size(mGlacMask),kEnd-kStart+1]); %S (sensible?) (W/m^2)
% m3Q_L       = zeros([size(mGlacMask),kEnd-kStart+1]); %LE (latent?) (W/m^2)
% m3Q_P       = zeros([size(mGlacMask),kEnd-kStart+1]); %Energy from precip? (W/m^2)
% m3T_s_Act_s = zeros([size(mGlacMask),kEnd-kStart+1]); %Surface temperature (C)
% m3T_s_Act_2 = zeros([size(mGlacMask),kEnd-kStart+1]); %Temperature of 2nd surface layer (C)
% m3T_a       = zeros([size(mGlacMask),kEnd-kStart+1]); %Air temperature (C)
% m3Precipitation       = zeros([size(mGlacMask),kEnd-kStart+1]); %Total precip, both phases (mm w.e./time)
% m3Rain_Freeze_Amount  = zeros([size(mGlacMask),kEnd-kStart+1]); %Freezing rain (mm w.e./time?)
% m3Q_net     = zeros([size(mGlacMask),kEnd-kStart+1]); %Net surface energy (W/m^2)
% m3Q_m       = zeros([size(mGlacMask),kEnd-kStart+1]);     %Melt energy (W/m^2)
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t = kStart:kEnd 
%     if rem(t,11393)==0
%         disp(t)
%     end
    %% Extract timestep data

    mSW_in = vSW_in(t) .* mGlacMask;
    mSW_in = Calc_SW_in_dir_Sl_As1(mSW_in,mLat,mLong,mGlacAlt,mSlope,mAspect,stTime,mGlacMask,kTimeStep_days,sGCM,t);

    mU = vU(t) .* mGlacMask;

    %% Apply precipitation gradient (mm/m) up to ?    
    % no precipitation gradient used:
%     mPrecip = (vP_d(t) .* mGlacMask); 

% Monthly gradient (Precip increases with height and as a function of the amount of precip, unlike temp)
%   mPrecip = P @ HAR elev + elev dif * frac of average hourly P for month * P lapse rate for month (mm/m) .* mGlacMask;
    mPrecip = (vP_d(t) + (mGlacAlt - kAWS_Alt) * (vP_d(t)/polyval(mPolyfit_prcp(stTime.month(t),:),kAWS_Alt)) * vPrcpLapseRates_smooth(stTime.month(t)) / 1000) .* mGlacMask;  
%     save take_LR_local2.mat vP_d mGlacAlt kAWS_Alt vPrcpLapseRates_smooth stTime mPolyfit_prcp mGlacMask
    %% Calculate air pressure gradient with the hydrostatic equation (hPa/m)
    
    % Calculate day of the year 
    iDay       = datetime(stTime.year(t),stTime.month(t),stTime.day(t)); 
    iDay        = day(iDay,'dayofyear');
    iLapseRate  = vLapseRates_smooth(iDay);

    % Impose limits on lapse rates? (beta testing)
    if iLapseRate < 3
        iLapseRate = 3;
    elseif iLapseRate > 9.8
        iLapseRate = 9.8;
    end

    % Adjust air temperature for elevation (input lapse rate should be positive)
    mT_a = (vT_a(t) - (mGlacAlt - kAWS_Alt) * iLapseRate / 1000) .* ...
    mGlacMask;

    % Vapor pressure of the air
    mVP_a = (6.11 * exp((2.5e6 / 461) * (1 / 273 - 1 ./ (273 + mT_a)))) ...
    .* vRH(t) / 100 .* mGlacMask;

    % dP/dz = -P * g * (1 - e/p * (1-E))/(R_d*T_v)
    mPressGrad = -vP_a(t) .* 9.81 .* (1 - (mVP_a ./ ...
        vP_a(t)) * (1 - 0.622)) ./ (287.06 * (273.15 + mT_a)); %AG changed from 287

    % Apply air pressure gradient (input pressure gradient should be negative)
    mP_a = (vP_a(t) + mPressGrad .* (mGlacAlt - kAWS_Alt)) .* mGlacMask;         

    %% Precipitation
    
    % Update snow clock
    mSnowClock = (mSnowClock + kTimeStep_days) .* mGlacMask;
    
    % Add accumulation (m)
    mTotalAccum(mT_a <= kRainThreshold) = (mTotalAccum(mT_a <= kRainThreshold) + mPrecip(mT_a <= kRainThreshold));
    mTotalAccum = mTotalAccum .* mGlacMask;
    
    % Reset snow clock anywhere it snowed more than specified precitation threshold                      
    mSnowClock(mT_a <= kRainThreshold & mPrecip > kPrecipThresh) = 0;
    
    % If it snowed, add snow (only if T_a > kRainThreshold, i.e. snow, not rain) (m) - AG changed
    mSnowDepth(mT_a <= kRainThreshold) = (mSnowDepth(mT_a <= kRainThreshold) + mPrecip(mT_a <= kRainThreshold) ...
        * 1000 / kRho_snow) .* mGlacMask(mT_a <= kRainThreshold);
    %% Albedo
                     
    % Aging snow
    mAlpha_s    = kAlpha_fi + (kAlpha_fs - kAlpha_fi) .* exp(-mSnowClock / kE_t);                                              % From Molg and Hardy, 2004 (Eq. 3)
    % Snow depth albedo
    mAlpha      = (mAlpha_s + (kAlpha_i - mAlpha_s) .* exp(-mSnowDepth / ...
        kE_d)) .* mGlacMask;                                                % From Molg and Hardy, 2004 (Eq. 2)

    %% Stability correction for Q_L and Q_S
    
    % Distribute roughness lengths for snow vs ice
    mZ_0m = zeros(size(mGlacMask));
    mZ_0T = zeros(size(mGlacMask));
    mZ_0q = zeros(size(mGlacMask));
    
    mZ_0m(mSnowDepth > 0) = kZ_0m_snow;
    mZ_0T(mSnowDepth > 0) = kZ_0T_snow;
    mZ_0q(mSnowDepth > 0) = kZ_0q_snow;
    
    mZ_0m(mSnowDepth == 0) = kZ_0m_ice;
    mZ_0T(mSnowDepth == 0) = kZ_0T_ice;
    mZ_0q(mSnowDepth == 0) = kZ_0q_ice;

    % R_b (Richardson stability constant)
    mR_b = kG .* (mT_a - mT_s_Act_s) .* (kZ - mZ_0m) ./ ((273.15 + ...
        mT_a) .* mU.^2);                                                    % From Anderson and Others, 2010 (Eq. 17)
    % Stability correction (only applied for sub-daily time steps)
    if kTimeStep_days < 1
        mStabCorr = (1 - 5 * mR_b).^2;
        % Only use stability corrections if 0.2 > R_b > 0.01
        mStabCorr(mR_b < 0.01)  = 1;                                                 
        mStabCorr(mR_b > 0.2)   = 0;                                                  
        mStabCorr(isinf(mStabCorr) == 1) = 1;
    else
        mStabCorr = 1;                                                          % Remove stability correction for daily timestep
    end
    % k_H (Exchange coefficient for sensible heat)
    mK_H = kK_0^2 .* mStabCorr ./ (log(10 ./ mZ_0m) .* (log(kZ ./ mZ_0T)));           % From Anderson and Others, 2010 (Eq. 16); Braithwaite, 1995
    % k_E (exchange coefficient for latent heat)
    mK_E = kK_0^2 .* mStabCorr ./ (log(10 ./ mZ_0m) .* (log(kZ ./ mZ_0q)));           % From Anderson and Others, 2010 (Eq. 16); Braithwaite, 1995
    
    %% Sublimation vs evaporation
    
    % Latent heat of vaporization/fusion
    mL_vf = zeros(size(mGlacMask));
    mL_vf(mT_s_Act_s  >= 0) = kL_v; % Evaporation
    mL_vf(mT_s_Act_s < 0)   = kL_s; % Sublimation

    %% Saturation vapor pressure of the surface

    mSVP_s = 6.11 * exp((2.5e6 / 461) .* (1 / 273.15 - 1 ./ (273.15 + mT_s_Act_s)));

    %% Layer thicknesses, densities, and diffusivities
    
    % Calculate snow thicknesses in the top layer (m)
    mSnowThickness_s = mSnowDepth .* mGlacMask;
    mSnowThickness_s(mSnowThickness_s > kLayerThick_s) = kLayerThick_s;
    % Calculate snow thickness in the bottom layer (m)
    mSnowThickness_2 = (mSnowDepth - kLayerThick_s) .* mGlacMask;
    mSnowThickness_2(mSnowThickness_2 < 0) = 0;
    mSnowThickness_2(mSnowThickness_2 > kLayerThick_2) = kLayerThick_2;
    
    % Calculate amount of ice in each layer (m)
    mIceThickness_s = kLayerThick_s - mSnowThickness_s;
    mIceThickness_s(mIceThickness_s < 0) = 0;
    mIceThickness_2 = kLayerThick_2 - mSnowThickness_2;
    mIceThickness_2(mIceThickness_2 < 0) = 0;
    % Calculate surface densities
    mRho_s = (mSnowThickness_s .* kRho_snow + mIceThickness_s * kRho_ice) / (kLayerThick_s);
    mRho_2 = (mSnowThickness_2 .* kRho_snow + mIceThickness_2 * kRho_ice) / (kLayerThick_2);
    
    % Calculate net surface diffusivity (m^2 s^-1)
    mK_s = (mSnowThickness_s * kK_snow + mIceThickness_s * kK_ice) / (kLayerThick_s);
    mK_2 = (mSnowThickness_2 * kK_snow + mIceThickness_2 * kK_ice) / (kLayerThick_2);
    
    %% Calculate effective emissivity
    
    kE_eff_a = vLW_in(t) / (kSigma * (273.15 + vT_a(t)).^4);
        
    %% Energy balance

    % S(1-alpha) (W m^-2)
    mS_net = (mSW_in .* (1 - mAlpha)) .* mGlacMask;
    % L_out (W m^-2)
    mL_out = -(kSigma * kE_i * (273.15 + mT_s_Act_s).^4) .* mGlacMask;
    % L_in (W m^-2)
    mL_in = kSigma .* kE_eff_a .* ((273.15 + mT_a) .^ 4) .* mGlacMask;
    % Q_S
    mQ_S = (kRho_a_sl .* kC_p .* mK_H .* (mP_a ./ kP_a_sl) .* mU .* (mT_a - mT_s_Act_s)) ...
        .* mGlacMask;
    % Q_L
    mQ_L = 0.622 * kRho_a_sl .* mK_E .* mU .* mL_vf .* (mVP_a - ...
        mSVP_s) ./ mP_a .* mGlacMask;
    % Q_P (W m^-2)
    mQ_P = 1000 * kC_w * mPrecip .* (mT_a) / kTimeStep_sec .* mGlacMask;   % From Hock 2005
    % Precipitation heat flux only from rain (i.e. when T_a > kRainThreshold) - AG changed
    mQ_P(mT_a <= kRainThreshold) = 0;
    % Q_G (W m^-2)
    mQ_G = ( kK * mK_s .* (mT_s_Act_2 - mT_s_Act_s) ./ ...
        kLayerThick_s ) ./ kLayerThick_s .* kTimeStep_sec .* mGlacMask;     % From Paterson, 1994; The Physics of Glaciers, 3rd ed., pg. 206 (Eq. 5)
    % Remove NaNs and stuff
    mQ_L(isnan(mQ_L) == 1) = 0;

    % Total (AG) Energy (W m^-2)
    mQ_net = (mS_net + mL_out + mL_in + mQ_S + mQ_L + mQ_P + mQ_G) ...
        .* mGlacMask;

    %% Surface temperature: equations explained in equations 11 and 12 of Arnold et al., 2006
% The difference between the implicit method form (the one in the manuscript) and the 
% explicit form (the model that you're currently working with) is the incorporation of 
% equations 11 and 12 of Arnold et al., 2006.

    % Calculate change in temperature of top 2 layers for this timestep (K)   
    mDelta_T_s = (( mK_s .* (mT_s_Act_2 - mT_s_Act_s) ./ ...
        kLayerThick_s ) ./ kLayerThick_s + (mQ_net - mQ_G) ./ (kC_i .* ...
        mRho_s .* kLayerThick_s)) .* kTimeStep_sec .* mGlacMask;
    mDelta_T_2 = (  ( mK_2 .* (mT_s_Act_3 - mT_s_Act_2) ./ ...
        kLayerThick_2 ) - ( mK_s .* (mT_s_Act_2 - mT_s_Act_s) )...
        ./ kLayerThick_s  ) ./ kLayerThick_2 .* kTimeStep_sec .* mGlacMask;

    %%
    %Heat fluxes from changing surface temperature (AG) (W m^{-2})
    mQ_T_s = (( mK_s .* (mT_s_Act_2 - mT_s_Act_s) ./ ...
        kLayerThick_s )./kLayerThick_s + (mQ_net - mQ_G) ./ (kC_i .* ...
        mRho_s .* kLayerThick_s)) .* mGlacMask .*... 
        (kC_i .* mRho_s .* kLayerThick_s);
    mQ_T_2 = (( mK_2 .* (mT_s_Act_3 - mT_s_Act_2) ./ ...
        kLayerThick_2) - (mK_s .* (mT_s_Act_2 - mT_s_Act_s)) ...
        ./ kLayerThick_s  ) ./ kLayerThick_2 .* mGlacMask .*... 
        (kC_i .* mRho_s .* kLayerThick_2);
    
    % Theoretical surface temperature (K)
    mT_s_Theor_s = mT_s_Act_s + mDelta_T_s;
    mT_s_Theor_2 = mT_s_Act_2 + mDelta_T_2;
    mT_s_Theor_s(mT_s_Theor_s<-50) = -50;

    % Actual surface temperature (can't be above zero) (K)
    mT_s_Act_s = mT_s_Theor_s;
    mT_s_Act_2 = mT_s_Theor_2;
    mT_s_Act_s(mT_s_Act_s > 0) = 0;
    mT_s_Act_2(mT_s_Act_2 > 0) = 0;
        
    %% Freezing of rain  
    % Rain that falls when the cold content of the surface layer is zero 
    % is assumed to runoff the glacier.
    
    mRain_Freeze_Energy = zeros(size(mGlacMask));
    % Freezing of rain (W m^-2) %AG changed
    mRain_Freeze_Energy(mT_a > kRainThreshold & mPrecip > 0) = mPrecip(mT_a > kRainThreshold & mPrecip > 0) ...
        * 1000 * kL_f / kTimeStep_sec;
    % Cold content of the surface layer (W m^-2)
    mColdContent_s = kC_i * mRho_s .* (kLayerThick_s * mRho_s / 1000) .* ...
        (0 - mT_s_Act_s) / kTimeStep_sec;
    % Can't freeze more than cold content allows
    mRain_Freeze_Energy(mRain_Freeze_Energy > mColdContent_s) = ...
        mColdContent_s(mRain_Freeze_Energy > mColdContent_s);
    % Freezing water warms the surface layer
    mT_s_Act_s = mT_s_Act_s + ( mRain_Freeze_Energy ./ (kC_i .* mRho_s .* kLayerThick_s) .* kTimeStep_sec );
    % Add refreeze energy to Q_net (but not to Q_m)
    mQ_net = mQ_net + mRain_Freeze_Energy;
    % Calculate amount of frozen rain (m w.e.)
    mRain_Freeze_Amt = mRain_Freeze_Energy * kTimeStep_sec / (kL_f * 1000);
    % Add frozen rain as accumulation (m w.e.)
    mTotalAccum = (mTotalAccum + mRain_Freeze_Amt) .* mGlacMask;
    %% Melt energy (W m^-2)

    mQ_m = zeros(size(mGlacMask));
    mQ_m(mT_s_Theor_s > 0) = mT_s_Theor_s(mT_s_Theor_s > 0) * kC_i .* ...
        mRho_s(mT_s_Theor_s > 0) * kLayerThick_s / kTimeStep_sec;

    %% Melting
    % Snow melt (m)
    mSnowMelt = mQ_m * kTimeStep_sec / (kL_f * kRho_snow);

    % Adjust snow depth for melt (m)
    mSnowDepth = (mSnowDepth - mSnowMelt) .* mGlacMask;            
    
    %% If no snow, melt ice
    
    % Find negative snow depths
    mXS_SnowMelt = zeros(size(mGlacMask));
    mXS_SnowMelt(mSnowDepth < 0) = mSnowDepth(mSnowDepth < 0);
    % Set negative snow depths to zero
    mSnowDepth(mSnowDepth < 0) = 0;
    % Ice melt (m)
    mIceMelt = -mXS_SnowMelt * kRho_snow / kRho_ice;
    % Cumulative ice melt (m)
    mIceSurface = mIceSurface - mIceMelt;
    
    %% Sublimation and evaporation
    
    % Sublimation (m)
    mSub = zeros(size(mGlacMask));
    mSub(mQ_L < 0 & mT_s_Act_s < 0) = -mQ_L(mQ_L < 0 & mT_s_Act_s < 0) .* kTimeStep_sec ./ (kL_s * kRho_snow) .* mGlacMask(mQ_L < 0 & mT_s_Act_s < 0);  
    
    % Evaporation (m)
    mEvap = zeros(size(mGlacMask));
    mEvap(mQ_L < 0 & mT_s_Act_s >= 0) = -mQ_L(mQ_L < 0 & mT_s_Act_s >= 0) .* kTimeStep_sec ./ (kL_v * kRho_snow) .* mGlacMask(mQ_L < 0 & mT_s_Act_s >= 0);  
    
    % Adjust snow depth for sublimation (m)
    mSnowDepth = (mSnowDepth - mSub) .* mGlacMask; 
    
    %% If no snow, melt ice
    
    % Find negative snow depths
    mXS_SnowSub = zeros(size(mGlacMask));
    mXS_SnowSub(mSnowDepth < 0) = mSnowDepth(mSnowDepth < 0);
    % Set negative snow depths to zero
    mSnowDepth(mSnowDepth < 0) = 0;
    % Ice melt (m)
    mIceSub = -mXS_SnowSub * kRho_snow / kRho_ice;
    % Cumulative ice melt (m)
    mIceSurface = mIceSurface - mIceSub;
    
    %% Melt/refreeze in the second layer
    
%     mRefreeze_Energy = zeros(size(mGlacMask));
%     % Melt/refreeze energy (W m^-2)
%     mRefreeze_Energy(mNewMelt > 0) = mPrecip(mNewMelt > 0) * 1000 * kL_f / (kTS_mins * 60);
%     % Cold content of the surface layer (W m^-2)
%     mColdContent_2 = kC_i * mRho_s .* (kLayerThick_2 * mRho_2 / 1000) .* ...
%         (0 - mT_s_Act_2) / (kTS_mins * 60);
%     % Energy difference (between cold content and potential refreeze) (W m^-2)
%     mEn_Diff = zeros(size(mGlacMask));
%     mEn_Diff(mRefreeze_Energy > mColdContent_2) = ...
%         mRefreeze_Energy(mRefreeze_Energy > mColdContent_2) - ...
%         mColdContent_2(mRefreeze_Energy > mColdContent_2);
%     % Can't freeze more than cold content allows
%     mRefreeze_Energy(mRefreeze_Energy > mColdContent_2) = ...
%         mColdContent_2(mRefreeze_Energy > mColdContent_2);
%     % Freezing water warms the surface layer
%     mT_s_Act_s = mT_s_Act_s + ( mRefreeze_Energy ./ (kC_i .* mRho_s .* kLayerThick_s) .* (kTS_mins * 60) );
%     % Excess water when cold content is zero runs off as melt (m w.e.)
%     mRunoff = mEn_Diff * (kTS_mins * 60) / (kL_f * 1000);

    
    %% Net changes
    
    % Total average glacier melt (m w.eq.)
    mTotalMelt = mTotalMelt - (((mSnowMelt * kRho_snow / 1000) + ...
        (mIceMelt * kRho_ice / 1000)));
    % Total sublimation (m w.e.)
    mTotalSub = mTotalSub - (mSub .* kRho_snow / 1000);
    % Total evaporation (m w.e.)
    mTotalEvap = mTotalEvap - (mEvap .* kRho_snow / 1000);
    % Total mass balance (m w.eq.)
    mDeltaM_net = (mTotalMelt + mTotalAccum + mTotalSub);                   % Evaporation not included because melt is assumed to runoff currently    
    
    %% Average variables across glacier surface
    RAve.Q_L(t) = nanmean(mQ_L(mGlacMask==1));
    RAve.Q_S(t) = nanmean(mQ_S(mGlacMask==1));
    RAve.S_net(t) = nanmean(mS_net(mGlacMask==1));
    RAve.S_in(t) = vSW_in(t);
    RAve.L_in(t) = nanmean(mL_in(mGlacMask==1));
    RAve.L_out(t) = nanmean(mL_out(mGlacMask==1));
    RAve.Q_P(t) = nanmean(mQ_P(mGlacMask==1));
    RAve.Q_G(t) = nanmean(mQ_G(mGlacMask==1));
    RAve.Q_m(t) = nanmean(mQ_m(mGlacMask==1));
    RAve.Q_net(t) = nanmean(mQ_net(mGlacMask==1));
    RAve.Q_T_s(t) = nanmean(mQ_T_s(mGlacMask==1)); %AG
    RAve.Q_T_2(t) = nanmean(mQ_T_2(mGlacMask==1)); %AG
    RAve.T_s_1(t) = nanmean(mT_s_Act_s(mGlacMask==1));
    RAve.T_s_2(t) = nanmean(mT_s_Act_2(mGlacMask==1));
    RAve.Lapse_Rate(t) = iLapseRate;
    RAve.PLapse_Rate(t) = vPrcpLapseRates_smooth(stTime.month(t));
    RAve.SnowDepth(t) = nanmean(mSnowDepth(mGlacMask==1));
    RAve.SnowMelt(t) = nanmean(mSnowMelt(mGlacMask==1));
    RAve.Sub(t) = nanmean(mSub(mGlacMask==1));
    RAve.IceSurface(t) = nanmean(mIceSurface(mGlacMask==1));
    RAve.IceMelt(t) = nanmean(mIceMelt(mGlacMask==1));
    RAve.TotalMelt(t) = nanmean(mTotalMelt(mGlacMask==1));
    RAve.Alpha(t) = nanmean(mAlpha(mGlacMask==1));
    RAve.T_a(t) = nanmean(mT_a(mGlacMask==1));
    RAve.Precip(t) = nanmean(mPrecip(mGlacMask==1));
    RAve.P_a(t) = nanmean(mP_a(mGlacMask==1));
    RAve.Accumulation(t) = nansum(nansum(mPrecip(mT_a <= kRainThreshold & mGlacMask == 1))) / nansum(mGlacMask(:)); %AG changed
    RAve.Rain_Freeze_Energy(t) = nanmean(mRain_Freeze_Energy(mGlacMask==1));
    RAve.Rain_Freeze_Amount(t) = nanmean(mRain_Freeze_Amt(mGlacMask==1));
    
    % Precipitation fraction that falls as snow
    if nansum(nansum(mPrecip)) ~= 0
        RAve.PrecipFrac(t) = nansum(nansum(mGlacMask(mT_a <= kRainThreshold & mPrecip > 0))) / nansum(nansum(mGlacMask(mPrecip > kRainThreshold))); %AG changed
    elseif nansum(nansum(mPrecip)) == 0
        RAve.PrecipFrac(t) = nan;
    end

if threeDsave == 1
%FOR SAVING FULLY DISTRIBUTED VARIABLES
    m3TotalMelt(:,:,t) = mTotalMelt;
    m3TotalAccum(:,:,t) = mTotalAccum;
%     m3Alpha(:,:,t) = mAlpha;
%     m3S_net(:,:,t) = mS_net;
%     m3L_out(:,:,t) = mL_out;
%     m3L_in(:,:,t) = mL_in;
%     m3Q_S(:,:,t) = mQ_S;
%     m3Q_L(:,:,t) = mQ_L;
%     m3Q_P(:,:,t) = mQ_P;
%     m3T_s_Act_s(:,:,t) = mT_s_Act_s;
%     m3T_s_Act_2(:,:,t) = mT_s_Act_2;
%     m3T_a(:,:,t) = mT_a;
%     m3Precipitation(:,:,t) = mPrecip;
%     m3Rain_Freeze_Amount(:,:,t) = mRain_Freeze_Amt;
%     m3Q_net(:,:,t) = mQ_net;
%     m3Q_m(:,:,t) = mQ_m;
end


%     if (t == kMBDay1975) & (strcmp(sGCM,'FLOR')) %#ok
%        mTotalMelt = zeros(size(mTotalMelt));
%        mTotalAccum = zeros(size(mTotalAccum));
%     elseif (t == kMBDay2000) & (strcmp(sGCM,'FLOR')) %#ok
%        R1.TotalAveMelt = nanmean(mTotalMelt(mGlacMask == 1));
%        R1.TotalAveAccum = nanmean(mTotalAccum(mGlacMask == 1));
%        R1.TotalAveSub = nanmean(mTotalSub(mGlacMask == 1));
%        R1.TotalAveEvap = nanmean(mTotalEvap(mGlacMask == 1));
%        mTotalMelt_1975 = mTotalMelt;
%        mTotalAccum_1975 = mTotalAccum;
%        mTotalMelt = zeros(size(mTotalMelt));
%        mTotalAccum = zeros(size(mTotalAccum));
%        mTotalSub = zeros(size(mTotalSub));
%        mTotalEvap = zeros(size(mTotalEvap));
%     end

end    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Total Mass changes

if ~exist('R1','var')
    R1 = nan;
end

% Total average melt (m w.eq.)
R2.TotalAveMelt = nanmean(mTotalMelt(mGlacMask == 1));
% Total accumulation (m w.eq.)
R2.TotalAveAccum = nanmean(mTotalAccum(mGlacMask == 1));
% Total average sublimation (m w.eq.)
R2.TotalAveSub = nanmean(mTotalSub(mGlacMask == 1));
% Total average evaporation (m w.eq.)
R2.TotalAveEvap = nanmean(mTotalEvap(mGlacMask == 1));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Display total average melt (m)
% disp(['Total Average Melt = ',num2str(R.TotalAveMelt),' meters'])
% % Display total accumulation (m)
% disp(['Total average accumulation ',num2str(R.TotalAveAccum),' meters'])

clearvars -except R1 R2 RAve R_Geo GUI_Input iGlacierNumber sMicroPhysics ...
    mTotalMelt mTotalAccum m3Alpha m3S_net m3L_out m3L_in m3Q_S m3Q_L m3Q_P ...
    m3T_s_Act_s m3T_s_Act_2 m3T_a m3Precip m3Q_net m3Q_m m3T_a m3Precipitation ...
    m3Rain_Freeze_Amount m3TotalMelt m3TotalAccum mGlacMask mGlacAlt stTime ...
    mTotalMelt_1975 mTotalAccum_1975 mTotalSub mTotalEvap Data mDebMask kAWS_Alt

% These lines convert the 3D variables from double precision to single                                                                                   785
% precision to save space and time.  
% m3TotalMelt = single(m3TotalMelt);
% m3TotalAccum = single(m3TotalAccum);
% m3Alpha = single(m3Alpha);
% m3S_net = single(m3S_net);
% m3L_out = single(m3L_out);
% m3L_in = single(m3L_in);
% m3Q_S = single(m3Q_S);
% m3Q_L = single(m3Q_L);
% m3Q_P = single(m3Q_P);
% m3T_s_Act_s = single(m3T_s_Act_s);
% m3T_s_Act_2 = single(m3T_s_Act_2);
% m3T_a = single(m3T_a);
% m3Precipitation = single(m3Precipitation);
% m3Rain_Freeze_Amount = single(m3Rain_Freeze_Amount);

%         if strcmp(GUI_Input.sGCM,'HAR')
%             
%             save(['/uufs/chpc.utah.edu/common/home/u0929154/Desktop/Glacier_Model/Files_For_Ali/Results/Josh_MB_',GUI_Input.sGCM,'_',GUI_Input.HAR_temp_res,'_Glacier_Number_',num2str(iGlacierNumber),'.mat'])
%                 
%         elseif strcmp(GUI_Input.sGCM,'FLOR')
%             
%             save(['/uufs/chpc.utah.edu/common/home/u0929154/Desktop/Glacier_Model/Files_For_Ali/Results/Josh_MB_',GUI_Input.sGCM,'_Glacier_Number_',num2str(iGlacierNumber),'.mat'])
%                 
%         end

% else
%     
% R = 1;
% RAve = 1;
% 
%   clearvars -except R RAve R_Geo GUI_Input iTile_Number sMicroPhysics ...
%     mTotalMelt mTotalAccum

% save([GUI_Input.output_filename, sMicroPhysics,'_', num2str(iTile_Number),'.mat'])  

% end

% varinfo=whos('m3TotalMelt');
% saveopt='';
% if varinfo.bytes >= 2^31
%   saveopt='-v7.3';
% end

% IF BIG: save([GUI_Input.output_filename, num2str(iGlacierNumber),'_AG.mat'],saveopt)
% save([GUI_Input.output_filename num2str(iGlacierNumber),'_fluxtest.mat'])
% ej_var=nanmean(mGlacAlt(mGlacMask==1))-kAWS_Alt

end
