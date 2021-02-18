function [R1, R2, RAve ] = Full_Model_AG_edits( GUI_Input,y,g )

% close all
addpath('/uufs/chpc.utah.edu/common/home/u6027899/sebm');

%% Load files
mGlacNum    = load(GUI_Input.GlacNum_filename,'mGlacNum');
mGlacNum    = mGlacNum.mGlacNum;
% mDebrisMask = load('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_90m_Debris_Mask.mat','mDebrisMask'); 
% mDebrisMask = mDebrisMask.mDebrisMask; 

%% MOVED B/C OF PARFOR LOOP
glID  = GUI_Input.geo(:,1);
geoMB = GUI_Input.geo(:,11); % m w.e. / a
geoMB_sigma = GUI_Input.geo(:,12);
basin = GUI_Input.geo(:,49);

glacier_number = y(g); 
r = find (glID==y(g)); %index in geodetic MB spreadsheet

p1 = GUI_Input.PLR(basin(r),1);
p2 = GUI_Input.PLR(basin(r),2);

%% Specifics for larger region of interest

kColLn = size(mGlacNum,1); % NUMBER OF ROWS (misnamed)
kRowLn = size(mGlacNum,2); % NUMBER OF COL  (misnamed)

%% Glacier to extract
iGlacierNumber  = glacier_number; 

% Border (in matrix indices) necessary for sky view factor calculation
kBorder         = 180;

% Find indices of glacier edges
[row, col]  = find(mGlacNum==iGlacierNumber); 
iM_i        = min(row) - kBorder;
iM_f        = max(row) + kBorder;
iN_i        = min(col) - kBorder;
iN_f        = max(col) + kBorder;
if iM_i < 1 
    R2.TotalAveMelt = nan;
    R2.TotalAveAccum = nan;
    R2.TotalAveSub = nan;
    RAve = nan;
    fprintf('mGlacNum issue at: %d\n',g)
    return
end
if iN_i < 1
    R2.TotalAveMelt = nan;
    R2.TotalAveAccum = nan;
    R2.TotalAveSub = nan;
    RAve = nan;
    fprintf('mGlacNum issue at: %d\n',g)
    return
end
if iN_f >  17858   
    R2.TotalAveMelt = nan;
    R2.TotalAveAccum = nan;
    R2.TotalAveSub = nan;
    RAve = nan;
    fprintf('mGlacNum issue at: %d\n',g)
    return
end
if isempty(iM_i) || isempty(iN_i)
    R2.TotalAveMelt = nan;
    R2.TotalAveAccum = nan;
    R2.TotalAveSub = nan;
    RAve = nan;
    fprintf('mGlacNum issue at: %d\n',g)
    return
end

mGlacNum_Extract    = mGlacNum(iM_i:iM_f,iN_i:iN_f); %+1 removed
% mDebrisMask_Extract = mDebrisMask(min(row):max(row),min(col):max(col)); 

% Extract DEMs and Masks for given glacier
mDEM    = ncread('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_90m_DEM.nc','mGlacAlt',[iM_i, iN_i],[iM_f-iM_i+1, iN_f-iN_i+1]); %+1 needed for ncread 'count'
mMask   = ncread('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_90m_Glacier_Mask.nc','mGlacMask',[iM_i, iN_i],[iM_f-iM_i+1, iN_f-iN_i+1]); %+1 needed for ncread 'count'

% Remove other glaciers from the glacier mask
mMask(mGlacNum_Extract ~= iGlacierNumber) = 0;

mGlacAlt    = double(mDEM);
mGlacMask   = mMask;
% mDebMask    = mDebrisMask_Extract;

% Make lat/lon matrices
vLat    = (linspace(28.999861, 37.117361,kColLn))'; 
mLat    = flipud(repmat(vLat,1,kRowLn));
mLat    = mLat(iM_i:iM_f,iN_i:iN_f); %+1 removed
vLong   = linspace(67.618750, 82.500694,kRowLn); %67.6187499 82.50069444
mLong   = repmat(vLong,kColLn,1);
mLong   = mLong(iM_i:iM_f,iN_i:iN_f); %+1 removed
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
[mSVF]                  = SkyViewFactor1(mGlacMask,mGlacAlt,mLat,mLong,kRadius); % Gridcells 
mSVF(isnan(mSVF))       = 1;

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
    [ vSW_in, vU, vT_a, vP_d, vP_a, vRH, ~, vLW_in, stTime, Data  ] = ExtractHAR1( kLat, kLong, sTempRes, sDataDirectory ); %pulls in all HAR output
elseif strcmp(sGCM,'WRF')
    
end

kAWS_Alt = Data.AWS_Alt; 

vP_d_th = nan(size(vP_d)); 

for r = 1:length(vP_d)
    precip_series = reshape(vP_d,24,[]);
    daily_precip  = sum(precip_series); %m/d
    A = find(daily_precip <= GUI_Input.precip_threshold); % Precipitation event threshold (m day^-1)
    precip_series(:,A) = 0;
    vP_d_th = precip_series(:);    
end

vP_d = vP_d_th;

% FOR PRECIP ADJUSTMENT
hr_count = nan(13,1); 
for i = 1:13 %2000-2013
  hr_count(i) = sum(vP_d(stTime.year==2000+i)~=0);
end
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
% kDeltaT_a       = GUI_Input.temp_forcing;

% Roughness length
kZ_0m_snow  = GUI_Input.iZ_0m_snow;
kZ_0m_ice   = GUI_Input.iZ_0m_ice;
kZ_0T_snow  = GUI_Input.iZ_0T_snow;
kZ_0T_ice   = GUI_Input.iZ_0T_ice;
kZ_0q_snow  = GUI_Input.iZ_0q_snow;
kZ_0q_ice   = GUI_Input.iZ_0q_ice;

%% Lapse rate (C km^-1)
vLapseRates             = GCM_Surface_Lapse_Rates1( sGCM, sTempRes, kLat, kLong ); 

[ vLapseRates_smooth ]  = SmoothLapseRates1( vLapseRates  ); 

%% Scale Wind Speed from 10m to 2m (for HAR)
if strcmp(sGCM,'HAR')
    kZ_0    = nanmean([kZ_0m_snow,kZ_0m_ice]);
    vU      = vU .* (log(2/kZ_0) ./ log(10/kZ_0));
end

%% Average weather station data inputs- n/a
% Assume no precip where data is missing
vP_d(isnan(vP_d)) = 0;

% vT_a = vT_a + kDeltaT_a; %uniform T ch

%% Time Step

% Where to start and stop model (in elements in averaged vectors)
if strcmp(sGCM,'FLOR')
    [~, kStart] = min(abs(vTime-kStartDay*kN_per_Day-(kN_per_Day-1)));
    [~, kEnd]   = min(abs(vTime-kEndDay*kN_per_Day));
else %HAR!
    kStart  = kStartDay;  
    kEnd    = kEndDay;  
        if strcmp(sTempRes,'hourly')
                kEnd    = kEndDay+23;   %note midnight assigned to previous day in model
        end
end

%% Constants

% Latent heat of sublimation (J kg^-1)
kL_s        = 2.848e6;                                                          % From Rupper and Roe, 2008; Molg and Hardy, 2004
% Latent heat of vaporization for water (J kg^-1)
kL_v        = 2.514e6;                                                        	% From Rupper and Roe, 2008; Molg and Hardy, 2004
% Latent heat of fusion (J kg^-1)
kL_f        = 3.34e5;                                                           % From Rupper and Roe, 2008
% Albedo of fresh snow (unitless)
kAlpha_fs   = 0.75; % was 0.85                                                  % 0.9 from Molg and Hardy, 2004; 0.85 from Physics of Glaciers
% Albedo of firn (unitless)
kAlpha_fi   = 0.53; %AG. was 0.55, 0.40                                         % 0.53 from Molg and Hardy, 2004; 0.55 from Physics of Glaciers; 0.40 in Johnson & Rupper
% Albedo of glacier ice (unitless)
kAlpha_i    = 0.34; %AG. was 0.35, 0.30                                         % 0.45 from Molg and Hardy, 2004; 0.35 from Physics of Glaciers; 0.30 in Johnson & Rupper
% e-folding constant for aging snow (days)
kE_t        = 21.9;                                                          	% From Molg and Hardy, 2004
% e-folding constant for snow depth (m)
kE_d        = 0.032;                                                            % From Molg and Hardy, 2004                                             
% Density of glacier ice (kg m^-3)
kRho_ice    = 900;

%% Make glacier matrices

mTotalMelt  = zeros(size(mGlacMask));
mTotalSnowMelt  = zeros(size(mGlacMask));
mTotalIceMelt  = zeros(size(mGlacMask));

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

%% Surface Temperature Constants

% Emissivity
kE_i        = 1; %AG. was 0.98;                                                 % 0.98 in Johnson & Rupper                                                     
% Stefan-Boltzmann constant (W m^-2 K^-4)
kSigma      = 5.6703737e-8;                                                   	% Google
% Air Density at standard sea level (kg m^-3)
kRho_a_sl   = 1.29;                                                          	% From Rupper and Roe, 2008
% Air pressure at sea level (hPa)
kP_a_sl     = 1013;                                                             % From Rupper and Roe, 2008
% Specific heat capacity of air at constant pressure (J kg^-1 K^-1)
kC_a        = 1010;                                                           	% From Rupper and Roe, 2008
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
% Precipitation phase transition threshold (C) [AG added this to lines where hard-coded as 0: approx. lines 427, 522, 566, 692, 696]
kRainThreshold = 2;

%% Preallocate 2 layer transient model
mT_s = zeros(size(mGlacMask));
mT_2 = zeros(size(mGlacMask));

% Thickness of surface layers (m)
kLayerThick_s = 0.22;
% kLayerThick_2 = 2.78;
kLayerThick_2 = 9.78; %Also tried to have T_b at depth 10m (mean Tair, consistent w/ Paterson; -1.2 C, see below)

% % Calculate temperature of the main body (DEEP) of the glacier (mean annual air temp) (C)
mT_b = (nanmean(vT_a) + (kAWS_Alt - mGlacAlt) .* nanmean(vLapseRates_smooth) / 1000) .* mGlacMask;
mT_b(mT_b > 0) = 0; 
% mT_b = -1.2 .* mGlacMask; % Molg & Hardy have at -1.2C at 10m depth

% Thermal diffusivity of ice (m^2 s^-1)
kK_ice  = 1.16E-06;
% Thermal diffusivity of snow (m^2 s^-1)
kK_snow = 0.4E-06;

% Preallocate for loop
RAve.Q_L    = nan(length(vTime),1);
RAve.Q_S    = nan(length(vTime),1);
RAve.S_net  = nan(length(vTime),1);
RAve.S_in   = nan(length(vTime),1);
RAve.S_in_ds= nan(length(vTime),1);
RAve.L_in   = nan(length(vTime),1);
RAve.L_out  = nan(length(vTime),1);
RAve.Q_P    = nan(length(vTime),1);
RAve.Q_G    = nan(length(vTime),1);
RAve.Q_m    = nan(length(vTime),1);
RAve.Q_net  = nan(length(vTime),1);
RAve.T_a    = nan(length(vTime),1);
RAve.T_s    = nan(length(vTime),1);
RAve.T_2    = nan(length(vTime),1);
RAve.T_b    = nan(length(vTime),1);
RAve.Rain_Freeze_Energy = zeros(length(vTime),1);
RAve.Rain_Freeze_Amount = zeros(length(vTime),1);
RAve.Alpha      = nan(length(vTime),1);
RAve.Lapse_Rate = nan(length(vTime),1);
RAve.TotalMelt  = nan(length(vTime),1);
RAve.TotalSub   = nan(length(vTime),1);
RAve.TotalEvap  = nan(length(vTime),1);
RAve.TotalAccum = nan(length(vTime),1);

RAve.SnowDepth  = nan(length(vTime),1);
RAve.SnowMelt   = nan(length(vTime),1);
RAve.IceSurface = nan(length(vTime),1);
RAve.IceMelt    = nan(length(vTime),1);
RAve.Pr_HAR   = nan(length(vTime),1);  
RAve.Pr_glac  = nan(length(vTime),1); 
RAve.Pr_hravg = nan(length(vTime),1);  

RAve.DeltaTs = nan(length(vTime),1);
RAve.DeltaT2 = nan(length(vTime),1);

RAve.RhoSnow  = kRho_snow;
RAve.RhoIce   = kRho_ice;
RAve.GlacMask = size(mGlacMask);
RAve.res      = GUI_Input.GlacNum_filename;

TH_i.Accum = nan(size(mGlacMask));

    RAve.U = nan(length(vTime),1);
    RAve.Prcp = nan(length(vTime),1);
    RAve.P_a = nan(length(vTime),1);
    RAve.VP_a = nan(length(vTime),1);
    
% FOR SAVING FULLY DISTRIBUTED VARIABLES -----------------------------------------------------------------
if GUI_Input.threeDsave == 1
m3TotalMelt = nan([size(mGlacMask),kEnd-kStart+1]); %Total Melt (mm/time) *
m3TotalAccum= nan([size(mGlacMask),kEnd-kStart+1]); %Total Acc (mm/time)  *
m3Alpha     = nan([size(mGlacMask),kEnd-kStart+1]); %Albedo (-)
m3S_net     = nan([size(mGlacMask),kEnd-kStart+1]); %Net SW (W/m^2)
m3S_in_ds   = nan([size(mGlacMask),kEnd-kStart+1]); %Incoming, downscaled SW (W/m^2)
m3L_out     = nan([size(mGlacMask),kEnd-kStart+1]); %Emited LW (W/m^2)
m3L_in      = nan([size(mGlacMask),kEnd-kStart+1]); %Incoming LW (W/m^2)
m3Q_S       = nan([size(mGlacMask),kEnd-kStart+1]); %S (W/m^2)
m3Q_L       = nan([size(mGlacMask),kEnd-kStart+1]); %LE (W/m^2)
m3Q_P       = nan([size(mGlacMask),kEnd-kStart+1]); %Energy from precip (W/m^2)
m3T_s       = zeros([size(mGlacMask),kEnd-kStart+1]); %Surface temperature (C)
m3Precipitation       = zeros([size(mGlacMask),kEnd-kStart+1]); %Total precip, both phases (mm w.e./time)
m3Rain_Freeze_Amount  = zeros([size(mGlacMask),kEnd-kStart+1]); %Freezing rain (mm w.e./time?)
m3Q_net     = nan([size(mGlacMask),kEnd-kStart+1]); %Net surface energy (W/m^2)
m3Q_m       = nan([size(mGlacMask),kEnd-kStart+1]); %Melt energy (W/m^2)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for t = kStart:kEnd 
for t = [kStart+1:17544, kStart:kEnd] %SPINUP7NOV
    %% Extract timestep data
    mSW_in = vSW_in(t) .* mGlacMask;
    mSW_in = Calc_SW_in_dir_Sl_As1(mSW_in,mLat,mLong,mGlacAlt,mSlope,mAspect,stTime,mGlacMask,kTimeStep_days,sGCM,t);

    mU = vU(t) .* mGlacMask;

    %% Apply precipitation gradient up  (Precip increases with height and as a function of the amount of precip, unlike temp)

%   If no precipitation gradient used, uncomment:
%     mPrecip = (vP_d(t) .* mGlacMask);     

    hourly_avg = p1*mGlacAlt + p2; %@ alt., mm/hr
%   mPrecip = P @ HAR elev  + elev dif           * frac of average hourly P for yr * P lapse rate for yr (mm/m) .* mGlacMask;
    mPrecip = (vP_d(t) + (mGlacAlt - kAWS_Alt) .* (vP_d(t)./hourly_avg) * p1) .* mGlacMask;  %PLR in m P / m elev b/c vP_d in m/h
%       m     +          m                  m  * hr / mm       *        * mm / hr * m

    RAve.Pr_HAR(t)   = vP_d(t);
    RAve.Pr_glac(t)  = mean(mean((mPrecip(mGlacMask==1))));
    RAve.Pr_hravg(t) = mean(mean(hourly_avg));

% %  Apply precipitation correction (PC) to shift mean
if vP_d(t) > 0
    mPrecip = mPrecip + (13*GUI_Input.PC(g)/sum(hr_count)).* mGlacMask;
end  

    %% Calculate air pressure gradient with the hydrostatic equation (hPa/m)
    
    % Calculate day of the year 
    iDay        = datetime(stTime.year(t),stTime.month(t),stTime.day(t)); 
    iDay        = day(iDay,'dayofyear');
    iLapseRate  = vLapseRates_smooth(iDay);

    % Impose limits on lapse rates? (beta testing) - SHOULDN'T KICK IN BUT NOTE IT WAS LEFT UNCOMMENTED.
    if iLapseRate < 3
        iLapseRate = 3;
    elseif iLapseRate > 9.8
        iLapseRate = 9.8;
    end

    % Adjust air temperature for elevation (input lapse rate should be positive)
    mT_a = (vT_a(t) - (mGlacAlt - kAWS_Alt) * iLapseRate / 1000) .* mGlacMask;

    % Pressures
    mSVP_s = 6.11 * exp((2.5e6 / 461) .* (1 / 273.15 - 1 ./ (273.15 + mT_s))); % Saturation vapor pressure of the surface

    mVP_a = (6.11 * exp((2.5e6 / 461) * (1 / 273.15 - 1 ./ (273.15 + mT_a)))) .* vRH(t) / 100 .* mGlacMask; % Vapor pressure of the air 
%   BECAUSE: Vapor pressure of the air (hPa) from Clausius-Clapeyron, which relates saturation vapor pressure to temp (mbar = hPa)
%     vp_sat = (6.11 * exp((2.453e6 / 461) * (1 / 273.15 - 1 ./ (273.15 + mT_a))));
%     vp = vp_sat * rh
    
    % dP/dz = -P * g * (1 - e/p * (1-E))/(R_d*T_v)
    mPressGrad = -vP_a(t) .* 9.81 .* (1 - (mVP_a ./ vP_a(t)) * (1 - 0.622)) ./ (287.06 * (273.15 + mT_a));

    % Apply air pressure gradient (input pressure gradient should be negative)
    mP_a = (vP_a(t) + mPressGrad .* (mGlacAlt - kAWS_Alt)) .* mGlacMask;         

    %% Precipitation
    
    % Update snow clock (length of time since last snow event)
    mSnowClock = (mSnowClock + kTimeStep_days) .* mGlacMask;
    
    mTotalAccum(mT_a <= kRainThreshold) = (mTotalAccum(mT_a <= kRainThreshold) + mPrecip(mT_a <= kRainThreshold));
    mTotalAccum = mTotalAccum .* mGlacMask;

    % Reset snow clock anywhere it snowed more than specified precitation threshold                      
%     mSnowClock(mT_a <= kRainThreshold & mPrecip > kPrecipThresh) = 0;
    mSnowClock(mT_a <= kRainThreshold & mPrecip > 0) = 0;
    
    % If it snowed, add snow (only if T_a > kRainThreshold, i.e. snow, not rain) (m) - AG changed
    mSnowDepth(mT_a <= kRainThreshold) = (mSnowDepth(mT_a <= kRainThreshold) + mPrecip(mT_a <= kRainThreshold) ...
        * 1000 / kRho_snow) .* mGlacMask(mT_a <= kRainThreshold);

    %% Albedo
                     
    % Aging snow
    mAlpha_s    = kAlpha_fi + (kAlpha_fs - kAlpha_fi) .* exp(-ceil(mSnowClock) / kE_t);      % From Molg and Hardy, 2004 (Eq. 3)
    % Snow depth albedo
    mAlpha      = (mAlpha_s + (kAlpha_i - mAlpha_s) .* exp(-mSnowDepth / kE_d)).* mGlacMask; % From Molg and Hardy, 2004 (Eq. 2)

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
    mR_b = kG .* (mT_a - mT_s) .* (kZ - mZ_0m) ./ ((273.15 + mT_a) .* mU.^2);     % From Anderson and Others, 2010 (Eq. 17); same as Wagnon & al (2003)
    
    % ERIC'S Stability correction (only applied for sub-daily time steps)
%     if kTimeStep_days < 1
%         mStabCorr = (1 - 5 * mR_b).^2;
%         % Only use stability corrections if 0.2 > R_b > 0.01 (Braithwaite, 1995) <-- Rb defined differently AND 
%         % there should be a correction when R_b is negative!
%         mStabCorr(mR_b < 0.01)  = 1;                                                 
%         mStabCorr(mR_b > 0.2)   = 0;                                                  
%         mStabCorr(isinf(mStabCorr) == 1) = 1; %for low wind speeds, AG gathers
%     else
%         mStabCorr = 1*ones(size(mGlacMask)); % Remove stability correction for daily timestep
%     end

%Note: Wagnon & al (2003), cited by Molg & Hardy, stability correction (Phi functions) gave numerically non-sensical results

% From Anderson: "The stability correction (the final term in Equation (16)) is only applied when Rb > 0 and when U > 1, 
% the latter to prevent the bulk Richardson number from becoming unreasonably high when wind speeds are very low (Oke, 1987)."
mStabCorr = ones(size(mGlacMask)); 
       mStabCorr(mR_b > 0 & mU > 1 ) = (1 - 5 * mR_b(mR_b > 0 & mU > 1)).^2;                                              
       mStabCorr(mR_b < 0) = 1; 
       mStabCorr(mU < 1 ) = 1;      
       mStabCorr(isinf(mStabCorr)) = 1; 
       
    % k_H (Exchange coefficient for sensible heat)
    mK_H = kK_0^2 .* mStabCorr ./ (log(kZ ./ mZ_0m) .* (log(kZ ./ mZ_0T)));  % From Anderson and Others, 2010 (Eq. 16)

    % k_E (exchange coefficient for latent heat)
    mK_E = kK_0^2 .* mStabCorr ./ (log(kZ ./ mZ_0m) .* (log(kZ ./ mZ_0q)));          

    %% Sublimation vs evaporation
    
    % Latent heat of vaporization/fusion
    mL_vf = zeros(size(mGlacMask));
    mL_vf(mT_s  >= 0) = kL_v; % Evaporation
    mL_vf(mT_s < 0)   = kL_s; % Sublimation

    %% Layer thicknesses, densities, and diffusivities
    
    % Calculate snow thicknesses in the top layer (m)
    mSnowThickness_s = mSnowDepth .* mGlacMask;
    mSnowThickness_s(mSnowThickness_s > kLayerThick_s) = kLayerThick_s;
    % Calculate snow thickness under top 3 m [22 cm] (surface layer) (m)
    mSnowThickness_2 = (mSnowDepth - kLayerThick_s) .* mGlacMask;
    mSnowThickness_2(mSnowThickness_2 < 0) = 0;
    mSnowThickness_2(mSnowThickness_2 > kLayerThick_2) = kLayerThick_2;
    if mSnowThickness_2 > kLayerThick_2
        disp('Too much snow! Aborting.')
        returm
    end

%   Calculate amount of ice in both layers (m)
    mIceThickness_s = kLayerThick_s - mSnowThickness_s;
    mIceThickness_s(mIceThickness_s < 0) = 0;
    mIceThickness_2 = kLayerThick_2 - mSnowThickness_2;
    mIceThickness_2(mIceThickness_2 < 0) = 0; 
    
    % Calculate surface densities
    mRho_s = (mSnowThickness_s .* kRho_snow + mIceThickness_s * kRho_ice) / (kLayerThick_s);
    mRho_2 = (mSnowThickness_2 .* kRho_snow + mIceThickness_2 * kRho_ice) / (kLayerThick_2);
    
    % Calculate net surface diffusivity (m^2 s^-1)
%     mK_s = (mSnowThickness_s * kK_snow + mIceThickness_s * kK_ice) / (kLayerThick_s);
%     mK_2 = (mSnowThickness_2 * kK_snow + mIceThickness_2 * kK_ice) / (kLayerThick_2);
    
    %% Calculate effective emissivity
    
    kE_eff_a = vLW_in(t) / (kSigma * (273.15 + vT_a(t)).^4);
        
    %% Energy balance: Fluxes are positive towards surface

    % S(1-alpha) (W m^-2)
    mS_net = (mSW_in .* (1 - mAlpha)) .* mGlacMask;
    % L_out (W m^-2) - temp in K
    mL_out = -(kSigma * kE_i * (273.15 + mT_s).^4) .* mGlacMask;
    % L_in (W m^-2) - temp in K
    mL_in = kSigma .* kE_eff_a .* ((273.15 + mT_a) .^ 4) .* mGlacMask; % Braithwaite & Olesen (1990), though there are other methods
%     mL_in = kSigma .* ((273.15 + mT_a) .^ 4) .* (0.585 + 0.062 * mVP_a) .* mGlacMask; %Molg & Hardy - another option
    % Q_S (Sensible Heat, abbreviation "H" used above and in literature)
    mQ_S = (kRho_a_sl .* kC_a .* mK_H .* (mP_a ./ kP_a_sl) .* mU .* (mT_a - mT_s)) .* mGlacMask; 
    % Q_L (Latent Energy, abbreviation "E" used above and in literature)
    mQ_L = 0.622 * kRho_a_sl .* mK_E .* mU .* mL_vf .* (mVP_a - mSVP_s) ./ mP_a .* mGlacMask;
    % Q_P (W m^-2)
    mQ_P = 1000 * kC_w * mPrecip .* mT_a / kTimeStep_sec .* mGlacMask;   % From Hock 2005 
    % Precipitation heat flux only from rain (i.e. when T_a > kRainThreshold) 
    mQ_P(mT_a <= kRainThreshold) = 0;
    % Q_G (W m^-2)
        % ERIC'S: mQ_G = ( kK * mK_s .* (mT_2 - mT_s) ./ kLayerThick_s ) ./ kLayerThick_s .* kTimeStep_sec .* mGlacMask;     % From Paterson, 1994; The Physics of Glaciers, 3rd ed., pg. 206 (Eq. 5)
    mQ_G = (kK * (mT_2 - mT_s) ./ kLayerThick_s ) .* mGlacMask;  
    % Remove NaNs and stuff
    mQ_L(isnan(mQ_L)) = 0;
    % Net Energy (W m^-2) [AG] originally labeled "melt energy"
    mQ_net = (mS_net + mL_out + mL_in + mQ_S + mQ_L + mQ_P + mQ_G) .* mGlacMask;

%% Surface temperature: equations 11 and 12 of Arnold et al., 2006: roundabout approach
%  Calculate change in temperature of top 2 layers for this timestep (K)   
%     mDelta_T_s = (( mK_s .* (mT_2 - mT_s) ./ ...
%         kLayerThick_s ) ./ kLayerThick_s + (mQ_net - mQ_G) ./ (kC_i .* ...
%         mRho_s .* kLayerThick_s)) .* kTimeStep_sec .* mGlacMask;
%     mDelta_T_2 = (  ( mK_2 .* (mT_b - mT_2) ./ ...
%         kLayerThick_2 ) - ( mK_s .* (mT_2 - mT_s) )...
%         ./ kLayerThick_s  ) ./ kLayerThick_2 .* kTimeStep_sec .* mGlacMask;         
%     but in ^, mK_s ~= surf_diff = kK ./ (kC_i .* mRho_s)
    
%    mDelta_T_s = (   ( kK .* (mT_2 - mT_s) ./ kLayerThick_s ) ./ (kC_i .* mRho_s .* kLayerThick_s)...
%                 + (mQ_net - mQ_G)                            ./ (kC_i .* mRho_s .* kLayerThick_s)   ) ...
%                 .* kTimeStep_sec .* mGlacMask;    
    mDelta_T_s = (mQ_net ./ (kC_i .* mRho_s .* kLayerThick_s)) .* kTimeStep_sec .* mGlacMask;
%     mDelta_T_2 = (  kK .* (mT_b - mT_2) ./ (kC_i .* mRho_2 .* kLayerThick_2)  ...
%                     -  kK .* (mT_2 - mT_s) ./ (kC_i .* mRho_s .* kLayerThick_s)  ) ...
%                 ./ (kC_i .* mRho_s .* kLayerThick_2) .* kTimeStep_sec .* mGlacMask;     
    mDelta_T_2 = (  kK .* (mT_b - mT_2) ./ (kLayerThick_2)  ...
                    -  kK .* (mT_2 - mT_s) ./ (kLayerThick_s)  ) ...
                ./ (kC_i .* mRho_2 .* kLayerThick_2) .* kTimeStep_sec .* mGlacMask;  
            
    mT_s = mT_s + mDelta_T_s;
    mT_2 = mT_2 + mDelta_T_2;
    
%  mT_2 should never be above freezing
    if isempty(find((mT_2>0),1)) == 0
        [row, col]  = find(mT_2>0);
        fprintf('Surface Temp 2 exceeds freezing at row %d and col %d\n', row, col)
        return
    end
    
% Remove lines that don't conserve energy
    % mT_s(mT_s<-50) = -50; 
    % mT_s_Act_2(mT_s_Act_2 > 0) = 0; 
    %% Melt energy (W m^-2)
% ERIC'S version had Theoretical surface temperature (K) and Actual surface temperature (which can't be above zero) (K)
% The sole purpose of this "theoretical temperature" is the melt calculation.  
% Instead, compute melt & Ts from energy balance:
    mQ_m = zeros(size(mGlacMask));
%     mQ_m(mT_s_Theor_s > 0) = mT_s(mT_s > 0) * kC_i .* mRho_s(mT_s > 0) * kLayerThick_s / kTimeStep_sec;
    mQ_m(mT_s > 0) = mT_s(mT_s > 0) * kC_i .* ...
        mRho_s(mT_s > 0) * kLayerThick_s / kTimeStep_sec;
% % Surface temperature: these are equivalent within machine precision (checked then kept as == 0 to be consistent with rest of code)
    mT_s(mT_s > 0) = 0; %(C)
%     mT_s = mT_s - mQ_m * kTimeStep_sec ./ (kC_i .* mRho_s .* kLayerThick_s) .* mGlacMask;
%     mT_s (mT_s > -(1*10^(-12)) & mT_s < (1*10^(-12)) ) = 0; %reset to exactly zero... corrects for double precision

    %% Freezing of rain  
    % Rain that falls when the cold content of the surface layer is zero 
    % is assumed to runoff the glacier.
    
    mRain_Freeze_Energy = zeros(size(mGlacMask));
    % Freezing of rain (W m^-2) %AG changed
    mRain_Freeze_Energy(mT_a > kRainThreshold & mPrecip > 0) = mPrecip(mT_a > kRainThreshold & mPrecip > 0) ...
        * 1000 * kL_f / kTimeStep_sec;
    % Cold content of the surface layer (W m^-2)
    mColdContent_s = kC_i * mRho_s .* (kLayerThick_s * mRho_s / 1000) .* ...
        (0 - mT_s) / kTimeStep_sec;
    % Can't freeze more than cold content allows
    mRain_Freeze_Energy(mRain_Freeze_Energy > mColdContent_s) = ...
        mColdContent_s(mRain_Freeze_Energy > mColdContent_s);
    % Freezing water warms the surface layer
    mT_s = mT_s + ( mRain_Freeze_Energy ./ (kC_i .* mRho_s .* kLayerThick_s) .* kTimeStep_sec );
    % Add refreeze energy to Q_net (but not to Q_m)
    mQ_net = mQ_net + mRain_Freeze_Energy;
    % Calculate amount of frozen rain (m w.e.)
    mRain_Freeze_Amt = mRain_Freeze_Energy * kTimeStep_sec / (kL_f * 1000);
    % Add frozen rain as accumulation (m w.e.)
    mTotalAccum = (mTotalAccum + mRain_Freeze_Amt) .* mGlacMask;    %% Melting
  
    
    %% Snow melt; sublimation and evaporation
    % Snow melt (m)
    mSnowMelt = mQ_m * kTimeStep_sec / (kL_f * kRho_snow);

%     % Adjust snow depth for melt (m)
%     mSnowDepth = (mSnowDepth - mSnowMelt) .* mGlacMask;       
    
    % Sublimation (m)
    mSub = zeros(size(mGlacMask));
    mSub(mQ_L < 0 & mT_s < 0) = -mQ_L(mQ_L < 0 & mT_s < 0) .* kTimeStep_sec ./ (kL_s * kRho_snow) .* mGlacMask(mQ_L < 0 & mT_s < 0);  

    % Evaporation (m)
    mEvap = zeros(size(mGlacMask));
    mEvap(mQ_L < 0 & mT_s >= 0) = -mQ_L(mQ_L < 0 & mT_s >= 0) .* kTimeStep_sec ./ (kL_v * kRho_snow) .* mGlacMask(mQ_L < 0 & mT_s >= 0);  
    
%     % Adjust snow depth for sublimation (m)
%     mSnowDepth = (mSnowDepth - mSub) .* mGlacMask; 
    
    % Adjust snow depth for melt and sublimation (m)
    mSnowDepth = (mSnowDepth - mSnowMelt - mSub - mEvap) .* mGlacMask;   
    
    %% If no snow, melt ice
    
    % Find negative snow depths
    mXS_SnowMelt = zeros(size(mGlacMask));
    mXS_SnowMelt(mSnowDepth < 0) = mSnowDepth(mSnowDepth < 0);
    % Set negative snow depths to zero
    mSnowDepth(mSnowDepth < 0) = 0; %SHOULD BE: mSnowDepth(mSnowDepth < 0) + mXS_SnowMelt(mSnowDepth < 0)
    % Ice melt (m)
    mIceMelt = -mXS_SnowMelt * kRho_snow / kRho_ice;
    % Cumulative ice melt (m)
    mIceSurface = mIceSurface - mIceMelt;
    
    
    
    
    %% Net changes
    
    if t == kStart %SPINUP7NOV
        mTotalMelt  = zeros(size(mGlacMask));
        mTotalSnowMelt  = zeros(size(mGlacMask));
        mTotalIceMelt  = zeros(size(mGlacMask));
        mTotalSub   = zeros(size(mGlacMask));
        mTotalEvap  = zeros(size(mGlacMask));
        mTotalAccum = zeros(size(mGlacMask));
    end
        % Total average glacier melt (m w.eq.)
        mTotalMelt = mTotalMelt - (((mSnowMelt * kRho_snow / 1000) + ...
            (mIceMelt * kRho_ice / 1000)));
        mTotalSnowMelt = mTotalSnowMelt - (mSnowMelt * kRho_snow / 1000);
        mTotalIceMelt  = mTotalIceMelt  - (mIceMelt  * kRho_ice / 1000);
        % Total sublimation (m w.e.)
        mTotalSub = mTotalSub - (mSub .* kRho_snow / 1000);
        % Total evaporation (m w.e.)
        mTotalEvap = mTotalEvap - (mEvap .* kRho_snow / 1000);
    %     % Total mass balance (m w.eq.)
    %     mDeltaM_net = (mTotalMelt + mTotalAccum + mTotalSub);                   % Evaporation not included because melt is assumed to runoff currently   
     
    %% Average variables across glacier surface
    RAve.Q_L(t) = nanmean(mQ_L(mGlacMask==1));
    RAve.Q_S(t) = nanmean(mQ_S(mGlacMask==1));
    RAve.S_net(t) = nanmean(mS_net(mGlacMask==1));
    RAve.S_in(t) = vSW_in(t);
    RAve.S_in_ds(t) = nanmean(mSW_in(mGlacMask==1));
    RAve.L_in(t) = nanmean(mL_in(mGlacMask==1));
    RAve.L_out(t) = nanmean(mL_out(mGlacMask==1));
    RAve.Q_P(t) = nanmean(mQ_P(mGlacMask==1));
    RAve.Q_G(t) = nanmean(mQ_G(mGlacMask==1));
    RAve.Q_m(t) = nanmean(mQ_m(mGlacMask==1));
    RAve.Q_net(t) = nanmean(mQ_net(mGlacMask==1));   
    RAve.T_a(t) = nanmean(mT_a(mGlacMask==1));
    RAve.T_s(t) = nanmean(mT_s(mGlacMask==1));
    RAve.T_2(t) = nanmean(mT_2(mGlacMask==1));
    RAve.T_b(t)   = nanmean(mT_b(mGlacMask==1));
    RAve.Rain_Freeze_Energy(t) = nanmean(mRain_Freeze_Energy(mGlacMask==1));
    RAve.Rain_Freeze_Amount(t) = nanmean(mRain_Freeze_Amt(mGlacMask==1));
    RAve.Alpha(t) = nanmean(mAlpha(mGlacMask==1));
    RAve.Lapse_Rate(t) = iLapseRate; %to monitor the effect of the "beta testing" condition    
    RAve.TotalMelt(t) = nanmean(mTotalMelt(mGlacMask==1));
    RAve.TotalSub(t) = nanmean(mTotalSub(mGlacMask==1));
    RAve.TotalEvap(t) = nanmean(mTotalEvap(mGlacMask==1));
    RAve.TotalAccum(t) = nanmean(mTotalAccum(mGlacMask==1));
    
    RAve.SnowDepth(t) = nanmean(mSnowDepth(mGlacMask==1));
    RAve.SnowMelt(t)  = nanmean(mSnowMelt(mGlacMask==1));
    RAve.IceSurface(t) = nanmean(mIceSurface(mGlacMask==1));
    RAve.IceMelt(t) = nanmean(mIceMelt(mGlacMask==1));

    RAve.DeltaTs(t) = nanmean(nanmean(mQ_net(mGlacMask==1) ./ (kC_i .* mRho_s(mGlacMask==1)))); %\Delta T_s dz/dt
    RAve.DeltaT2(t) = nanmean(nanmean( (kK .* (mT_b(mGlacMask==1) - mT_2(mGlacMask==1))./(kLayerThick_2) - kK .* (mT_2(mGlacMask==1) - mT_s(mGlacMask==1)) ./ (kLayerThick_s)) ./ (kC_i .* mRho_2(mGlacMask==1)) )); %Delta T_2 dz/dt

    RAve.U(t) = nanmean(nanmean(mU(mGlacMask==1)));
    RAve.Prcp(t) = nanmean(nanmean(mPrecip(mGlacMask==1)));
    RAve.P_a(t) = nanmean(nanmean(mP_a(mGlacMask==1)));
    RAve.VP_a(t) = nanmean(nanmean(mVP_a(mGlacMask==1)));

    
% Save distributed variables for thickness tracking @ end of spin up
if t == 17544 
    if isnan(TH_i.Accum(1,1))
%     TH_i.SnowThickness_s = mSnowThickness_s;
%     TH_i.SnowThickness_2 = mSnowThickness_2;
%     TH_i.SnowDepth = mSnowDepth;
%     TH_i.SnowMelt  = mSnowMelt;
%     TH_i.IceSurface = mIceSurface;
%     TH_i.IceMelt = mIceMelt;
    TH_i.Melt = mTotalMelt;
    TH_i.SnowMelt = mTotalSnowMelt;
    TH_i.IceMelt = mTotalIceMelt;
    TH_i.Accum = mTotalAccum;
    TH_i.Sub = mTotalSub;
    TH_i.Evap = mTotalEvap;
    end
end
 % Save distributed variables for thickness tracking @ end of run
if t == kEnd
%     TH_f.SnowThickness_s = mSnowThickness_s;
%     TH_f.SnowThickness_2 = mSnowThickness_2;
%     TH_f.SnowDepth = mSnowDepth;
%     TH_f.SnowMelt  = mSnowMelt;
%     TH_f.IceSurface = mIceSurface;
%     TH_f.IceMelt = mIceMelt;
    TH_f.Melt = mTotalMelt;
    TH_f.SnowMelt = mTotalSnowMelt;
    TH_f.IceMelt = mTotalIceMelt;
    TH_f.Accum = mTotalAccum;
    TH_f.Sub = mTotalSub;
    TH_f.Evap = mTotalEvap;
end
       
    % Precipitation fraction that falls as snow
    if nansum(nansum(mPrecip)) ~= 0
        RAve.PrecipFrac(t) = nansum(nansum(mGlacMask(mT_a <= kRainThreshold & mPrecip > 0))) / nansum(nansum(mGlacMask(mPrecip > kRainThreshold))); %AG changed
    elseif nansum(nansum(mPrecip)) == 0
        RAve.PrecipFrac(t) = nan;
    end

%FOR SAVING FULLY DISTRIBUTED VARIABLES
if GUI_Input.threeDsave == 1
    m3TotalMelt(:,:,t) = mTotalMelt;
    m3TotalAccum(:,:,t) = mTotalAccum;
    m3Alpha(:,:,t) = mAlpha;
    m3S_net(:,:,t) = mS_net;
    m3S_in_ds(:,:,t) = mSW_in; 
    m3L_out(:,:,t) = mL_out;
    m3L_in(:,:,t) = mL_in;
    m3Q_S(:,:,t) = mQ_S;
    m3Q_L(:,:,t) = mQ_L;
    m3Q_P(:,:,t) = mQ_P;
    m3T_s(:,:,t) = mT_s;
    m3Precipitation(:,:,t) = mPrecip;
    m3Q_net(:,:,t) = mQ_net;
    m3Rain_Freeze_Amount(:,:,t) = mRain_Freeze_Amt;
    m3Q_m(:,:,t) = mQ_m;
end

% Plot for debugging
%     if rem(t,1)==0
%         figure(1);  subplot(411); plot(t,RAve.Q_G(t),'.'); ylabel('Q_G'); drawnow; hold on
%             subplot(412); plot(t,RAve.T_s(t),'.'); hold on; 
%                 plot(t,RAve.T_2(t),'.');
%                 plot(t,RAve.T_b(t),'.');
%                 drawnow; ylabel('Temps'); 
%             subplot(413); plot(t,nanmean(mQ_net ./ (kC_i .* mRho_s)),'.'); ylabel('\Delta T_s dz/dt'); drawnow; hold on
%             subplot(414); plot(t,nanmean( (kK .* (mT_b - mT_2)./(kLayerThick_2) - kK .* (mT_2 - mT_s) ./ (kLayerThick_s)) ./ (kC_i .* mRho_2) ),'.'); ylabel('\Delta T_2 dz/dt'); drawnow; hold on
%     (kK .* (mT_b - mT_2)./(kLayerThick_2) - kK .* (mT_2 - mT_s) ./ (kLayerThick_s)) ./ (kC_i .* mRho_2) 
%     end
    
%     if isnan(RAve.T_s(t))==1
%         keyboard
%     end
end    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Total Mass changes

% Total average melt (m w.eq.)
R2.TotalAveMelt = nanmean(mTotalMelt(mGlacMask == 1));
% Total accumulation (m w.eq.)
R2.TotalAveAccum = nanmean(mTotalAccum(mGlacMask == 1));
% Total average sublimation (m w.eq.)
R2.TotalAveSub = nanmean(mTotalSub(mGlacMask == 1));
% Total average evaporation (m w.eq.)
R2.TotalAveEvap = nanmean(mTotalEvap(mGlacMask == 1));

clear r
r = find (glID==y(g)); %index in geodedic MB spreadsheet
    R1.modelMB    = (R2.TotalAveMelt+R2.TotalAveAccum+R2.TotalAveSub)/13;
    R1.geodeticMB = geoMB(r);
    R1.sigma      = geoMB_sigma(r);

clearvars -except R1 R2 RAve R_Geo GUI_Input iGlacierNumber mTotalMelt mTotalSnowMelt mTotalIceMelt mTotalAccum mGlacMask mGlacAlt stTime ...
    mTotalSub mTotalEvap Data kAWS_Alt modelMB geodeticMB sigma TH_i TH_f y g

if GUI_Input.threeDsave == 1 % IF BIG
    save([GUI_Input.output_filename, num2str(iGlacierNumber),'_3D_90m_exploreSW.mat'],'-v7.3')
else
    save([GUI_Input.output_filename, num2str(iGlacierNumber),'_CALcheck.mat'])
end

end
