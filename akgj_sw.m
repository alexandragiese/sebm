addpath('/uufs/chpc.utah.edu/common/home/u6027899/sebm');

SW_akgj = nan(30,131496);
load 30gl_eachbasin.mat
    foo = ~isnan(glac_nums(:));
    y = glac_nums(foo); %glacier_number = y(121:150); 
    
for g = 1:30
    iGlacierNumber  = y(g+120); 
    
%% Load files
mGlacNum    = load('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_90m_Numbers.mat','mGlacNum');
mGlacNum    = mGlacNum.mGlacNum;

mDebrisMask = load('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_90m_Debris_Mask.mat','mDebrisMask'); 
mDebrisMask = mDebrisMask.mDebrisMask; 

%% Specifics for larger region of interest

kColLn = size(mGlacNum,1); 
kRowLn = size(mGlacNum,2); 

%% Glacier to extract


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
mDEM    = ncread('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_90m_DEM.nc','mGlacAlt',[iM_i, iN_i],[iM_f-iM_i+1, iN_f-iN_i+1]); %+1 needed for ncread 'count'
mMask   = ncread('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_90m_Glacier_Mask.nc','mGlacMask',[iM_i, iN_i],[iM_f-iM_i+1, iN_f-iN_i+1]); %+1 needed for ncread 'count'

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
kRadius                 = 98;
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

sGCM          = 'HAR';
sDataDirectory = '/uufs/chpc.utah.edu/common/home/steenburgh-group6/Eric/HAR_Data';
sTempRes  = 'hourly'; 

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
SW_akgj(g,:) = vSW_in;
g
end

save('SW_30gl_inAKGJ.mat','SW_akgj')
load SW_30gl_inAKGJ.mat
% SW_akgj((SW_akgj)==0)=nan;
figure; histogram(nanmean(SW_akgj,2)); xlabel('Mean SW (W m^{-2})'); print -dpng meanSWakgj %std = 12.4417
figure; histogram(nanstd(SW_akgj'));   xlabel('St. dev. SW (W m^{-2})'); print -dpng stdevSWakgj

load 30gl_eachbasin.mat
    foo = ~isnan(glac_nums(:));
    y = glac_nums(foo);   
% Make .mat file containing downscaled SW_in---------
k = 121;
for j = 1:30
    load (['../SEBM_output/Glacier_Number_',num2str(y(k)),'_res90_getSW.mat'],'RAve') 
    SW_in_ds(j,:) = RAve.S_in_ds;
    clear RAve
    k = k+1;
end

% save ('SW_downscaled.mat','SW_in_ds')
figure; histogram(nanmean(SW_in_ds,2)); xlabel('Mean SW_{ds} (W m^{-2})'); print -dpng meanSWakgj_ds %std = 9.0990
figure; histogram(nanstd(SW_in_ds'));   xlabel('St. dev. SW_{ds} (W m^{-2})'); print -dpng stdevSWakgj_ds