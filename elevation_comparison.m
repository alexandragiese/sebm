addpath('/uufs/chpc.utah.edu/common/home/u6027899/sebm')
A = csvread('combination.csv'); 
glac = A(:,1);

% Border (in matrix indices) necessary for sky view factor calculation
kBorder         = 180;

mGlacNum    = load('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_30m_Numbers.mat','mGlacNum');
mGlacNum    = mGlacNum.mGlacNum;
mDebrisMask = load('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_30m_Debris_Mask.mat','mDebrisMask'); 
mDebrisMask = mDebrisMask.mDebrisMask; 
kColLn = size(mGlacNum,1); %ACTUALLY NUMBER OF ROWS, EJ: 15600; 
kRowLn = size(mGlacNum,2); %ACTYALLY NUMBER OF COL,  EJ: 45600; 
c=0;
for j=8990:length(glac)
iGlacierNumber  = glac(j);

% Find indices of glacier edges
[row, col]  = find(mGlacNum==iGlacierNumber); 
iM_i        = min(row) - kBorder; if iM_i < 1, iM_i = 1; c = c+1; indx_toolow(c) = iM_i; end
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

% Make lat/lon matrices
vLat    = (linspace(28.999861, 37.117361,kColLn))'; %28.999861 37.117361
mLat    = flipud(repmat(vLat,1,kRowLn));
mLat    = mLat(iM_i:iM_f,iN_i:iN_f); %+1 removed
vLong   = linspace(67.618750, 82.500694,kRowLn); %67.6187499 82.50069444
mLong   = repmat(vLong,kColLn,1);
mLong   = mLong(iM_i:iM_f,iN_i:iN_f); %+1 removed

kLat    = nanmean(mLat(mMask==1)); % Find average latitude of glacier
kLong   = nanmean(mLong(mMask==1)); % Find average longitude of glacier

sDataDirectory = '/uufs/chpc.utah.edu/common/home/steenburgh-group6/Eric/HAR_Data';
sTempRes  = 'hourly'; % 'daily' or 'hourly'
[ vSW_in, vU, vT_a, vP_d, vP_a, vRH, ~, vLW_in, stTime, Data  ] = ExtractHAR1( kLat, kLong, sTempRes, sDataDirectory );

glac_num(j)     = iGlacierNumber;
DEM_min_elev(j) = nanmin(mGlacAlt(mGlacMask==1));
DEM_max_elev(j) = nanmax(mGlacAlt(mGlacMask==1));
HAR_elev(j) = Data.AWS_Alt;
j
clear mDEM mMask mGlacAlt mGlacMask kLat kLong vSW_in vU vT_a vP_d vP_a vRH vLW_in stTime Data
end
save('elevations3.mat','DEM_min_elev','DEM_max_elev','HAR_elev','glac_num') %through 7675, 8989, end
