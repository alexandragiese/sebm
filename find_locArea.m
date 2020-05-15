%Get glacier locations and areas to evaluate by sub-basin
clear
M = csvread('shean_UIB.csv');

mGlacNum = load('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_30m_Numbers.mat','mGlacNum');
mGlacNum = mGlacNum.mGlacNum;
kColLn = size(mGlacNum,1); %ACTUALLY NUMBER OF ROWS, EJ: 15600; 
kRowLn = size(mGlacNum,2); %ACTYALLY NUMBER OF COL,  EJ: 45600; 

%% Glacier to extract
for j = 1:length(M) 
    iGlacierNumber  = M(j,1);

% Border (in matrix indices) necessary for sky view factor calculation
kBorder         = 0;

% Find indices of glacier edges
[row, col]  = find(mGlacNum==iGlacierNumber); 
iM_i        = min(row) - kBorder;
iM_f        = max(row) + kBorder;
iN_i        = min(col) - kBorder;
iN_f        = max(col) + kBorder;
mGlacNum_Extract    = mGlacNum(iM_i:iM_f,iN_i:iN_f); %+1 removed
% mDebrisMask_Extract = mDebrisMask(min(row):max(row),min(col):max(col)); 

% Extract DEMs and Masks for given glacier
% mDEM    = ncread('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_30m_DEM.nc','mGlacAlt',[iM_i, iN_i],[iM_f-iM_i+1, iN_f-iN_i+1]); %+1 needed for ncread 'count'
mMask   = ncread('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_30m_Glacier_Mask.nc','mGlacMask',[iM_i, iN_i],[iM_f-iM_i+1, iN_f-iN_i+1]); %+1 needed for ncread 'count'

% Remove other glaciers from the glacier mask
mMask(mGlacNum_Extract ~= iGlacierNumber) = 0;

% mGlacAlt    = double(mDEM);
mGlacMask   = mMask;
% mDebMask    = mDebrisMask_Extract;

% Make lat/lon matrices
vLat    = (linspace(28.999861, 37.117361,kColLn))'; %28.999861 37.117361
mLat    = flipud(repmat(vLat,1,kRowLn));
mLat    = mLat(iM_i:iM_f,iN_i:iN_f); %+1 removed
vLong   = linspace(67.618750, 82.500694,kRowLn); %67.6187499 82.50069444
mLong   = repmat(vLong,kColLn,1);
mLong   = mLong(iM_i:iM_f,iN_i:iN_f); %+1 removed

kLat    = nanmean(mLat(mMask==1)); % Find average latitude of glacier
kLong   = nanmean(mLong(mMask==1)); % Find average longitude of glacier

lat(j) = kLat;
lon(j) = kLong;
area(j)= M(j,13); %m^2

clearvars -except j M mGlacNum kColLn kRowLn lat lon area
if rem(j,100)==0
    disp(j)
end
if j==627, keyboard, end
end
save location_area.mat

