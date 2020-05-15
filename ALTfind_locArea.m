%Get glacier locations and areas to evaluate by sub-basin
clear

M = csvread('shean_UIB.csv');

mGlacNum = load('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_30m_Numbers.mat','mGlacNum');
mGlacNum = mGlacNum.mGlacNum;
mDEM    = ncread('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_30m_DEM.nc','mGlacAlt');
mMask   = ncread('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_30m_Glacier_Mask.nc','mGlacMask');

kColLn = size(mGlacNum,1); %ACTUALLY NUMBER OF ROWS, EJ: 15600; 
kRowLn = size(mGlacNum,2); %ACTYALLY NUMBER OF COL,  EJ: 45600; 
%% Glacier to extract
for j = 1:length(M) %ROWS
    iGlacierNumber  = M(j,1);

% Find indices of glacier edges
[row, col]  = find(mGlacNum==iGlacierNumber); 
iM_i        = min(row);
iM_f        = max(row);
iN_i        = min(col);
iN_f        = max(col);
mGlacNum_Extract    = mGlacNum(iM_i:iM_f,iN_i:iN_f); %+1 removed

% Extract DEMs and Masks for given glacier
% mDEM    = ncread('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_30m_DEM.nc','mGlacAlt',[iM_i, iN_i],[iM_f-iM_i+1, iN_f-iN_i+1]); %+1 needed for ncread 'count'
% mMask   = ncread('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_30m_Glacier_Mask.nc','mGlacMask',[iM_i, iN_i],[iM_f-iM_i+1, iN_f-iN_i+1]); %+1 needed for ncread 'count'
mDEM_Extract    = mDEM(iM_i:iM_f,iN_i:iN_f); 
mMask_Extract   = mMask(iM_i:iM_f,iN_i:iN_f); 

% Remove other glaciers from the glacier mask
mMask_Extract(mGlacNum_Extract ~= iGlacierNumber) = 0;
% mMask(mGlacNum_Extract ~= iGlacierNumber) = 0;

mGlacAlt    = double(mDEM_Extract);
mGlacMask   = mMask_Extract;
% mGlacAlt    = double(mDEM);
% mGlacMask   = mMask;

% Make lat/lon matrices
vLat    = (linspace(28.999861, 37.117361,kColLn))'; %28.999861 37.117361
mLat    = flipud(repmat(vLat,1,kRowLn));
mLat    = mLat(iM_i:iM_f,iN_i:iN_f); %+1 removed
vLong   = linspace(67.618750, 82.500694,kRowLn); %67.6187499 82.50069444
mLong   = repmat(vLong,kColLn,1);
mLong   = mLong(iM_i:iM_f,iN_i:iN_f); %+1 removed
% kLat    = nanmean(mLat(mMask==1)); % Find average latitude of glacier
% kLong   = nanmean(mLong(mMask==1)); % Find average longitude of glacier
%lat long where alt lowest
kLat    = mLat(mGlacAlt==min(min(mGlacAlt))); % Find average latitude of glacier
kLong   = mLong(mGlacAlt==min(min(mGlacAlt))); % Find average longitude of glacier

lat(j) = kLat(1);
lon(j) = kLong(1);
area(j)= M(j,13); %m^2

if rem(j,100)==0
    j
end

% lat(j), lon(j), area(j), keyboard
% clearvars -except j M mGlacNum mDEM mMask kColLn kRowLn

end

save low_alt_location_area.mat
