A = csvread('combination.csv'); 
glac = A(:,1);

% Border (in matrix indices) necessary for sky view factor calculation
kBorder         = 180;

mGlacNum    = load('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_30m_Numbers.mat','mGlacNum');
mGlacNum    = mGlacNum.mGlacNum;
mDebrisMask = load('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_30m_Debris_Mask.mat','mDebrisMask'); 
mDebrisMask = mDebrisMask.mDebrisMask; 

for j=1:length(glac)
iGlacierNumber  = glac(j);

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

glac_num(j)      = iGlacierNumber;
DEM_min_elev(j) = nanmin(mGlacAlt(mGlacMask==1));
DEM_max_elev(j) = nanmax(mGlacAlt(mGlacMask==1));
HAR_elev(j) = kAWS_Alt;
j
clear mGlacNum mDEM mMask mGlacAlt mGlacMask
end
save('elevations.mat','DEM_min_elev','DEM_max_elev','HAR_elev','glac_num')
