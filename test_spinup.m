clear
GUI_Input.glacier_number = 1421894; %first glacier in AKGJ y(g=121)
GUI_Input.sGCM          = 'HAR';
GUI_Input.HAR_temp_res  = 'hourly'; % 'daily' or 'hourly'
addpath('/uufs/chpc.utah.edu/common/home/u6027899/sebm');

%% Load files
mGlacNum    = load('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_90m_Numbers.mat','mGlacNum');
mGlacNum    = mGlacNum.mGlacNum;
mDebrisMask = load('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_90m_Debris_Mask.mat','mDebrisMask'); 
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


sDataDirectory = '/uufs/chpc.utah.edu/common/home/steenburgh-group6/Eric/HAR_Data';
sTempRes = GUI_Input.HAR_temp_res;

[ vSW_in, vU, vT_a, vP_d, vP_a, vRH, ~, vLW_in, stTime, Data  ] = ExtractHAR1( kLat, kLong, sTempRes, sDataDirectory );

% 2001 indices: 
    t1 = 8785;
    t2 = 17544;
% cut vectors:
cSW_in = vSW_in(t1:t2);
cU   = vU(t1:t2);
cT_a = vT_a(t1:t2);
cP_d = vP_d(t1:t2);
cP_a = vP_a(t1:t2);
cRH  = vRH(t1:t2);
cLW_in = vLW_in(t1:t2);

% full vectors:
fSW_in  = repmat(cSW_in',1,13);   
fU      = repmat(cU',1,13); 
fT_a    = repmat(cT_a',1,13); 
fP_d    = repmat(cP_d',1,13); 
fP_a    = repmat(cP_a',1,13); 
fRH     = repmat(cRH',1,13); 
fLW_in  = repmat(cLW_in',1,13); 

% insert leap years
ly_2004 = 36457; %Feb 28 1am 2004
ly_2008 = 71521; %Feb 28 1am 2008
ly_2012 = 106585; %Feb 28 1am 2012

suSW_in = [nan(1,8784) fSW_in(1:ly_2004+23) fSW_in(ly_2004:ly_2004+23) fSW_in(ly_2004+24:ly_2008+23) ...
                                            fSW_in(ly_2008:ly_2008+23) fSW_in(ly_2008+24:ly_2012+23) ...
                                            fSW_in(ly_2012:ly_2012+23) fSW_in(ly_2012+24:end) nan(1,8760)]; %131496 
suU = [nan(1,8784) fU(1:ly_2004+23) fU(ly_2004:ly_2004+23) fU(ly_2004+24:ly_2008+23) ...
                                    fU(ly_2008:ly_2008+23) fU(ly_2008+24:ly_2012+23) ...
                                    fU(ly_2012:ly_2012+23) fU(ly_2012+24:end) nan(1,8760)]; %131496                                         
suT_a = [nan(1,8784) fT_a(1:ly_2004+23) fT_a(ly_2004:ly_2004+23) fT_a(ly_2004+24:ly_2008+23) ...
                                        fT_a(ly_2008:ly_2008+23) fT_a(ly_2008+24:ly_2012+23) ...
                                        fT_a(ly_2012:ly_2012+23) fT_a(ly_2012+24:end) nan(1,8760)]; %131496                                         
suP_d = [nan(1,8784) fP_d(1:ly_2004+23) fP_d(ly_2004:ly_2004+23) fP_d(ly_2004+24:ly_2008+23) ...
                                        fP_d(ly_2008:ly_2008+23) fP_d(ly_2008+24:ly_2012+23) ...
                                        fP_d(ly_2012:ly_2012+23) fP_d(ly_2012+24:end) nan(1,8760)]; %131496                                         
suP_a = [nan(1,8784) fP_a(1:ly_2004+23) fP_a(ly_2004:ly_2004+23) fP_a(ly_2004+24:ly_2008+23) ...
                                        fP_a(ly_2008:ly_2008+23) fP_a(ly_2008+24:ly_2012+23) ...
                                        fP_a(ly_2012:ly_2012+23) fP_a(ly_2012+24:end) nan(1,8760)]; %131496     
suRH = [nan(1,8784) fRH(1:ly_2004+23) fRH(ly_2004:ly_2004+23) fRH(ly_2004+24:ly_2008+23) ...
                                      fRH(ly_2008:ly_2008+23) fRH(ly_2008+24:ly_2012+23) ...
                                      fRH(ly_2012:ly_2012+23) fRH(ly_2012+24:end) nan(1,8760)]; %131496                                     
suLW_in = [nan(1,8784) fLW_in(1:ly_2004+23) fLW_in(ly_2004:ly_2004+23) fLW_in(ly_2004+24:ly_2008+23) ...
                                            fLW_in(ly_2008:ly_2008+23) fLW_in(ly_2008+24:ly_2012+23) ...
                                            fLW_in(ly_2012:ly_2012+23) fLW_in(ly_2012+24:end) nan(1,8760)]; %131496      
return                                   
save ('HAR_spinup.mat','suSW_in','suU','suT_a','suP_d','suP_a','suRH','suLW_in')