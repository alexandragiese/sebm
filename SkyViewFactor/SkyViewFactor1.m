function [mSVF] = SkyViewFactor1(mGlacMask,mGlacAlt,mLat,mLong,kRadius)
% Calculate the sky view factor for each location on a given glacier 

[mLong_UTM,mLat_UTM,~] = ll2utm(mLat,mLong);

vGlacMaskIdx = find(mGlacMask==1);

% Define circle stuff
% kRadius = 100; % Gridcells (1 gridcell = 90m)
kNumAnglesPerPoint=36;
vAngle=linspace(-pi,pi,kNumAnglesPerPoint+1);

mSVF = zeros(size(mGlacMask));

for t = 1:nansum(mGlacMask(:))

    iLat_UTM = mLat_UTM(vGlacMaskIdx(t));
    iLong_UTM = mLong_UTM(vGlacMaskIdx(t));
    iGlacAlt = mGlacAlt(vGlacMaskIdx(t));
    
    R = sqrt((mLong_UTM-iLong_UTM).^2 + (mLat_UTM-iLat_UTM).^2);
    Z = (mGlacAlt-iGlacAlt);

    mTangents = Z./R;

    %% Create cirle of radius kRadius with kN points in the circle
    
    [iY,iX] = ind2sub(size(mGlacMask),vGlacMaskIdx(t));
    
    vX_idx = round(kRadius * cos(vAngle) + iX);
    vY_idx = round(kRadius * sin(vAngle) + iY);
    
    mCircle = zeros(size(mGlacMask));
      
    for g = 1:kNumAnglesPerPoint
    
        mCircle(vY_idx(g),vX_idx(g)) = 1;

    end

    %% Bresenham 
    
    vCircle_idx = find(mCircle == 1);
    [vCircle_idx_Y, vCircle_idx_X] = ind2sub(size(mGlacMask), vCircle_idx); 
    
    vMaxTangents = zeros(kNumAnglesPerPoint,1);
    
    for g = 1:kNumAnglesPerPoint
        
        [vX, vY] = bresenham(iX, iY, vCircle_idx_X(g), vCircle_idx_Y(g));
        
        vTangents_idx = sub2ind(size(mGlacMask),vX,vY);
        vMaxTangents(g) = max(mTangents(vTangents_idx));

    end
    
    vHorizonAngles = atand(vMaxTangents);
    vHorizonAngles(vHorizonAngles < 0) = 0;
    iAveHorizonAngle = nanmean(vHorizonAngles);

    mSVF(vGlacMaskIdx(t)) = (cosd(iAveHorizonAngle))^2;

end




















% 
% % Used to define the circle of radius kRadius around the point of interest
% x = 1:size(mGlacMask,2);
% y = 1:size(mGlacMask,1);
% [mX, mY] = meshgrid(x,y);
% 
%     [iY,iX] = ind2sub(size(mGlacMask),vGlacMaskIdx(t));
% 
%     mCircle = zeros(size(mGlacMask));
%     mCircle(((mX-iX).^2+(mY-iY).^2) < kRadius^2) = 1;
%     mBoundaryIdx = bwboundaries(mCircle,8,'noholes');
%     mBoundaryIdx = mBoundaryIdx{:};
%     mCircleOutline = zeros(size(mGlacMask));

%     for g = 1:length(mBoundaryIdx)
%         mCircleOutline(mBoundaryIdx(g,1),mBoundaryIdx(g,2)) = 1;
%     end





    
%     [mVis_temp,~] = viewshed(mGlacAlt,R_Geo,iLat,iLong);
%     [B,L] = bwboundaries(mVis_temp,8,'noholes');
%     
%     mHorizonMask = mHorizonPts;
%     mHorizonMask(mHorizonMask ~= 0) = 1;
%     
    
%     for g = 1:nansum(mHorizonMask(:))
%         
%         
%         
%         
%         
%     end











































end

