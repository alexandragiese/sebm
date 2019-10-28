function [ vLapseRates ] = CalcLapseRates1(  )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here




% kLapseRates = GUI_Input.lapse_rate; 
addpath(genpath('/uufs/chpc.utah.edu/common/home/steenburgh-group6/Eric/HAR_Data'))
mLat_HAR = double(ncread('har_d10km_static_hgt.nc', 'lat'));
mLat_HAR = flipud(mLat_HAR');
mLong_HAR = double(ncread('har_d10km_static_hgt.nc','lon'));
mLong_HAR = flipud(mLong_HAR');
mLat_HAR = mLat_HAR(:);
mLong_HAR = mLong_HAR(:);
[~, iLatLongIdx] = min((mLat_HAR - kLat).^2 + (mLong_HAR - kLong).^2);
[m_Idx, n_Idx] = ind2sub([180, 270], iLatLongIdx);

m3LapseRates = load('m3LapseRates_HAR.mat','m3LapseRates');
m3LapseRates = m3LapseRates.m3LapseRates;
vLapseRates = squeeze(m3LapseRates(m_Idx,n_Idx,:));

% Check for NaNs
if isnan(vLapseRates(1)) == 1
   m3LapseRates_temp = m3LapseRates(m_Idx-1:m_Idx+1,n_Idx-1:n_Idx+1,:);
   vLapseRates = squeeze(nanmean(nanmean(m3LapseRates_temp,1),2));
    if isnan(vLapseRates(1)) == 1
        m3LapseRates_temp = m3LapseRates(m_Idx-2:m_Idx+2,n_Idx-2:n_Idx+2,:);
        vLapseRates = squeeze(nanmean(nanmean(m3LapseRates_temp,1),2));
        if isnan(vLapseRates(1)) == 1
            m3LapseRates_temp = m3LapseRates(m_Idx-3:m_Idx+3,n_Idx-3:n_Idx+3,:);
            vLapseRates = squeeze(nanmean(nanmean(m3LapseRates_temp,1),2));
        end
    end
   clearvars m3LapseRates_temp     
end







































end

