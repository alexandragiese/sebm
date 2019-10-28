function [ SW_in, WindSpeed, Temp, Precip, Pressure, RelHumidity, LW_in, Data ] = ExtractWRF1( Lat, Long, sFileDir  )
%UNTITLED2 Summary of this function goes here
%   Extract climate data for a single 10km x 10km grid point for use with
%   surface and energy mass balance model
%
%
% http://www2.mmm.ucar.edu/wrf/users/docs/user_guide/users_guide_chap5.html
%
% http://www.meteo.unican.es/wiki/cordexwrf/OutputVariables

addpath(genpath(sFileDir))
cd(sFileDir)

% mLat = double(ncread('har_d10km_h_2d_t2_2001.nc', 'lat'));
% mLong = double(ncread('har_d10km_h_2d_t2_2001.nc','lon'));

mLat = mLat(:);
mLong = mLong(:);

[~, iLatLongIdx] = min((mLat - Lat).^2 + (mLong - Long).^2);
% [m_Idx, n_Idx] = ind2sub([270, 180], iLatLongIdx);

%% Add WRF data to the path

% '2000', '2000_mp7', '2001','2002', or '2003'
MicroPhysics = '2000';


rmpath('/uufs/chpc.utah.edu/common/home/iutah1-0/adamk/HIMAT/WRF_output/2000')
rmpath('/uufs/chpc.utah.edu/common/home/iutah1-0/adamk/HIMAT/WRF_output/2000_mp7')

if strcmp(MicroPhysics,'2000_mp7')
    addpath('/uufs/chpc.utah.edu/common/home/iutah1-0/adamk/HIMAT/WRF_output/2000_mp7')
elseif strcmp(MicroPhysics,'2000')
    addpath('/uufs/chpc.utah.edu/common/home/iutah1-0/adamk/HIMAT/WRF_output/2000/output')
elseif strcmp(MicroPhysics,'2001')
    addpath('/uufs/chpc.utah.edu/common/home/iutah1-0/adamk/HIMAT/WRF_output/2001/output')
elseif strcmp(MicroPhysics,'2002')
    addpath('/uufs/chpc.utah.edu/common/home/iutah1-0/adamk/HIMAT/WRF_output/2002/output')
elseif strcmp(MicroPhysics,'2003')
    addpath('/uufs/chpc.utah.edu/common/home/iutah1-0/adamk/HIMAT/WRF_output/2003/output')
end

%% Find matrix indices best associated with given lat/long coordinates

% mLat = single(ncread('wrfout_d03_2000-01-01_00:00:00', 'XLAT',[1 1 1],[525 402 1],[1 1 1]));
% mLong = single(ncread('wrfout_d03_2000-01-01_00:00:00','XLONG',[1 1 1],[525 402 1],[1 1 1]));
% 
% mLat = mLat(:);
% mLong = mLong(:);
% 
% [~, iLatLongIdx] = min((mLat - Lat).^2 + (mLong - Long).^2);
% [m_Idx, n_Idx] = ind2sub([525, 402], iLatLongIdx);

%% Grid Elevation

mAWS_Alt = double(ncread('wrfout_d03_2000-01-01_00:00:00','HGT',[1 1 1],[525 402 1],[1 1 1]));
Data.AWS_Alt = mAWS_Alt(m_Idx, n_Idx);

%% Temperature at 2m above the surface (C)

Temp1 = squeeze(ncread('wrfout_d03_2000-01-01_00:00:00','T2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Temp2 = squeeze(ncread('wrfout_d03_2000-01-31_13:00:00','T2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Temp3 = squeeze(ncread('wrfout_d03_2000-03-02_01:00:00','T2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Temp4 = squeeze(ncread('wrfout_d03_2000-04-01_13:00:00','T2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Temp5 = squeeze(ncread('wrfout_d03_2000-05-02_01:00:00','T2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Temp6 = squeeze(ncread('wrfout_d03_2000-06-01_13:00:00','T2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Temp7 = squeeze(ncread('wrfout_d03_2000-07-02_01:00:00','T2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Temp8 = squeeze(ncread('wrfout_d03_2000-08-01_13:00:00','T2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Temp9 = squeeze(ncread('wrfout_d03_2000-09-01_01:00:00','T2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Temp10 = squeeze(ncread('wrfout_d03_2000-10-01_13:00:00','T2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Temp11 = squeeze(ncread('wrfout_d03_2000-11-01_01:00:00','T2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Temp12 = squeeze(ncread('wrfout_d03_2000-12-01_13:00:00','T2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));

Temp = [Temp1; Temp2; Temp3; Temp4; Temp5; Temp6; Temp7; Temp8; Temp9;...
    Temp10; Temp11; Temp12];

clearvars Temp1 Temp2 Temp3 Temp4 Temp5 Temp6 Temp7 Temp8 Temp9 Temp10 ...
    Temp11 Temp12

% Convert Kelvin to Celsius
Temp = Temp - 273.15;

%% SW_in (W m^-2)

SW_in1 = squeeze(ncread('wrfout_d03_2000-01-01_00:00:00','SWDOWN',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
SW_in2 = squeeze(ncread('wrfout_d03_2000-01-31_13:00:00','SWDOWN',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
SW_in3 = squeeze(ncread('wrfout_d03_2000-03-02_01:00:00','SWDOWN',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
SW_in4 = squeeze(ncread('wrfout_d03_2000-04-01_13:00:00','SWDOWN',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
SW_in5 = squeeze(ncread('wrfout_d03_2000-05-02_01:00:00','SWDOWN',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
SW_in6 = squeeze(ncread('wrfout_d03_2000-06-01_13:00:00','SWDOWN',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
SW_in7 = squeeze(ncread('wrfout_d03_2000-07-02_01:00:00','SWDOWN',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
SW_in8 = squeeze(ncread('wrfout_d03_2000-08-01_13:00:00','SWDOWN',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
SW_in9 = squeeze(ncread('wrfout_d03_2000-09-01_01:00:00','SWDOWN',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
SW_in10 = squeeze(ncread('wrfout_d03_2000-10-01_13:00:00','SWDOWN',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
SW_in11 = squeeze(ncread('wrfout_d03_2000-11-01_01:00:00','SWDOWN',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
SW_in12 = squeeze(ncread('wrfout_d03_2000-12-01_13:00:00','SWDOWN',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));

SW_in = [SW_in1; SW_in2; SW_in3; SW_in4; SW_in5; SW_in6; SW_in7; SW_in8;...
    SW_in9; SW_in10; SW_in11; SW_in12];

clearvars SW_in1 SW_in2 SW_in3 SW_in4 SW_in5 SW_in6 SW_in7 SW_in8 SW_in9...
    SW_in10 SW_in11 SW_in12

%% WindSpeed @ 10m (m s^-1)

U1 = squeeze(ncread('wrfout_d03_2000-01-01_00:00:00','U10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
U2 = squeeze(ncread('wrfout_d03_2000-01-31_13:00:00','U10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
U3 = squeeze(ncread('wrfout_d03_2000-03-02_01:00:00','U10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
U4 = squeeze(ncread('wrfout_d03_2000-04-01_13:00:00','U10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
U5 = squeeze(ncread('wrfout_d03_2000-05-02_01:00:00','U10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
U6 = squeeze(ncread('wrfout_d03_2000-06-01_13:00:00','U10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
U7 = squeeze(ncread('wrfout_d03_2000-07-02_01:00:00','U10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
U8 = squeeze(ncread('wrfout_d03_2000-08-01_13:00:00','U10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
U9 = squeeze(ncread('wrfout_d03_2000-09-01_01:00:00','U10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
U10 = squeeze(ncread('wrfout_d03_2000-10-01_13:00:00','U10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
U11 = squeeze(ncread('wrfout_d03_2000-11-01_01:00:00','U10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
U12 = squeeze(ncread('wrfout_d03_2000-12-01_13:00:00','U10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));

WindSpeed_U = [U1; U2; U3; U4; U5; U6; U7; U8; U9; U10; U11; U12];

V1 = squeeze(ncread('wrfout_d03_2000-01-01_00:00:00','V10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
V2 = squeeze(ncread('wrfout_d03_2000-01-31_13:00:00','V10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
V3 = squeeze(ncread('wrfout_d03_2000-03-02_01:00:00','V10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
V4 = squeeze(ncread('wrfout_d03_2000-04-01_13:00:00','V10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
V5 = squeeze(ncread('wrfout_d03_2000-05-02_01:00:00','V10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
V6 = squeeze(ncread('wrfout_d03_2000-06-01_13:00:00','V10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
V7 = squeeze(ncread('wrfout_d03_2000-07-02_01:00:00','V10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
V8 = squeeze(ncread('wrfout_d03_2000-08-01_13:00:00','V10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
V9 = squeeze(ncread('wrfout_d03_2000-09-01_01:00:00','V10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
V10 = squeeze(ncread('wrfout_d03_2000-10-01_13:00:00','V10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
V11 = squeeze(ncread('wrfout_d03_2000-11-01_01:00:00','V10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
V12 = squeeze(ncread('wrfout_d03_2000-12-01_13:00:00','V10',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));

WindSpeed_V = [V1; V2; V3; V4; V5; V6; V7; V8; V9; V10; V11; V12];

WindSpeed = sqrt(WindSpeed_U.^2 + WindSpeed_V.^2);

clearvars U1 U2 U3 U4 U5 U6 U7 U8 U9 U10 U11 U12 V1 V2 V3 V4 V5 V6 V7 ...
    V8 V9 V10 V11 V12 WindSpeed_U WindSpeed_V

%% Precipitation (mm/timestep)

% Accumulated total cumulus precipitation (kg m^-2 s^-1)
RAINC1 = squeeze(ncread('wrfout_d03_2000-01-01_00:00:00','RAINC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINC2 = squeeze(ncread('wrfout_d03_2000-01-31_13:00:00','RAINC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINC3 = squeeze(ncread('wrfout_d03_2000-03-02_01:00:00','RAINC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINC4 = squeeze(ncread('wrfout_d03_2000-04-01_13:00:00','RAINC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINC5 = squeeze(ncread('wrfout_d03_2000-05-02_01:00:00','RAINC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINC6 = squeeze(ncread('wrfout_d03_2000-06-01_13:00:00','RAINC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINC7 = squeeze(ncread('wrfout_d03_2000-07-02_01:00:00','RAINC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINC8 = squeeze(ncread('wrfout_d03_2000-08-01_13:00:00','RAINC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINC9 = squeeze(ncread('wrfout_d03_2000-09-01_01:00:00','RAINC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINC10 = squeeze(ncread('wrfout_d03_2000-10-01_13:00:00','RAINC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINC11 = squeeze(ncread('wrfout_d03_2000-11-01_01:00:00','RAINC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINC12 = squeeze(ncread('wrfout_d03_2000-12-01_13:00:00','RAINC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));

RAINC = [RAINC1; RAINC2; RAINC3; RAINC4; RAINC5; RAINC6; RAINC7; RAINC8; RAINC9; RAINC10; RAINC11; RAINC12];

% Accumulated total gridscale precipitation (kg m^-2 s^-1)
RAINNC1 = squeeze(ncread('wrfout_d03_2000-01-01_00:00:00','RAINNC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINNC2 = squeeze(ncread('wrfout_d03_2000-01-31_13:00:00','RAINNC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINNC3 = squeeze(ncread('wrfout_d03_2000-03-02_01:00:00','RAINNC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINNC4 = squeeze(ncread('wrfout_d03_2000-04-01_13:00:00','RAINNC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINNC5 = squeeze(ncread('wrfout_d03_2000-05-02_01:00:00','RAINNC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINNC6 = squeeze(ncread('wrfout_d03_2000-06-01_13:00:00','RAINNC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINNC7 = squeeze(ncread('wrfout_d03_2000-07-02_01:00:00','RAINNC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINNC8 = squeeze(ncread('wrfout_d03_2000-08-01_13:00:00','RAINNC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINNC9 = squeeze(ncread('wrfout_d03_2000-09-01_01:00:00','RAINNC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINNC10 = squeeze(ncread('wrfout_d03_2000-10-01_13:00:00','RAINNC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINNC11 = squeeze(ncread('wrfout_d03_2000-11-01_01:00:00','RAINNC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
RAINNC12 = squeeze(ncread('wrfout_d03_2000-12-01_13:00:00','RAINNC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));

RAINNC = [RAINNC1; RAINNC2; RAINNC3; RAINNC4; RAINNC5; RAINNC6; RAINNC7; RAINNC8; RAINNC9; RAINNC10; RAINNC11; RAINNC12];

% Total precipitation (kg m^-2 s^-1)
Precip = (RAINC + RAINNC);
Precip = [0; diff(Precip)];

% Convert precipitation to mm per timestep (mm/hr)
% Precip = Precip / 1000 * (60 * 60) / 1000;

clearvars RAINC1 RAINC2 RAINC3 RAINC4 RAINC5 RAINC6 RAINC7 RAINC8 ...
    RAINC9 RAINC10 RAINC11 RAINC12 RAINNC1 RAINNC2 RAINNC3 RAINNC4 ...
    RAINNC5 RAINNC6 RAINNC7 RAINNC8 RAINNC9 RAINNC10 RAINNC11 RAINNC12

%% Pressure (hPa)

PSFC1 = squeeze(ncread('wrfout_d03_2000-01-01_00:00:00','PSFC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
PSFC2 = squeeze(ncread('wrfout_d03_2000-01-31_13:00:00','PSFC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
PSFC3 = squeeze(ncread('wrfout_d03_2000-03-02_01:00:00','PSFC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
PSFC4 = squeeze(ncread('wrfout_d03_2000-04-01_13:00:00','PSFC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
PSFC5 = squeeze(ncread('wrfout_d03_2000-05-02_01:00:00','PSFC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
PSFC6 = squeeze(ncread('wrfout_d03_2000-06-01_13:00:00','PSFC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
PSFC7 = squeeze(ncread('wrfout_d03_2000-07-02_01:00:00','PSFC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
PSFC8 = squeeze(ncread('wrfout_d03_2000-08-01_13:00:00','PSFC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
PSFC9 = squeeze(ncread('wrfout_d03_2000-09-01_01:00:00','PSFC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
PSFC10 = squeeze(ncread('wrfout_d03_2000-10-01_13:00:00','PSFC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
PSFC11 = squeeze(ncread('wrfout_d03_2000-11-01_01:00:00','PSFC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
PSFC12 = squeeze(ncread('wrfout_d03_2000-12-01_13:00:00','PSFC',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));

Pressure = [PSFC1; PSFC2; PSFC3; PSFC4; PSFC5; PSFC6; PSFC7; PSFC8; PSFC9; PSFC10; PSFC11; PSFC12];

% Convert Pascals to hPa
Pressure = Pressure / 100;

clearvars PSFC1 PSFC2 PSFC3 PSFC4 PSFC5 PSFC6 PSFC7 PSFC8 PSFC9 PSCFC10 ...
    PSFC11 PSFC12

%% Specific Humidity (kg kg^-1)

Q2_1 = squeeze(ncread('wrfout_d03_2000-01-01_00:00:00','Q2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Q2_2 = squeeze(ncread('wrfout_d03_2000-01-31_13:00:00','Q2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Q2_3 = squeeze(ncread('wrfout_d03_2000-03-02_01:00:00','Q2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Q2_4 = squeeze(ncread('wrfout_d03_2000-04-01_13:00:00','Q2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Q2_5 = squeeze(ncread('wrfout_d03_2000-05-02_01:00:00','Q2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Q2_6 = squeeze(ncread('wrfout_d03_2000-06-01_13:00:00','Q2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Q2_7 = squeeze(ncread('wrfout_d03_2000-07-02_01:00:00','Q2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Q2_8 = squeeze(ncread('wrfout_d03_2000-08-01_13:00:00','Q2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Q2_9 = squeeze(ncread('wrfout_d03_2000-09-01_01:00:00','Q2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Q2_10 = squeeze(ncread('wrfout_d03_2000-10-01_13:00:00','Q2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Q2_11 = squeeze(ncread('wrfout_d03_2000-11-01_01:00:00','Q2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
Q2_12 = squeeze(ncread('wrfout_d03_2000-12-01_13:00:00','Q2',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));

SpecHumidity = [Q2_1; Q2_2; Q2_3; Q2_4; Q2_5; Q2_6; Q2_7; Q2_8; Q2_9; ...
    Q2_10; Q2_11; Q2_12];

clearvars Q2_1 Q2_2 Q2_3 Q2_4 Q2_5 Q2_6 Q2_7 Q2_8 Q2_9 Q2_10 Q2_11 Q2_12


%% Cloudiness (fraction)

CLDFRA1 = squeeze(ncread('wrfout_d03_2000-01-01_00:00:00','CLDFRA',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
CLDFRA2 = squeeze(ncread('wrfout_d03_2000-01-31_13:00:00','CLDFRA',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
CLDFRA3 = squeeze(ncread('wrfout_d03_2000-03-02_01:00:00','CLDFRA',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
CLDFRA4 = squeeze(ncread('wrfout_d03_2000-04-01_13:00:00','CLDFRA',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
CLDFRA5 = squeeze(ncread('wrfout_d03_2000-05-02_01:00:00','CLDFRA',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
CLDFRA6 = squeeze(ncread('wrfout_d03_2000-06-01_13:00:00','CLDFRA',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
CLDFRA7 = squeeze(ncread('wrfout_d03_2000-07-02_01:00:00','CLDFRA',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
CLDFRA8 = squeeze(ncread('wrfout_d03_2000-08-01_13:00:00','CLDFRA',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
CLDFRA9 = squeeze(ncread('wrfout_d03_2000-09-01_01:00:00','CLDFRA',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
CLDFRA10 = squeeze(ncread('wrfout_d03_2000-10-01_13:00:00','CLDFRA',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
CLDFRA11 = squeeze(ncread('wrfout_d03_2000-11-01_01:00:00','CLDFRA',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
CLDFRA12 = squeeze(ncread('wrfout_d03_2000-12-01_13:00:00','CLDFRA',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));

Cloudiness = [CLDFRA1; CLDFRA2; CLDFRA3; CLDFRA4; CLDFRA5; CLDFRA6; ...
    CLDFRA7; CLDFRA8; CLDFRA9; CLDFRA10; CLDFRA11; CLDFRA12];

clearvars CLDFRA1 CLDFRA2 CLDFRA3 CLDFRA4 CLDFRA5 CLDFRA6 CLDFRA7 ...
    CLDFRA8 CLDFRA9 CLDFRA10 CLDFRA11 CLDFRA12

%% Relative Humidity (%)

% Humidity = 100 * 0.263 * Pressure .* SpecHumidity .* (exp(17.67 * Temp ...
%     ./ (Temp + 273 - 29.65))) .^ -1;                                        % http://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
% Humidity(Humidity > 100) = 100;

% Saturation Vapor Pressure
SVP = 6.11 .* exp( 17.2694 .* (Temp) ./ (Temp + 238.3) );                   % http://www.conservationphysics.org/atmcalc/atmoclc2.pdf
% Saturation Mixing Ratio
SMR = 0.622 * (SVP ./ (Pressure - SVP));                                    % http://forum.wrfforum.com/viewtopic.php?f=12&t=9518

RelHumidity = 100 * SpecHumidity ./ SMR;
RelHumidity(RelHumidity > 100) = 100;

%% LW_in (W m^-2)

LW_1 = squeeze(ncread('wrfout_d03_2000-01-01_00:00:00','GLW',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
LW_2 = squeeze(ncread('wrfout_d03_2000-01-31_13:00:00','GLW',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
LW_3 = squeeze(ncread('wrfout_d03_2000-03-02_01:00:00','GLW',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
LW_4 = squeeze(ncread('wrfout_d03_2000-04-01_13:00:00','GLW',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
LW_5 = squeeze(ncread('wrfout_d03_2000-05-02_01:00:00','GLW',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
LW_6 = squeeze(ncread('wrfout_d03_2000-06-01_13:00:00','GLW',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
LW_7 = squeeze(ncread('wrfout_d03_2000-07-02_01:00:00','GLW',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
LW_8 = squeeze(ncread('wrfout_d03_2000-08-01_13:00:00','GLW',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
LW_9 = squeeze(ncread('wrfout_d03_2000-09-01_01:00:00','GLW',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
LW_10 = squeeze(ncread('wrfout_d03_2000-10-01_13:00:00','GLW',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
LW_11 = squeeze(ncread('wrfout_d03_2000-11-01_01:00:00','GLW',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
LW_12 = squeeze(ncread('wrfout_d03_2000-12-01_13:00:00','GLW',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));

LW_in = [LW_1; LW_2; LW_3; LW_4; LW_5; LW_6; LW_7; LW_8; LW_9; ...
    LW_10; LW_11; LW_12];

clearvars LW_1 LW_2 LW_3 LW_4 LW_5 LW_6 LW_7 LW_8 LW_9 LW_10 LW_11 LW_12

%% Emissivity (W m^-2)

% Emiss_1 = squeeze(ncread('wrfout_d03_2000-01-01_00:00:00','EMISS',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
% Emiss_2 = squeeze(ncread('wrfout_d03_2000-01-31_13:00:00','EMISS',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
% Emiss_3 = squeeze(ncread('wrfout_d03_2000-03-02_01:00:00','EMISS',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
% Emiss_4 = squeeze(ncread('wrfout_d03_2000-04-01_13:00:00','EMISS',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
% Emiss_5 = squeeze(ncread('wrfout_d03_2000-05-02_01:00:00','EMISS',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
% Emiss_6 = squeeze(ncread('wrfout_d03_2000-06-01_13:00:00','EMISS',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
% Emiss_7 = squeeze(ncread('wrfout_d03_2000-07-02_01:00:00','EMISS',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
% Emiss_8 = squeeze(ncread('wrfout_d03_2000-08-01_13:00:00','EMISS',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
% Emiss_9 = squeeze(ncread('wrfout_d03_2000-09-01_01:00:00','EMISS',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
% Emiss_10 = squeeze(ncread('wrfout_d03_2000-10-01_13:00:00','EMISS',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
% Emiss_11 = squeeze(ncread('wrfout_d03_2000-11-01_01:00:00','EMISS',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
% Emiss_12 = squeeze(ncread('wrfout_d03_2000-12-01_13:00:00','EMISS',[m_Idx n_Idx 1],[1 1 732],[1 1 1]));
% 
% Emissivity = [Emiss_1; Emiss_2; Emiss_3; Emiss_4; Emiss_5; Emiss_6; ...
%     Emiss_7; Emiss_8; Emiss_9; Emiss_10; Emiss_11; Emiss_12];
% 
% clearvars Emiss_1 Emiss_2 Emiss_3 Emiss_4 Emiss_5 Emiss_6 Emiss_7 ...
%     Emiss_8 Emiss_9 Emiss_10 Emiss_11 Emiss_12

%% Other Data

% Temporal resolution of the raw data (minutes)
Data.Resolution = 60; % 1 hour
% Data.AWS_Lat = Lat;
% Data.AWS_Long = Long;
Data.TimeZone = 0;
Data.FirstDay = 0;









