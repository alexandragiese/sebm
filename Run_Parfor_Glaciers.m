clear
% addpath('/uufs/chpc.utah.edu/common/home/u6027899/sebm/SkyViewFactor');
tic

GUI_Input.sGCM          = 'HAR';
GUI_Input.HAR_temp_res  = 'hourly'; % 'daily' or 'hourly'
GUI_Input.start_date    = '01/01/2001';
GUI_Input.end_date      = '12/31/2013';
GUI_Input.kRadius       = 98;
% GUI_Input.output_filename   = ['/uufs/chpc.utah.edu/common/home/steenburgh-group6/Ali2/SEBM_output/Glacier_Number_'];
GUI_Input.output_filename   = ('/uufs/chpc.utah.edu/common/home/u6027899/SEBM_output/Glacier_Number_');
GUI_Input.cal_filename      = ('/uufs/chpc.utah.edu/common/home/u6027899/SEBM_output/');
GUI_Input.GlacNum_filename  = ('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_90m_Numbers.mat');
GUI_Input.ensemble_number   = 8;
GUI_Input.precip_threshold  = 250e-6; %(m/d) OLD: 1e-2; %(m day^-1) 
GUI_Input.snow_density      = 330; %Arnold. Eric: 250;
% GUI_Input.lapse_rate        = 6.500;
GUI_Input.iZ_0m_snow        = 0.001;
GUI_Input.iZ_0m_ice         = 0.016;
GUI_Input.iZ_0T_snow        = 0.001;
GUI_Input.iZ_0T_ice         = 0.004;
GUI_Input.iZ_0q_snow        = 0.001;
GUI_Input.iZ_0q_ice         = 0.004;
GUI_Input.threeDsave        = 0;

% Import Geodedic Mass Balance Data----------------------------------------
GUI_Input.geo = csvread('DS_wRGI.csv');

% Glacier IDs:-------------------------------------------------------------
% load gapfilled_subset.mat
load UIB_under15pct_debris.mat; glac_nums = subbasin_sub;
    foo = ~isnan(glac_nums(:));
    y = glac_nums(foo);
    
% Precipitation lapse rate:------------------------------------------------
load bisquare_fits_threshold25.mat %.00025/day variable: PLR
GUI_Input.PLR = PLR;

% Precipitation Correction (Calibration) factor:---------------------------
% cal_PC = zeros(length(y),1); % CAL
% load([GUI_Input.cal_filename, 'calibration_march.mat'],'cal_PC')
%   GUI_Input.PC = cal_PC;
% load PC_fullUIB.mat          %<-- updated one basin at a time with different radius thresholds (perbasin.m local)
load PC_20km10.mat
    GUI_Input.PC = PC_calc;
% GUI_Input.PC = zeros(13971,1);

% poolobj = parpool('local',4);
% par
for g = 2302:4960 % g is the index in the LIST of glacier numbers; y(g) is the glacier number %CAL
    disp(g)
    [R1, ~, ~ ] = Full_Model( GUI_Input,y,g);
%     mdl(g)  = R1.modelMB;
%     geo(g)  = R1.geodeticMB;
%     area(g) = R1.area;
%     melt(g) = R2.TotalAveMelt/13; 
%     accum(g)= R2.TotalAveAccum/13; 
%     sub(g)  = R2.TotalAveSub/13; 
end
toc
delete(gcp('nocreate'))