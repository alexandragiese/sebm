clear
tic

GUI_Input.sGCM          = 'HAR';
GUI_Input.HAR_temp_res  = 'hourly'; % 'daily' or 'hourly'
GUI_Input.start_date    = '01/01/2001';
GUI_Input.end_date      = '12/31/2013';
GUI_Input.kRadius       = 98;
GUI_Input.output_filename   = ('/uufs/chpc.utah.edu/common/home/u6027899/SEBM_output/Glacier_Number_');
GUI_Input.cal_filename      = ('/uufs/chpc.utah.edu/common/home/u6027899/SEBM_output/');
GUI_Input.GlacNum_filename  = ('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_90m_Numbers.mat');
GUI_Input.ensemble_number   = 8;
GUI_Input.precip_threshold  = 250e-6; %(m/d)  
GUI_Input.snow_density      = 330; %Arnold
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
load UIB_under15pct_debris.mat; glac_nums = subbasin_sub;
    foo = ~isnan(glac_nums(:));
    y = glac_nums(foo);
    
% Precipitation lapse rate:------------------------------------------------
load bisquare_fits_threshold25.mat %.00025/day. variable: PLR
GUI_Input.PLR = PLR;

% Precipitation Correction (Calibration) factor:---------------------------
load PC_20km10.mat
    GUI_Input.PC = PC_calc;
% GUI_Input.PC = zeros(13971,1); %for no PCF applied

poolobj = parpool('local',4);
parfor g = 1:4960 % g is the index in the LIST of glacier numbers; y(g) is the glacier number. ADJUST THIS AS NECESSARY
    disp(g)
    [R1, ~, ~ ] = Full_Model( GUI_Input,y,g);
%     mdl(g)  = R1.modelMB;
%     geo(g)  = R1.geodeticMB;
%     melt(g) = R2.TotalAveMelt/13; 
%     accum(g)= R2.TotalAveAccum/13; 
%     sub(g)  = R2.TotalAveSub/13; 
end
toc
delete(gcp('nocreate'))