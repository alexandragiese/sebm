clear
% addpath('/uufs/chpc.utah.edu/common/home/u6027899/sebm/SkyViewFactor');
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
GUI_Input.precip_threshold  = 250e-6; %(m/d) OLD: 1e-2; %(m day^-1) 
GUI_Input.snow_density      = 330;
GUI_Input.lapse_rate        = 6.500;
GUI_Input.iZ_0m_snow        = 0.001;
GUI_Input.iZ_0m_ice         = 0.016;
GUI_Input.iZ_0T_snow        = 0.001;
GUI_Input.iZ_0T_ice         = 0.004;
GUI_Input.iZ_0q_snow        = 0.001;
GUI_Input.iZ_0q_ice         = 0.004;
GUI_Input.threeDsave        = 0;

% Import Geodedic Mass Balance Data----------------------------------------
GUI_Input.geo = csvread('DS_wRGI.csv');
% glID  = GUI_Input.geo(:,1);
% geoMB = GUI_Input.geo(:,11); % m w.e. / a
% geoMB_sigma = GUI_Input.geo(:,12);
% basin = GUI_Input.geo(:,49);

% Glacier IDs:-------------------------------------------------------------
load gapfilled_subset.mat
    foo = ~isnan(glac_nums(:));
    y = glac_nums(foo);
    
% Precipitation lapse rate:------------------------------------------------
load bisquare_fits_threshold25.mat %.00025/day variable: PLR
GUI_Input.PLR = PLR;

% Precipitation Correction (Calibration) factor:---------------------------
% cal_PC = zeros(length(y),1); 
load([GUI_Input.cal_filename, 'calibration_march.mat'],'cal_PC')
GUI_Input.PC = cal_PC;

% MB_corr = nan(length(y),1);
% save([GUI_Input.cal_filename, 'calibration_march.mat'],'cal_PC','MB_corr') 

poolobj = parpool('local',4);
parfor g = 151:324 %DONT OVERWRITE CALIBRATION .MAT FILE WHEN 59:length(y) % g is the index in the LIST of glacier numbers; y(g) is the glacier number %CAL
%     disp(g)
    g
    calibration(GUI_Input,y,g)
end
toc
delete(gcp('nocreate'))