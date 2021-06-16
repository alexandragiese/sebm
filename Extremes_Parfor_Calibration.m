clear
addpath(genpath('/uufs/chpc.utah.edu/common/home/u6027899/sebm'));

GUI_Input.sGCM          = 'HAR';
GUI_Input.HAR_temp_res  = 'hourly'; 
GUI_Input.start_date    = '01/01/2001';
GUI_Input.end_date      = '12/31/2013';
GUI_Input.kRadius       = 98;
GUI_Input.output_filename   = ('/uufs/chpc.utah.edu/common/home/cryosphere/agiese/uncertainty/Glacier_Number_');
GUI_Input.cal_filename      = ('/uufs/chpc.utah.edu/common/home/u6027899/SEBM_output/');
GUI_Input.GlacNum_filename  = ('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_90m_Numbers.mat');
GUI_Input.ensemble_number   = 8;
GUI_Input.precip_threshold  = 250e-6; %(m/d) OLD: 1e-2; %(m day^-1) 
GUI_Input.snow_density      = 330;
GUI_Input.threeDsave        = 0;

%     cal_PC    = nan(12,1);
%     MB_corr   = nan(12,1); 
%     MELT_corr = nan(12,1);
%     geoMB     = nan(12,1);
%     varied_PC   = nan(12,1);
%     varied_MB   = nan(12,1);
%     varied_MELT = nan(12,1);

% Import Geodedic Mass Balance Data----------------------------------------
GUI_Input.geo = csvread('DS_wRGI.csv');

% Glacier IDs:-------------------------------------------------------------
load UIB_under15pct_debris.mat; glac_nums = subbasin_sub;
    foo = ~isnan(glac_nums(:));
    y = glac_nums(foo);
%     g = 736; %glacier number
        
% Precipitation lapse rate:------------------------------------------------
load bisquare_fits_threshold25.mat %.00025/day variable: PLR
GUI_Input.PLR = PLR;

sp12 = [736 6355 4036 11031 8629 12850 5107 2893 6346 12188 10593 2880]; %the special 12 for uncertainty analysis (via extremes)

% Create results file-------------------------
% x = length(sp12);
%     cal_PC      = nan(x,1);
%     MB_corr     = nan(x,1);
%     MELT_corr   = nan(x,1);
%     geoMB       = nan(x,1);
%     varied_PC   = nan(x,1);
%     varied_MB   = nan(x,1);
%     varied_MELT = nan(x,1);
%     iGlacierNumber = y(g);
% save([GUI_Input.cal_filename,'Extremes.mat'],'cal_PC','MB_corr','MELT_corr','geoMB','varied_PC','varied_MB','varied_MELT')

j = [3, 5, 7, 12];
poolobj = parpool('local',4);
parfor k = 1:4 
    g = sp12(j(k));
    y(g)
    Extremes_calibration(GUI_Input,y,g,k)
end

delete(gcp('nocreate'))