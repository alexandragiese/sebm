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
GUI_Input.GlacNum_filename  = ('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_90m_Numbers.mat');
GUI_Input.ensemble_number   = 8;
GUI_Input.precip_threshold  = 250e-6; %(m/d) OLD: 1e-2; %(m day^-1) 
GUI_Input.snow_depth        = 0;
GUI_Input.snow_density      = 250;
% GUI_Input.temp_forcing      = 0;
GUI_Input.lapse_rate        = 6.500;
GUI_Input.iZ_0m_snow        = 0.001;
GUI_Input.iZ_0m_ice         = 0.016;
GUI_Input.iZ_0T_snow        = 0.001;
GUI_Input.iZ_0T_ice         = 0.004;
GUI_Input.iZ_0q_snow        = 0.001;
GUI_Input.iZ_0q_ice         = 0.004;
GUI_Input.threeDsave        = 0;

% Import Geodedic Mass Balance Data
GUI_Input.geo = csvread('DS_wRGI.csv');
% glID  = GUI_Input.geo(:,1);
% geoMB = GUI_Input.geo(:,11); % m w.e. / a
% geoMB_sigma = GUI_Input.geo(:,12);
% basin = GUI_Input.geo(:,49);

% Glacier IDs:
% load deb10size10th5HAR.mat 
% load Indus1_rndm30.mat
% load 30gl_eachbasin_90mRES.mat
% load clean_median.mat
% load gapfilled_subset.mat
load UIB_under15pct_debris.mat; glac_nums = subbasin_sub;
    foo = ~isnan(glac_nums(:));
    y = glac_nums(foo);
    
% Precipitation lapse rate:
% load linear_fits.mat %contains 11 x 3 PLR
% load bisquare_fits.mat
% load bisquare_fits_10basins.mat
% load bisquare_fits_threshold.mat %.000125/day
load bisquare_fits_threshold25.mat %.00025/day variable: PLR
GUI_Input.PLR = PLR;

% load precip_corr_10b.mat
% load PC_by_glacier_th25.mat
%     foo = ~isnan(C(:));
%     GUI_Input.PC = C(foo);
load PC_fullUIB.mat
    GUI_Input.PC = PC_calc;
% load subset_PC_latlon.mat
% GUI_Input.PC = subset(:,2); %<-- for checking calibration!
%     

% % S: solution vector: 
% S = nan(length(y)*3,2);
% i = -1;

poolobj = parpool('local',4);
parfor g = 10505:11859
    disp(g)
    Full_Model_AG_edits( GUI_Input,y,g);
end

toc
delete(gcp('nocreate'))