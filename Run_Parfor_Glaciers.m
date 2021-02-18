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

% Import Geodedic Mass Balance Data----------------------------------------
GUI_Input.geo = csvread('DS_wRGI.csv');
% glID  = GUI_Input.geo(:,1);
% geoMB = GUI_Input.geo(:,11); % m w.e. / a
% geoMB_sigma = GUI_Input.geo(:,12);
% basin = GUI_Input.geo(:,49);

% Glacier IDs:-------------------------------------------------------------
% % load deb10size10th5HAR.mat 
% % load Indus1_rndm30.mat
% % load 30gl_eachbasin_90mRES.mat
% % load clean_median.mat
load gapfilled_subset.mat
% load UIB_under15pct_debris.mat; glac_nums = subbasin_sub;
    foo = ~isnan(glac_nums(:));
    y = glac_nums(foo);
    
% Precipitation lapse rate:------------------------------------------------
% % load linear_fits.mat %contains 11 x 3 PLR
% % load bisquare_fits.mat
% % load bisquare_fits_10basins.mat
% % load bisquare_fits_threshold.mat %.000125/day
load bisquare_fits_threshold25.mat %.00025/day variable: PLR
GUI_Input.PLR = PLR;

% Precipitation Correction (Calibration) factor:---------------------------
% % load precip_corr_10b.mat
% % load PC_by_glacier_th25.mat
% %     foo = ~isnan(C(:));
% %     GUI_Input.PC = C(foo);
% load subset_PC_latlon.mat %calibration subset: glacier number, PC, lat, lon
load PC_iterated.mat
GUI_Input.PC = subset(:,2); %<-- for checking calibration! 
% % cal_PC = zeros(length(y),1); % CAL
% % GUI_Input.PC = cal_PC;
% load PC_fullUIB.mat          %<-- updated one basin at a time with different radius thresholds (perbasin.m local)
%     GUI_Input.PC = PC_calc;

% MB_corr = nan(length(y),1); % CAL
% redo = [46    55    56    62    64    69    75    80   101   106   111   136   137   159   196   200   210   212   225   229   233   237   238   244   253   285   292   297   300   303   307   309   312   318   320   321   322   324];
% ct = 1; %CAL counter
mdl = zeros(1,length(y));

poolobj = parpool('local',4);
parfor g = 1:length(y) % g is the index in the LIST of glacier numbers; y(g) is the glacier number %CAL
    disp(g)
%     load([GUI_Input.cal_filename, 'calibration_fix.mat'],'cal_PC','MB_corr')
    [R1, ~, ~ ] = Full_Model_AG_edits( GUI_Input,y,g);
% while R1.geodeticMB < 0
%     if R1.modelMB > 0.9* R1.geodeticMB
%         ct = ct+1;
%         GUI_Input.PC(g) = GUI_Input.PC(g) - abs( (R1.geodeticMB-R1.modelMB)/ct );
%         [R1, ~, ~ ] = Full_Model_AG_edits( GUI_Input,y,g); 
%         mdl = R1.modelMB
%         ct
%         corr_fac = GUI_Input.PC(g)  
%     elseif R1.modelMB < 1.1* R1.geodeticMB
%         ct = ct+1;
%         GUI_Input.PC(g) = GUI_Input.PC(g) + abs( (R1.geodeticMB-R1.modelMB)/ct );
%         [R1, ~, ~ ] = Full_Model_AG_edits( GUI_Input,y,g); R1.modelMB
%         mdl = R1.modelMB
%         ct
%         corr_fac = GUI_Input.PC(g)  
%     elseif R1.modelMB >= 1.1* R1.geodeticMB && R1.modelMB <= 0.9* R1.geodeticMB
%         break
%     end
%  end
% while R1.geodeticMB > 0
%     if R1.modelMB < 0.9* R1.geodeticMB
%         ct = ct+1
%         GUI_Input.PC(g) = GUI_Input.PC(g) + abs( (R1.geodeticMB-R1.modelMB)/ct );
%         [R1, ~, ~ ] = Full_Model_AG_edits( GUI_Input,y,g);
%         mdl = R1.modelMB
%         corr_fac = GUI_Input.PC(g)  
%     elseif R1.modelMB > 1.1* R1.geodeticMB
%         ct = ct+1
%         GUI_Input.PC(g) = GUI_Input.PC(g) - abs( (R1.geodeticMB-R1.modelMB)/ct );
%         [R1, ~, ~ ] = Full_Model_AG_edits( GUI_Input,y,g);
%         mdl = R1.modelMB
%         corr_fac = GUI_Input.PC(g)         
%     elseif R1.modelMB <= 1.1* R1.geodeticMB && R1.modelMB >= 0.9* R1.geodeticMB
%         break
%     end   
% end
% while R1.geodeticMB == 0
%     if R1.modelMB < -0.0194191
%         ct = ct+1
%         GUI_Input.PC(g) = GUI_Input.PC(g) + abs( (R1.geodeticMB-R1.modelMB)/ct );
%         [R1, ~, ~ ] = Full_Model_AG_edits( GUI_Input,y,g);
%         mdl = R1.modelMB
%         corr_fac = GUI_Input.PC(g)         
%     elseif R1.modelMB > 0.0194191
%         ct = ct+1
%         GUI_Input.PC(g) = GUI_Input.PC(g) - abs( (R1.geodeticMB-R1.modelMB)/ct );
%         [R1, ~, ~ ] = Full_Model_AG_edits( GUI_Input,y,g);
%         mdl = R1.modelMB
%         corr_fac = GUI_Input.PC(g)         
%     elseif R1.modelMB >= -0.0194191 && R1.modelMB <= 0.0194191
%         break
%     end
%  end

mdl(g) = R1.modelMB
% corr_fac = GUI_Input.PC(g)
%     cal_PC(g) = GUI_Input.PC(g); % <-- SAVE THISb  %CAL
%     MB_corr(g) = R1.modelMB;
%     save([GUI_Input.cal_filename, 'calibration_fix.mat'],'cal_PC','MB_corr') 
end
save([GUI_Input.cal_filename, 'calibration_check.mat'],'mdl') 
toc
delete(gcp('nocreate'))