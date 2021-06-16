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
% GUI_Input.lapse_rate        = 6.500;
GUI_Input.threeDsave        = 0;

% Import Geodedic Mass Balance Data----------------------------------------
GUI_Input.geo = csvread('DS_wRGI.csv');

% Glacier IDs:-------------------------------------------------------------
load UIB_under15pct_debris.mat; glac_nums = subbasin_sub;
    foo = ~isnan(glac_nums(:));
    y = glac_nums(foo);
    g = 736; %glacier number
        % Coordinates calculated in Full_Model.  Pasted here for each of
        % the 12 glaciers so that vP_d with the threshold is calculated
        % only once
%         kLat = 
%         kLong = 
        
% Precipitation lapse rate:------------------------------------------------
load bisquare_fits_threshold25.mat %.00025/day variable: PLR
GUI_Input.PLR = PLR;

% Create results file-------------------------
x = 1000;
    cal_PC      = nan(x,1);
    MB_corr     = nan(x,1);
    MELT_corr   = nan(x,1);
    geoMB       = nan(x,1);
    varied_PC   = nan(x,1);
    varied_MB   = nan(x,1);
    varied_MELT = nan(x,1);
    iGlacierNumber = y(g);
save([GUI_Input.cal_filename,'Glacier_Number_',num2str(iGlacierNumber),'MonteCarlo.mat'],'cal_PC','MB_corr','MELT_corr','geoMB','varied_PC','varied_MB','varied_MELT')

%Slow application of precipitation threshold (do once when parfor loop is multiple runs for the same glacier!)
sDataDirectory = '/uufs/chpc.utah.edu/common/home/steenburgh-group6/Eric/HAR_Data'; %for all HAR
% if strcmp(GUI_Input.sGCM,'HAR')
%     sTempRes    = GUI_Input.HAR_temp_res;
%     [ ~, ~, ~, vP_d, ~, ~, ~, ~, ~, ~  ] = ExtractHAR1( kLat, kLong, sTempRes, sDataDirectory ); %pulls in all HAR output
% end
vP_d_th = nan(size(vP_d)); 
for r = 1:length(vP_d) %42 seconds (3.6 s before this)
    precip_series = reshape(vP_d,24,[]);
    daily_precip  = sum(precip_series); %m/d
    A = find(daily_precip <= GUI_Input.precip_threshold); % Precipitation event threshold (m day^-1)
    precip_series(:,A) = 0;
    vP_d_th = precip_series(:);    
end
% checked saved file for R_k
% reinstate PARFOR, change to random selection of parameters, remove keyboard & tic/toc
% poolobj = parpool('local',4);
% par
for mc = 1:50 %x 
%     mc
    MCcalibration(GUI_Input,y,g,mc,vP_d_th)
end

delete(gcp('nocreate'))