clear

tic

GUI_Input.sGCM          = 'HAR';
GUI_Input.HAR_temp_res  = 'daily'; % 'daily' or 'hourly'
GUI_Input.start_date    = '01/01/2001';
GUI_Input.end_date      = '12/31/2013';
% GUI_Input.start_date    = '01/01/2001';
% GUI_Input.end_date      = '12/31/2001';
GUI_Input.kRadius       = 98;
% GUI_Input.output_filename   = ['/uufs/chpc.utah.edu/common/home/steenburgh-group6/Ali2/SEBM_output/Glacier_Number_'];
GUI_Input.output_filename   = ['/uufs/chpc.utah.edu/common/home/u6027899/SEBM_output/Glacier_Number_'];
GUI_Input.ensemble_number   = 8;
GUI_Input.mb_date_1975      = '01/01/1975';
GUI_Input.mb_date_2000      = '12/31/1999';
GUI_Input.precip_threshold  = 1e-2; % (m day^-1)
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
G     = csvread('DS_wRGI.csv');
glID  = G(:,1);
geoMB = G(:,11); % m w.e. / a
geoMB_sigma = G(:,12);
basin = G(:,49);

% Glacier IDs to run
% % load deb10size10th5HAR.mat 
% % load Indus1_rndm30.mat
load 30gl_eachbasin.mat
% load clean_median.mat
    foo = ~isnan(glac_nums(:));
    y = glac_nums(foo);
% load linear_fits.mat %contains 11 x 3 PLR
load bisquare_fits.mat
% load precip_corr.mat
%     foo = ~isnan(C(:));
%     GUI_Input.PC = C(foo);
    
% S: solution vector: 
S = nan(length(y),2);
i = -1;

for g = 121:150 %1:length(y) %r = 5006 
    GUI_Input.glacier_number = y(g); %= glID(r); 
    r = find (glID==y(g));

    GUI_Input.glacier_basin = basin(r);
%     GUI_Input.p1 = PLR(basin(r),1);
%     GUI_Input.p2 = PLR(basin(r),2);
%     GUI_Input.r2 = PLR(basin(r),3);
    
    [R1, R2, RAve] = Full_Model_AG3( GUI_Input,g );
%     modelMB    = (R2.TotalAveMelt+R2.TotalAveAccum+R2.TotalAveSub)/13;
%     geodeticMB = geoMB(r);
%     sigma      = geoMB_sigma(r);
        S(g*2+i:g*2+i+2,1) = g;
        S(g*2+i,2)   = (R2.TotalAveMelt+R2.TotalAveAccum+R2.TotalAveSub)/13; 
        S(g*2+i+1,2) = geoMB(r);
        S(g*2+i+2,2) = geoMB_sigma(r);
    save([GUI_Input.output_filename, num2str(GUI_Input.glacier_number),'_TEST.mat'])
    keyboard
    i = i+1;
    disp(g)
end


% IF BIG: save inside full_model script

toc




%Attributes csv---------------
% 1)  RGIId
% 2)  x          
% 3)  y          
% 4)  z_med     
% 5)  z_min      
% 6)  z_max      
% 7)  z_slope     
% 8)  z_aspect    
% 9)  dhdt_ma      
% 10) dhdt_ma_sigma  
% 11) mb_mwea     
% 12) mb_mwea_sigma  
% 13) area_m2        
% 14) mb_m3wea      
% 15) mb_m3wea_sigma  
% 16) t1           
% 17) t2        
% 18) dt        
% 19) valid_area_perc 
% 20) H_m        
% 21) debris_m    
% 22) perc_debris 
% 23) perc_pond  
% 24) perc_clean  
% 25) vm_ma    
% 26) BgnDate   
% 27) EndDate  
% 28) CenLon   
% 29) CenLat   
% 30) O1Region 
% 31) O2Region 
% 32) Area     
% 33) Zmin     
% 34) Zmax      
% 35) Zmed    
% 36) Slope    
% 37) Aspect 
% 38) Lmax    
% 39) Status  
% 40) Connect     
% 41) Form        
% 42) TermType    
% 43) Surging     
% 44) Linkages    
% 45) DC_Area     
% 46) DC_BgnDate	
% 47) DC_EndDate	
% 48) DC_CTSmean  
% 49) Subbasin 1-11

% plot(RAve.MB_track(2:end)-RAve.MB_track(1:end-1),'.','MarkerSize',20)
% ylim([-.27 -.265])
% ylabel('Change in MB (m w.e.)')
% xlabel('Year after initial (2001)')