tic

GUI_Input.sGCM          = 'HAR';
GUI_Input.HAR_temp_res  = 'hourly'; % 'daily' or 'hourly'
GUI_Input.start_date    = '01/01/2001';
GUI_Input.end_date      = '12/31/2013';
GUI_Input.kRadius       = 98;
GUI_Input.output_filename   = ['/uufs/chpc.utah.edu/common/home/u6027899/SEBM_output/Glacier_Number_'];
GUI_Input.ensemble_number   = 8;
GUI_Input.mb_date_1975      = '01/01/1975';
GUI_Input.mb_date_2000      = '12/31/1999';
GUI_Input.precip_threshold  = 1e-2; % (m day^-1)
GUI_Input.snow_depth        = 0;
GUI_Input.snow_density      = 250;
GUI_Input.temp_forcing      = 0;
GUI_Input.lapse_rate        = 6.500;
GUI_Input.precip_maxa_alt   = 5500; % m a.s.l. AG
GUI_Input.iZ_0m_snow        = 0.001;
GUI_Input.iZ_0m_ice         = 0.016;
GUI_Input.iZ_0T_snow        = 0.001;
GUI_Input.iZ_0T_ice         = 0.004;
GUI_Input.iZ_0q_snow        = 0.001;
GUI_Input.iZ_0q_ice         = 0.004;
GUI_Input.threeDsave        = 1;
GUI_Input.cf                = 1; %correction factor for PLR

% Import Geodedic Mass Balance Data
% load deb10size10th5HAR.mat %IDs to run
load Indus1_rndm30.mat
G     = csvread('DS_wRGI.csv');
glID  = G(:,1);
geoMB = G(:,11); % m w.e. / a
geoMB_sigma = G(:,12);

%these while statements may not work when HAR elev is above glacier?
for g = 3 %  r = 5006 %
    GUI_Input.glacier_number    = y(g); %glID(r); %
    r = find (glID==y(g));
    [R1, R2, RAve] = Full_Model_AG2( GUI_Input );
        modelMB    = (R2.TotalAveMelt+R2.TotalAveAccum+R2.TotalAveSub)/13
        geodeticMB = geoMB(r)
        sigma      = geoMB_sigma(r)
%     while modelMB < geodeticMB - sigma %model melts too much --> precipitate more (increase lapse rate)
%         f = abs(abs(geodeticMB-modelMB)/geodeticMB); %factor, 1 would be exact match
%             if f >= 10
%                 disp(['Resetting f=',num2str(f),' to 9.5.'])
%                 f = 9.5;
%             end
%         GUI_Input.cf   = GUI_Input.cf*(1+f/10); %adjust by how much off; 0.1 may need to be adjusted
%         [R1, R2, RAve] = Full_Model_AG2( GUI_Input );
%         modelMB    = (R2.TotalAveMelt+R2.TotalAveAccum+R2.TotalAveSub)/13
%     end
%     while modelMB > geodeticMB + sigma %model melts too little --> precipitate less (decrease lapse rate)
%         f = abs(abs(geodeticMB-modelMB)/geodeticMB);
%             if f >= 10
%                 disp(['Resetting f=',num2str(f),' to 9.5.'])
%                 f = 9.5;
%             end
%         GUI_Input.cf   = GUI_Input.cf*(1-f/10); %adjust by how much off
%         [R1, R2, RAve] = Full_Model_AG2( GUI_Input );
%         modelMB    = (R2.TotalAveMelt+R2.TotalAveAccum+R2.TotalAveSub)/13
%     end
%     if modelMB < geodeticMB+sigma && modelMB > geodeticMB-sigma %this if statement *shouldn't* be necessary
%         save([GUI_Input.output_filename, num2str(GUI_Input.glacier_number),'3D_PLRnoI.mat'])
%     end
%     g
end
varinfo=whos('m3TotalMelt');
saveopt='';
if varinfo.bytes >= 2^31
  saveopt='-v7.3';
end

% IF BIG: save([GUI_Input.output_filename, num2str(iGlacierNumber),'_AG.mat'],saveopt)
save([GUI_Input.output_filename, num2str(GUI_Input.glacier_number),'3D_PLRnoI.mat'],saveopt)

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
% 44) TermType    
% 45) Surging     
% 46) Linkages    
% 47) DC_Area     
% 48) DC_BgnDate	
% 49) DC_EndDate	
% 50) DC_CTSmean  

  