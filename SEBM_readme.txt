Input data files are:
DS_wRGI.csv**
UIB_under15pct_debris.mat
bisquare_fits_threshold25.mat
PC_20km10.mat
calibration_march.mat for calibration process (manually remove glacier 232 as it is not modeled in full UIB population)

ExtractHAR will need to be edited to reflect the location and format of your downloaded HAR data.  ExtractHAR.m is provided as an example and to show where it enters the broader code.

Main code is Run_Parfor_Glaciers.m -> Full_Model.m*
variations are:
FOR CALIBRATING THE PCF: Run_Parfor_Calibration.m -> calibrations.m -> full_model.m
FOR 12-GLACIER UNCERTAINTY ANALYSIS: Uncertainty_Parfor_Calibration.m -> Uncertainty_calibration.m -> Uncertainty_model.m

In all top-level (parfor) scrips, redesignate directory locations for the following:
GUI_Input.output_filename   = ('/uufs/chpc.utah.edu/common/home/u6027899/SEBM_output/Glacier_Number_');
GUI_Input.cal_filename      = ('/uufs/chpc.utah.edu/common/home/u6027899/SEBM_output/');
GUI_Input.GlacNum_filename  = ('/uufs/chpc.utah.edu/common/home/u6027899/ASTER/ASTER_90m_Numbers.mat');

Full_Model.m (and variations) notes:
- variables have prefixes: m (matrix), k (constant), v (vector), st (structure), etc.
- contains a structure for time; use can apply the date time built-in matlab function with datetime(stTime.year(t),stTime.month(t),stTime.day(t))
- NOTE: in time variables, the first hour of a day (hour 0) is considered part of the previous day (hour 24).  Days change at 1 am instead of midnight.
- the time vector is 131496 long, running Jan 1, 2000 0:00 -- Jan 1, 2015 0:00
but the HAR forcing is only t = 8785 (Jan 1, 2001 01:00) through t = 122713 (Dec 31, 2013 1:00) +23h added in code by A. Giese

Code necessary for generating Figure 1 in Giese & al was written by D. Keeler (durbank777@gmail.com) and is available at  <https://drive.google.com/drive/folders/1S3xxoWpfewmEHPkc9ZqJo6ACrK75gkXC>.  The basemap comes from an ESRI layer of USGS topo maps. The Python functions pull from this server: https://server.arcgisonline.com/arcgis/rest/services/USA_Topo_Maps/MapServer



** Columns are: (field names from Shean & al., 2020 and RGI)
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