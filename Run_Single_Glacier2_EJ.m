clock

tic

GUI_Input.sGCM = 'HAR';
GUI_Input.HAR_temp_res = 'hourly'; % 'daily' or 'hourly'
GUI_Input.start_date = '01/01/2001';
GUI_Input.end_date = '12/31/2013';
GUI_Input.output_filename   = ('/uufs/chpc.utah.edu/common/home/u6027899/SEBM_output/Glacier_Number_');
GUI_Input.ensemble_number = 8;
GUI_Input.mb_date_1975 = '01/01/1975';
GUI_Input.mb_date_2000 = '12/31/1999';
GUI_Input.precip_threshold = 1e-2; % (m day^-1)
GUI_Input.snow_depth = 0;
GUI_Input.snow_density = 250;
GUI_Input.temp_forcing = 0;
GUI_Input.lapse_rate = 6.500;
GUI_Input.iZ_0m_snow = 0.001;
GUI_Input.iZ_0m_ice = 0.016;
GUI_Input.iZ_0T_snow = 0.001;
GUI_Input.iZ_0T_ice = 0.004;
GUI_Input.iZ_0q_snow = 0.001;
GUI_Input.iZ_0q_ice = 0.004;

addpath '/uufs/chpc.utah.edu/common/home/steenburgh-group6/Ali/SEBM_from_Eric'
vJosh_MB_Data = load('Josh_Geo_MB_Data1.mat');
vJosh_Object_ID = vJosh_MB_Data.vJosh_Object_ID;
vcJosh_Category = vJosh_MB_Data.vcJosh_Category;
vJosh_PctDeb = vJosh_MB_Data.vJosh_PctDeb;
vJosh_GeoMassBal = vJosh_MB_Data.vJosh_GeoMassBal;

vModel_MB_1975_to_2000 = NaN(length(vJosh_Object_ID),1);
vModel_MB_2000_to_2016 = NaN(length(vJosh_Object_ID),1);
vModel_MB = NaN(length(vJosh_Object_ID),1);


GUI_Input.glacier_number = 56076;

[R1, R2, RAve] = Full_Model_Nov15_2019email( GUI_Input );

save([GUI_Input.output_filename,'Eric_gl',num2str(GUI_Input.glacier_number),'_QGfixed_6h.mat'])
toc



  