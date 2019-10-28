tic

GUI_Input.sGCM = 'HAR';
% GUI_Input.start_date = '01/01/1970';
% GUI_Input.end_date = '12/31/2016';
GUI_Input.start_date = '01/01/2001';
GUI_Input.end_date = '12/31/2013';
GUI_Input.output_filename = '/uufs/chpc.utah.edu/common/home/u0929154/Desktop/Glacier_Model/Indus_Glacier_Sensitivity/Results_FLOR_hourly/Results_Glacier_Number_';
% GUI_Input.lapse_rate = 6.5;
GUI_Input.ensemble_number = 8;
GUI_Input.mb_date_1975 = '01/01/1975';
GUI_Input.mb_date_2000 = '12/31/1999';
GUI_Input.precip_threshold = 1e-2; % (m day^-1)
GUI_Input.model_resolution = 1; % Days
GUI_Input.snow_depth = 0;
% GUI_Input.snow_gradient_0 = 0.00025;
% GUI_Input.min_snow_alt = 3000; % (m)
GUI_Input.snow_density = 250;
GUI_Input.temp_forcing = 0;
GUI_Input.lapse_rate = 6.500;
GUI_Input.iZ_0m_snow = 0.001;
GUI_Input.iZ_0m_ice = 0.016;
GUI_Input.iZ_0T_snow = 0.001;
GUI_Input.iZ_0T_ice = 0.004;
GUI_Input.iZ_0q_snow = 0.001;
GUI_Input.iZ_0q_ice = 0.004;
% Glacier number (Chhota Shigri = 14673, Hamtah = 14219)
% GUI_Input.glacier_number = 14673; 


vJosh_MB_Data = load('Josh_Geo_MB_Data1.mat');
vJosh_Object_ID = vJosh_MB_Data.vJosh_Object_ID;
vcJosh_Category = vJosh_MB_Data.vcJosh_Category;
vJosh_PctDeb = vJosh_MB_Data.vJosh_PctDeb;
vJosh_GeoMassBal = vJosh_MB_Data.vJosh_GeoMassBal;

vModel_MB_1975_to_2000 = NaN(length(vJosh_Object_ID),1);
vModel_MB_2000_to_2016 = NaN(length(vJosh_Object_ID),1);
vModel_MB = NaN(length(vJosh_Object_ID),1);

for r = 377:length(vJosh_Object_ID)

%     try
        
        if strcmp(vcJosh_Category{r}(2:end-1),'Clean') && vJosh_PctDeb(r) <= 10
            
%         if vJosh_PctDeb(r) <= 10

%         for g = 2:30
        
%             GUI_Input.ensemble_number = g;
        
            GUI_Input.glacier_number = vJosh_Object_ID(r); % RGI50-14.18948

            [R1, R2, RAve] = Full_Model_FLOR_Hourly1( GUI_Input );
            
            if strcmp(GUI_Input.sGCM,'FLOR')
                
                try
                    vModel_MB_1975_to_2000(r) = (R1.TotalAveAccum + R1.TotalAveMelt + R1.Total) / (2000-1975+1);
                    vModel_MB_2000_to_2016(r) = (R2.TotalAveAccum + R2.TotalAveMelt) / (2016-2000+1);
                    disp(['Glacier ',num2str(r),' of ',num2str(length(vJosh_Object_ID)),': MB1975 = ',num2str(vModel_MB_1975_to_2000(r)),', MB2000 = ',num2str(vModel_MB_2000_to_2016(r)),' (',num2str(vJosh_GeoMassBal(r)),')'])
                catch
                    iYears = str2double(GUI_Input.end_date(7:10)) - str2double(GUI_Input.start_date(7:10)) + 1;
                    vModel_MB(r) = (R2.TotalAveAccum + R2.TotalAveMelt) / iYears;
                    disp(['Glacier ',num2str(r),' of ',num2str(length(vJosh_Object_ID)),': MB = ',num2str(vModel_MB(r)),' (',num2str(vJosh_GeoMassBal(r)),')'])
                end
                
            elseif strcmp(GUI_Input.sGCM,'HAR')
                iYears = str2double(GUI_Input.end_date(7:10)) - str2double(GUI_Input.start_date(7:10)) + 1;
                vModel_MB(r) = (R2.TotalAveAccum + R2.TotalAveMelt + R2.TotalAveSub) / iYears;
                disp(['Glacier ',num2str(r),' of ',num2str(length(vJosh_Object_ID)),': MB = ',num2str(vModel_MB(r)),' (',num2str(vJosh_GeoMassBal(r)),')'])
            end
            
        
            
            
            toc
            
        
%         end
        
        end

%     catch
%        
%         disp(['Error at r = ',num2str(r),', Glacier: ',num2str(vJosh_Object_ID(r))])
%         
%     end

    if rem(r,30)

        save Josh_MB_temp_save1.mat

    end

end




  