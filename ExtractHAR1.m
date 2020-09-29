function [ SW_in, WindSpeed, Temp, Precip, Pressure, Humidity, Emissivity, LW_in, stTime, Data ] = ExtractHAR1( Lat, Long, sRes, sFileDir )
%UNTITLED2 Summary of this function goes here
%   Extract climate data for a single 10km x 10km grid point for use with
%   surface and energy mass balance model

addpath(genpath(sFileDir))
cd(sFileDir)


mLat = double(ncread('har_d10km_h_2d_t2_2001.nc', 'lat'));
% mLat = flipud(mLat');
mLong = double(ncread('har_d10km_h_2d_t2_2001.nc','lon'));
% mLong = flipud(mLong');

vLat = mLat(:);
vLong = mLong(:);

[~, iLatLongIdx] = min((vLat - Lat).^2 + (vLong - Long).^2);
[m_Idx, n_Idx] = ind2sub([270, 180], iLatLongIdx); %[270, 180] is size of mLong & mLat

%% Grid Elevation

ncid = netcdf.open('har_d10km_static_hgt.nc');

mAWS_Alt = double(ncread('har_d10km_static_hgt.nc', 'hgt'));
Data.AWS_Alt = mAWS_Alt(m_Idx, n_Idx);

netcdf.close(ncid)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hourly Data

if strcmp(sRes,'hourly')

    %% LW_in (W m^-2)

    LW_in = [];

    for i = 0:14

        if i < 10
            sFileName = ['har_d10km_h_2d_lwdown_200', num2str(i), '_rechunked.nc'];
        else
            sFileName = ['har_d10km_h_2d_lwdown_20', num2str(i), '_rechunked.nc'];
        end

        ncid = netcdf.open(sFileName);

            if rem(i,4) == 0
                Leap = 1;
            else
                Leap = 0;
            end
        vLW_in_Data = double(squeeze(ncread(sFileName, 'lwdown',[m_Idx,n_Idx,1], [1,1,8760+24*Leap], [1,1,1])));
    %         SW_in_temp = mSW_in_Data(m_Idx, n_Idx, :);
    %         SW_in_temp = double(SW_in_temp(:));
        LW_in = [LW_in;vLW_in_Data];

        netcdf.close(ncid)

    end


    %% LW_out (W m^-2)
    % 
    % for i = 0:14
    %     
    %     if i < 10
    %         sFileName = ['har_d10km_d_2d_lwup_200', num2str(i), '.nc'];
    %     else
    %         sFileName = ['har_d10km_d_2d_lwup_20', num2str(i), '.nc'];
    %     end
    %     
    %     mLW_out_Data = ncread(sFileName, 'lwup');
    %     LW_out_temp = mLW_out_Data(m_Idx, n_Idx, :);
    %     LW_out_temp = double(LW_out_temp(:));
    %     LW_out(i*365+ceil((i)/4)+1:(i+1)*365+ceil((i+1)/4)) = LW_out_temp;
    % 
    % end

    %% SW_in (W m^-2)

    SW_in = [];

    for i = 0:14

        if i < 10
            sFileName = ['har_d10km_h_2d_swdown_200', num2str(i), '_rechunked.nc'];
        else
            sFileName = ['har_d10km_h_2d_swdown_20', num2str(i), '_rechunked.nc'];
        end

        ncid = netcdf.open(sFileName);

            if rem(i,4) == 0
                Leap = 1;
            else
                Leap = 0;
            end
            vSW_in_Data = double(squeeze(ncread(sFileName, 'swdown',[m_Idx,n_Idx,1], [1,1,8760+24*Leap], [1,1,1])));
    %         SW_in_temp = mSW_in_Data(m_Idx, n_Idx, :);
    %         SW_in_temp = double(SW_in_temp(:));
            SW_in = [SW_in;vSW_in_Data];

        netcdf.close(ncid)

    end

    %% WindSpeed @ 10m (m s^-1)

    WindSpeed = [];

    for i = 0:14

        if i < 10
            sFileName = ['har_d10km_h_2d_ws10_200', num2str(i), '_rechunked.nc'];
        else
            sFileName = ['har_d10km_h_2d_ws10_20', num2str(i), '_rechunked.nc'];
        end

        ncid = netcdf.open(sFileName);

        if rem(i,4) == 0
            Leap = 1;
        else
            Leap = 0;
        end

        vWindSpeed_Data = double(squeeze(ncread(sFileName, 'ws10',[m_Idx,n_Idx,1], [1,1,8760+24*Leap], [1,1,1])));
    %     WindSpeed_temp = mWindSpeed_Data(m_Idx, n_Idx, :);
    %     WindSpeed_temp = double(WindSpeed_temp(:));
        WindSpeed = [WindSpeed;vWindSpeed_Data];

        netcdf.close(ncid)

    end

    %% Temperature (C)

    Temp = [];

    for i = 0:14

        if i < 10
            sFileName = ['har_d10km_h_2d_t2_200', num2str(i), '_rechunked.nc'];
        else
            sFileName = ['har_d10km_h_2d_t2_20', num2str(i), '_rechunked.nc'];
        end

        ncid = netcdf.open(sFileName);

        if rem(i,4) == 0
            Leap = 1;
        else
            Leap = 0;
        end

        vTemp_Data = double(squeeze(ncread(sFileName, 't2',[m_Idx,n_Idx,1], [1,1,8760+24*Leap], [1,1,1])));
    %     Temp_temp = mTemp_Data(m_Idx, n_Idx, :);
    %     Temp_temp = double(Temp_temp(:));
        Temp = [Temp; vTemp_Data];

        netcdf.close(ncid)

    end

    % Convert Kelvin to Celsius
    Temp = Temp - 273.15;

    %% Precipitation (m/h)

    Precip = [];

    for i = 0:14

        if i < 10
            sFileName = ['har_d10km_h_2d_prcp_200', num2str(i), '_rechunked.nc'];
        else
            sFileName = ['har_d10km_h_2d_prcp_20', num2str(i), '_rechunked.nc'];
        end

        ncid = netcdf.open(sFileName);

        if rem(i,4) == 0
            Leap = 1;
        else
            Leap = 0;
        end


        vPrecip_Data = double(squeeze(ncread(sFileName, 'prcp',[m_Idx,n_Idx,1], [1,1,8760+24*Leap], [1,1,1])));
    %     Precip_temp = mPrecip_Data(m_Idx, n_Idx, :);
    %     Precip_temp = double(Precip_temp(:));
        Precip = [Precip;vPrecip_Data];

        netcdf.close(ncid)

    end

    % Convert to m / hr
    Precip = Precip / 1000;

    %% Pressure (hPa)

    Pressure = [];

    for i = 0:14

        if i < 10
            sFileName = ['har_d10km_h_2d_psfc_200', num2str(i), '_rechunked.nc'];
        else
            sFileName = ['har_d10km_h_2d_psfc_20', num2str(i), '_rechunked.nc'];
        end

        ncid = netcdf.open(sFileName);

        if rem(i,4) == 0
            Leap = 1;
        else
            Leap = 0;
        end

        vPressure_Data = double(squeeze(ncread(sFileName, 'psfc',[m_Idx,n_Idx,1], [1,1,8760+24*Leap], [1,1,1]))) / 100;
    %     Pressure_temp = mPressure_Data(m_Idx, n_Idx, :);
    %     Pressure_temp = double(Pressure_temp(:)) / 100;
        Pressure = [Pressure;vPressure_Data];

        netcdf.close(ncid)

    end

    %% Specific Humidity (kg kg^-1)

    SpecHumidity = [];

    for i = 0:14

        if i < 10
            sFileName = ['har_d10km_h_2d_q2_200', num2str(i), '_rechunked.nc'];
        else
            sFileName = ['har_d10km_h_2d_q2_20', num2str(i), '_rechunked.nc'];
        end

        ncid = netcdf.open(sFileName);

        if rem(i,4) == 0
            Leap = 1;
        else
            Leap = 0;
        end

        vSpecHum_Data = double(squeeze(ncread(sFileName, 'q2',[m_Idx,n_Idx,1], [1,1,8760+24*Leap], [1,1,1])));
    %     SpecHum_temp = mSpecHum_Data(m_Idx, n_Idx, :);
    %     SpecHum_temp = double(SpecHum_temp(:));
        SpecHumidity = [SpecHumidity;vSpecHum_Data];

        netcdf.close(ncid)

    end


    %% Cloudiness (fraction)

    % Cloudiness = [];
    % 
    % for i = 0:14
    %     
    %     if i < 10
    %         sFileName = ['har_d10km_h_2d_scldfra_200', num2str(i), '.nc'];
    %     else
    %         sFileName = ['har_d10km_h_2d_scldfra_20', num2str(i), '.nc'];
    %     end
    %     
    %     mCloudiness_Data = ncread(sFileName, 'scldfra');
    %     Cloudiness_temp = mCloudiness_Data(m_Idx, n_Idx, :);
    %     Cloudiness_temp = double(Cloudiness_temp(:));
    %     Cloudiness = [Cloudiness;Cloudiness_temp];
    % 
    % end

    %% Emissivity

    Emissivity = [];

    for i = 0:14

        if i < 10
            sFileName = ['har_d10km_h_2d_emiss_200', num2str(i), '_rechunked.nc'];
        else
            sFileName = ['har_d10km_h_2d_emiss_20', num2str(i), '_rechunked.nc'];
        end

        ncid = netcdf.open(sFileName);

        if rem(i,4) == 0
            Leap = 1;
        else
            Leap = 0;
        end

        vEmissivity_Data = double(squeeze(ncread(sFileName, 'emiss',[m_Idx,n_Idx,1], [1,1,8760+24*Leap], [1,1,1])));
    %     Emissivity_temp = mEmissivity_Data(m_Idx, n_Idx, :);
    %     Emissivity_temp = double(Emissivity_temp(:));
        Emissivity = [Emissivity;vEmissivity_Data];

        netcdf.close(ncid)

    end


    %% Humidity (%)

    % vSVP = 611.2 .* exp(17.67 * (Temp) ./ (Temp - 29.65));                      % http://www.atmos.washington.edu/epic/soundings/qc_description.html
    % Humidity = 100 * SpecHumidity ./ (0.622 * (vSVP ./ Pressure));               % http://www.atmos.washington.edu/epic/soundings/qc_description.html
    % Humidity(Humidity > 100) = 100;

    Humidity = 100 * 0.263 * Pressure .* SpecHumidity .* (exp(17.67 * Temp ...
        ./ (Temp + 273 - 29.65))) .^ -1;                                        % http://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
    Humidity(Humidity > 100) = 100;

    %% Other Data

    % Temporal resolution of the raw data (minutes)
    Data.Resolution_mins = 60;
    Data.AWS_Lat = Lat;
    Data.AWS_Long = Long;
    Data.TimeZone = 0;
    Data.FirstDay = 0;

    %% Calculate time vectors

    stTime.year = [ones(8784,1)*2000; ones(8760,1)*2001; ones(8760,1)*2002;...
        ones(8760,1)*2003;ones(8784,1)*2004;ones(8760,1)*2005;ones(8760,1)*2006;...
        ones(8760,1)*2007;ones(8784,1)*2008;ones(8760,1)*2009;ones(8760,1)*2010;...
        ones(8760,1)*2011;ones(8784,1)*2012;ones(8760,1)*2013;ones(8760,1)*2014];
    stTime.month = repmat([repmat(1,[31,1]);repmat(2,[28,1]);repmat(3,[31,1]);repmat(4,[30,1]);...
        repmat(5,[31,1]);repmat(6,[30,1]);repmat(7,[31,1]);repmat(8,[31,1]);...
        repmat(9,[30,1]);repmat(10,[31,1]);repmat(11,[30,1]);repmat(12,[31,1])],15,1);
    stTime.month = repelem(stTime.month,24);
    stTime.month = [stTime.month(1:1416);repmat(2,[24,1]);...
        stTime.month(1416+1:4*365*24+1416);repmat(2,[24,1]);...
        stTime.month(4*365*24+1416+1:8*365*24+1416);repmat(2,[24,1]);...
        stTime.month(8*365*24+1416+1:12*365*24+1416);repmat(2,[24,1]);...
        stTime.month(12*365*24+1416+1:end)];
    stTime.day = [1:31, 1:28, 1:31, 1:30, 1:31, 1:30, 1:31, 1:31, 1:30, 1:31, 1:30, 1:31]';
    stTime.day = repelem(stTime.day,24);
    stTime.day = repmat(stTime.day,[15,1]);
    stTime.day = [stTime.day(1:1416);repmat(29,[24,1]);...
        stTime.day(1416+1:4*365*24+1416);repmat(29,[24,1]);...
        stTime.day(4*365*24+1416+1:8*365*24+1416);repmat(29,[24,1]);...
        stTime.day(8*365*24+1416+1:12*365*24+1416);repmat(29,[24,1]);...
        stTime.day(12*365*24+1416+1:end)];
    stTime.hour = repmat((1:24)',length(stTime.year)/24,1);
    stTime.min = zeros(size(stTime.year));
    stTime.sec = zeros(size(stTime.year));
    stTime.UTC = 5;
    
    

elseif strcmp(sRes,'daily')
    
    %% LW_in (W m^-2)

    LW_in = [];

    for i = 0:14

        if i < 10
            sFileName = ['har_d10km_d_2d_lwdown_200', num2str(i), '_rechunked.nc'];
        else
            sFileName = ['har_d10km_d_2d_lwdown_20', num2str(i), '_rechunked.nc'];
        end

        ncid = netcdf.open(sFileName);

            if rem(i,4) == 0
                Leap = 1;
            else
                Leap = 0;
            end
            
        vLW_in_Data = double(squeeze(ncread(sFileName, 'lwdown',[m_Idx,n_Idx,1], [1,1,365+Leap], [1,1,1])));
        LW_in = [LW_in;vLW_in_Data];

        netcdf.close(ncid)

    end

    %% SW_in (W m^-2)

    SW_in = [];

    for i = 0:14

        if i < 10
            sFileName = ['har_d10km_d_2d_swdown_200', num2str(i), '_rechunked.nc'];
        else
            sFileName = ['har_d10km_d_2d_swdown_20', num2str(i), '_rechunked.nc'];
        end

        ncid = netcdf.open(sFileName);

            if rem(i,4) == 0
                Leap = 1;
            else
                Leap = 0;
            end
            
            vSW_in_Data = double(squeeze(ncread(sFileName, 'swdown',[m_Idx,n_Idx,1], [1,1,365+Leap], [1,1,1])));
            SW_in = [SW_in;vSW_in_Data];

        netcdf.close(ncid)

    end

    %% WindSpeed @ 10m (m s^-1)

    WindSpeed = [];

    for i = 0:14

        if i < 10
            sFileName = ['har_d10km_d_2d_ws10_200', num2str(i), '_rechunked.nc'];
        else
            sFileName = ['har_d10km_d_2d_ws10_20', num2str(i), '_rechunked.nc'];
        end

        ncid = netcdf.open(sFileName);

        if rem(i,4) == 0
            Leap = 1;
        else
            Leap = 0;
        end

        vWindSpeed_Data = double(squeeze(ncread(sFileName, 'ws10',[m_Idx,n_Idx,1], [1,1,365+Leap], [1,1,1])));
        WindSpeed = [WindSpeed;vWindSpeed_Data];

        netcdf.close(ncid)

    end

    %% Temperature (C)

    Temp = [];

    for i = 0:14

        if i < 10
            sFileName = ['har_d10km_d_2d_t2_200', num2str(i), '_rechunked.nc'];
        else
            sFileName = ['har_d10km_d_2d_t2_20', num2str(i), '_rechunked.nc'];
        end

        ncid = netcdf.open(sFileName);

        if rem(i,4) == 0
            Leap = 1;
        else
            Leap = 0;
        end

        vTemp_Data = double(squeeze(ncread(sFileName, 't2',[m_Idx,n_Idx,1], [1,1,365+Leap], [1,1,1])));
        Temp = [Temp; vTemp_Data];

        netcdf.close(ncid)

    end

    % Convert Kelvin to Celsius
    Temp = Temp - 273.15;
    %% Precipitation (m/day)

    Precip = [];

    for i = 0:14

        if i < 10
            sFileName = ['har_d10km_d_2d_prcp_200', num2str(i), '_rechunked.nc'];
        else
            sFileName = ['har_d10km_d_2d_prcp_20', num2str(i), '_rechunked.nc'];
        end

        ncid = netcdf.open(sFileName);

        if rem(i,4) == 0
            Leap = 1;
        else
            Leap = 0;
        end


        vPrecip_Data = double(squeeze(ncread(sFileName, 'prcp',[m_Idx,n_Idx,1], [1,1,365+Leap], [1,1,1])));
        Precip = [Precip;vPrecip_Data];

        netcdf.close(ncid)

    end

    % Convert to m / hr
    Precip = Precip * 24 / 1000; %24h is necessary b/c HAR gives hourly avg. precip for the day --> calculate a daily total


    %% Pressure (hPa)

    Pressure = [];

    for i = 0:14

        if i < 10
            sFileName = ['har_d10km_d_2d_psfc_200', num2str(i), '_rechunked.nc'];
        else
            sFileName = ['har_d10km_d_2d_psfc_20', num2str(i), '_rechunked.nc'];
        end

        ncid = netcdf.open(sFileName);

        if rem(i,4) == 0
            Leap = 1;
        else
            Leap = 0;
        end

        vPressure_Data = double(squeeze(ncread(sFileName, 'psfc',[m_Idx,n_Idx,1], [1,1,365+Leap], [1,1,1]))) / 100;
        Pressure = [Pressure;vPressure_Data];

        netcdf.close(ncid)

    end

    %% Specific Humidity (kg kg^-1)

    SpecHumidity = [];

    for i = 0:14

        if i < 10
            sFileName = ['har_d10km_d_2d_q2_200', num2str(i), '_rechunked.nc'];
        else
            sFileName = ['har_d10km_d_2d_q2_20', num2str(i), '_rechunked.nc'];
        end

        ncid = netcdf.open(sFileName);

        if rem(i,4) == 0
            Leap = 1;
        else
            Leap = 0;
        end

        vSpecHum_Data = double(squeeze(ncread(sFileName, 'q2',[m_Idx,n_Idx,1], [1,1,365+Leap], [1,1,1])));
        SpecHumidity = [SpecHumidity;vSpecHum_Data];

        netcdf.close(ncid)

    end

    %% Emissivity

    Emissivity = [];

    for i = 0:14

        if i < 10
            sFileName = ['har_d10km_d_2d_emiss_200', num2str(i), '_rechunked.nc'];
        else
            sFileName = ['har_d10km_d_2d_emiss_20', num2str(i), '_rechunked.nc'];
        end

        ncid = netcdf.open(sFileName);

        if rem(i,4) == 0
            Leap = 1;
        else
            Leap = 0;
        end

        vEmissivity_Data = double(squeeze(ncread(sFileName, 'emiss',[m_Idx,n_Idx,1], [1,1,365+Leap], [1,1,1])));
        Emissivity = [Emissivity;vEmissivity_Data];

        netcdf.close(ncid)

    end

    %% Humidity (%)

    % vSVP = 611.2 .* exp(17.67 * (Temp) ./ (Temp - 29.65));                      % http://www.atmos.washington.edu/epic/soundings/qc_description.html
    % Humidity = 100 * SpecHumidity ./ (0.622 * (vSVP ./ Pressure));               % http://www.atmos.washington.edu/epic/soundings/qc_description.html
    % Humidity(Humidity > 100) = 100;

    Humidity = 100 * 0.263 * Pressure .* SpecHumidity .* (exp(17.67 * Temp ...
        ./ (Temp + 273 - 29.65))) .^ -1;                                        % http://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
    Humidity(Humidity > 100) = 100;

    %% Other Data

    % Temporal resolution of the raw data (minutes)
    Data.Resolution_mins = 1440;
    Data.AWS_Lat = Lat;
    Data.AWS_Long = Long;
    Data.TimeZone = 0;
    Data.FirstDay = 0;

    %% Calculate time vectors

    stTime.year = zeros(length(Temp),1);
    stTime.month = zeros(length(Temp),1);
    stTime.day = zeros(length(Temp),1);

    kYearLength = 365;

    for iYear = 2000:2014

        if rem(iYear,4) > 0
            iLeapYear = 0;
        elseif rem(iYear,4) == 0
            iLeapYear = 1;
        end

        iYearLength = kYearLength + iLeapYear;
        iLast_idx = find(stTime.year,1,'last');

        if isempty(iLast_idx)
            iLast_idx = 0;
        end

        stTime.year(iLast_idx+1:iLast_idx+iYearLength) = iYear;
        vMonth_temp = [repmat(1,[31,1]);repmat(2,[28+iLeapYear,1]);repmat(3,[31,1]);repmat(4,[30,1]);...
            repmat(5,[31,1]);repmat(6,[30,1]);repmat(7,[31,1]);repmat(8,[31,1]);...
            repmat(9,[30,1]);repmat(10,[31,1]);repmat(11,[30,1]);repmat(12,[31,1])];
        stTime.month(iLast_idx+1:iLast_idx+iYearLength) = vMonth_temp;
        vDay_temp = [1:31, 1:28+iLeapYear, 1:31, 1:30, 1:31, 1:30, 1:31, 1:31, 1:30, 1:31, 1:30, 1:31]';
        stTime.day(iLast_idx+1:iLast_idx+iYearLength) = vDay_temp;

    end

    stTime.hour = zeros(size(stTime.year));
    stTime.min = zeros(size(stTime.year));
    stTime.sec = zeros(size(stTime.year));
    stTime.UTC = 0;
    
end








































































end

