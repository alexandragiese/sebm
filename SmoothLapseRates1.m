function [ vLapseRates_smooth ] = SmoothLapseRates1( vLapseRates  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

warning ('off','all');


%% Smooth annual lapse rates

FitOptions = fitoptions('poly7', 'Robust', 'Bisquare');
warning('off','curvefit:fit:equationBadlyConditioned');
LinModel = fit((1:365)',vLapseRates,'poly7',FitOptions); 
vLapseRates_smooth = feval(LinModel,1:365);
vLapseRates_smooth = [vLapseRates_smooth;vLapseRates_smooth(end)]; % Add one day for leap years
vLapseRates_smooth = abs(vLapseRates_smooth);



warning ('on','all');


end

