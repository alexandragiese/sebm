%Randomly designate subsets of glaciers for which to run SEBM
clear
subbasin_sub = nan(250,11);
load glnum_bybasin.mat %glacier numbers by basin (col)
R = csvread('DS_wRGI.csv'); %attributes
    RGIId   = R(:,1);
    area_m2 = R(:,13);
    perc_debris = R(:,22);
    CenLon = R(:,28);
    z_med  = R(:,4);
    slope  = R(:,36);
    aspect = R(:,37);
    basin  = R(:,49);
    
for l=1:11; 
C = ~isnan(basin_shean_UIB(:,l)); %glacier IDs 1:1211 (Indus1)
subbasin = basin_shean_UIB(C,l);

k=1;
for j = 1:length(subbasin)
    A(j) = find(RGIId == subbasin(j)); %
    if (area_m2(A(j)) > 100000) && (area_m2(A(j)) < 5000000) && (perc_debris(A(j))<15) && (CenLon(A(j)) > 71.5) 
%     if (area_m2(A(j)) > 100000) && (area_m2(A(j)) < 5000000) && (perc_debris(A(j))<20) && (CenLon(A(j)) > 71) %<-- change restrictions to get 30 in each!!!
        subbasin_sub(k,l) = subbasin(j); %subbasin subset
        k=k+1;
    end
end
% subbasin_m = R(A,:);  
end
% --> make a 30 x 11 table testing_subset.csv
writematrix(subbasin_sub,'clean_median.csv')

glac_nums = subbasin_sub;
save ('clean_median.mat','glac_nums')
