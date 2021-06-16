function[R1, GUI_Input] = Extremes_calibration( GUI_Input,y,g,k)
iGlacierNumber = y(g);
GUI_Input.PC = 0;

%RENAME MINIMUM EXTREMES before starting max!!!!!!!!!!!!!!!!!!!!!!!
ext = 1; %2; %1) minimum or 2) maximum DESIGNATE THE EXTREME

% Constant Parameters: full UIB run
GUI_Input.Alpha_i    = 0.34;            
GUI_Input.Alpha_fs   = 0.75;            
GUI_Input.RainThreshold = 2;
GUI_Input.lra = 0; %deg/km (-1 --> 1)
GUI_Input.Z_0m_snow  = 0.001;
GUI_Input.Z_0m_ice   = 0.016;
GUI_Input.Z_0T_snow  = 0.001;
GUI_Input.Z_0T_ice   = 0.004;
GUI_Input.Z_0q_snow  = 0.001;
GUI_Input.Z_0q_ice   = 0.004;

% Constant Parameters: MINIMUM MELT
if ext ==1
    GUI_Input.Alpha_i    = 0.4;            
    GUI_Input.Alpha_fs   = 0.85;            
    GUI_Input.RainThreshold = 0;
    GUI_Input.lra = 1.1;
    GUI_Input.Z_0m_snow  = 0.0005;
    GUI_Input.Z_0m_ice   = 0.00689;
    GUI_Input.Z_0T_snow  = 0.0008;
    GUI_Input.Z_0T_ice   = 0.0038;
    GUI_Input.Z_0q_snow  = 0.0008;
    GUI_Input.Z_0q_ice   = 0.0038;
end
% Constant Parameters: MAXIMUM
if ext ==2
    GUI_Input.Alpha_i    = 0.28;            
    GUI_Input.Alpha_fs   = 0.69;            
    GUI_Input.RainThreshold = 2;
    GUI_Input.lra = 0.9;
    GUI_Input.Z_0m_snow  = 0.0015;
    GUI_Input.Z_0m_ice   = 0.016;
    GUI_Input.Z_0T_snow  = 0.0012;
    GUI_Input.Z_0T_ice   = 0.0042;
    GUI_Input.Z_0q_snow  = 0.0012;
    GUI_Input.Z_0q_ice   = 0.0042;   
end

% Other parameters
ct = 1; %CAL counter
fl = 'cal'; %calibrating
[R1, ~, ~ ] = Extremes_Model( GUI_Input,y,g,ct,fl);
geo_MB = R1.geodeticMB; %disp('^ make sure this is not zero')

while R1.geodeticMB < 0
    if R1.modelMB > 0.9* R1.geodeticMB
        ct = ct+1;
        GUI_Input.PC = GUI_Input.PC - abs( (R1.geodeticMB-R1.modelMB)/ct );
        [R1, ~, ~] = Extremes_Model (GUI_Input,y,g,ct,fl );  
    elseif R1.modelMB < 1.1* R1.geodeticMB
        ct = ct+1;
        GUI_Input.PC = GUI_Input.PC + abs( (R1.geodeticMB-R1.modelMB)/ct );
        [R1, ~, ~] = Extremes_Model (GUI_Input,y,g,ct,fl ); 
    elseif R1.modelMB >= 1.1* R1.geodeticMB && R1.modelMB <= 0.9* R1.geodeticMB
        break
    end
end
while R1.geodeticMB > 0
    if R1.modelMB < 0.9* R1.geodeticMB
        ct = ct+1;
        GUI_Input.PC = GUI_Input.PC + abs( (R1.geodeticMB-R1.modelMB)/ct );
        [R1, ~, ~] = Extremes_Model (GUI_Input,y,g,ct,fl );
    elseif R1.modelMB > 1.1* R1.geodeticMB
        ct = ct+1;
        GUI_Input.PC = GUI_Input.PC - abs( (R1.geodeticMB-R1.modelMB)/ct );
        [R1, ~, ~] = Extremes_Model (GUI_Input,y,g,ct,fl );     
    elseif R1.modelMB <= 1.1* R1.geodeticMB && R1.modelMB >= 0.9* R1.geodeticMB
        break
    end   
end
% while R1.geodeticMB == 0
%     if R1.modelMB < -0.0194191 %This is given specifically for the one calibration subset glacier with a 0 geoMB.  If one of the MC glaciers has a 0 geoMB, rethink this.
%         ct = ct+1;
%         GUI_Input.PC = GUI_Input.PC + abs( (R1.geodeticMB-R1.modelMB)/ct );
%         [R1, ~, ~ ] = Extremes_Model (GUI_Input,y,g,ct,fl );       
%     elseif R1.modelMB > 0.0194191
%         ct = ct+1;
%         GUI_Input.PC = GUI_Input.PC - abs( (R1.geodeticMB-R1.modelMB)/ct );
%         [R1, ~, ~ ] = Extremes_Model (GUI_Input,y,g,ct,fl );       
%     elseif R1.modelMB >= -0.0194191 && R1.modelMB <= 0.0194191
%         break
%     end
% end

tmp_cal_PC   = GUI_Input.PC;
tmp_cal_MB   = R1.modelMB;
tmp_cal_MELT = R1.modelMelt;

% Give a PC from within 1 std of the calibration distribution to use for a final iteration-----------
% load ('/uufs/chpc.utah.edu/common/home/cryosphere/agiese/calibration/calibration_march.mat','cal_PC')
% cal_PC(232) = nan; %(remember to remove glacier 232)
% S = nanstd(cal_PC);
S = 0.5645;
if ext == 1
    s = S; %min 
elseif ext == 2
    s = -S; %max
elseif ext == 0
    s = 0;
end
ct = ct+1;
GUI_Input.PC = tmp_cal_PC + s;
fl = 'sigmaPC';
[R1, ~, ~] = Extremes_Model (GUI_Input,y,g,ct,fl );  

load([GUI_Input.cal_filename,'MinimumExtremes.mat'],'cal_PC','MB_corr','MELT_corr','geoMB','varied_PC','varied_MB','varied_MELT')  
    cal_PC(k)    = tmp_cal_PC;
    MB_corr(k)   = tmp_cal_MB; 
    MELT_corr(k) = tmp_cal_MELT;
    geoMB(k)     = geo_MB;
    varied_PC(k)   = GUI_Input.PC;
    varied_MB(k)   = R1.modelMB;
    varied_MELT(k) = R1.modelMelt;
save([GUI_Input.cal_filename,'MinimumExtremes.mat'],'cal_PC','MB_corr','MELT_corr','geoMB','varied_PC','varied_MB','varied_MELT') %will be missing ones calculated during other calculations if simultaneous