function[R1, GUI_Input] = MCcalibration( GUI_Input,y,g,mc,vP_d_th)

GUI_Input.PC = 0;
mc
% Constant Parameters
GUI_Input.Alpha_i    = 0.34;                                                             % 0.45 from Molg and Hardy, 2004; 0.35 from Physics of Glaciers; 0.30 in Johnson & Rupper
GUI_Input.Alpha_fs   = 0.75;                                                             % 0.9 from Molg and Hardy, 2004; 0.85 from Physics of Glaciers
GUI_Input.RainThreshold = 2;
GUI_Input.lra = 0; %deg/km (-1 --> 1)
GUI_Input.Z_0m_snow  = 0.01; %0.001;
GUI_Input.Z_0m_ice   = 0.018; %0.016;
GUI_Input.Z_0T_snow  = 0.002; %0.001;
GUI_Input.Z_0T_ice   = 0.006; %0.004;
GUI_Input.Z_0q_snow  = 0.002; %0.001;
GUI_Input.Z_0q_ice   = 0.006; %0.004;


% Randomly selected parameters-----------
% pd = makedist('Normal',0.75,0.004);
%     t = truncate(pd,0.6,0.9);
%     GUI_Input.Alpha_fs = random(t);
% GUI_Input.Alpha_i    = random('Normal',0.34,0.025);  
% 
% pd = makedist('Lognormal','mu',log(1.2),'sigma',.4);
%     t = truncate(pd,0,4);
%     GUI_Input.RainThreshold = 3-random(t);                                                          
% pd = makedist('Normal',0,0.3);
%     t = truncate(pd,-1,1);
%     GUI_Input.lra = random(t);
%     
% pd = makedist('Lognormal','mu',log(0.001),'sigma',.55);
%     GUI_Input.Z_0m_snow = random(pd);
% GUI_Input.Z_0T_snow = random('Normal',0.001,.0002) ;
% GUI_Input.Z_0q_snow  = random('Normal',0.001,.0002);
% 
% pd = makedist('Lognormal','mu',log(0.0007),'sigma',.9); 
%     t = truncate(pd,-1,0.0167);
%     GUI_Input.Z_0m_ice = .0167-random(t);
% GUI_Input.Z_0T_ice = random('Normal',0.004,.0002);
% GUI_Input.Z_0q_ice = random('Normal',0.004,.0002);

% Other parameters
ct = 1; %CAL counter
fl = 'cal'; %calibrating
[R1, ~, ~ ] = MC_Model( GUI_Input,vP_d_th,y,g,mc,ct,fl);
geo_MB = R1.geodeticMB, disp('make sure this is not zero')

keyboard
while R1.geodeticMB < 0
    if R1.modelMB > 0.9* R1.geodeticMB
        ct = ct+1;
        GUI_Input.PC = GUI_Input.PC - abs( (R1.geodeticMB-R1.modelMB)/ct );
        [R1, ~, ~, ~] = MC_Model( GUI_Input,vP_d_th,y,g,mc,ct,fl);  
    elseif R1.modelMB < 1.1* R1.geodeticMB
        ct = ct+1;
        GUI_Input.PC = GUI_Input.PC + abs( (R1.geodeticMB-R1.modelMB)/ct );
        [R1, ~, ~, ~ ] = MC_Model( GUI_Input,vP_d_th,y,g,mc,ct,fl); 
    elseif R1.modelMB >= 1.1* R1.geodeticMB && R1.modelMB <= 0.9* R1.geodeticMB
        break
    end
end
while R1.geodeticMB > 0
    if R1.modelMB < 0.9* R1.geodeticMB
        ct = ct+1;
        GUI_Input.PC = GUI_Input.PC + abs( (R1.geodeticMB-R1.modelMB)/ct );
        [R1, ~, ~, ~ ] = MC_Model( GUI_Input,vP_d_th,y,g,mc,ct,fl);
    elseif R1.modelMB > 1.1* R1.geodeticMB
        ct = ct+1;
        GUI_Input.PC = GUI_Input.PC - abs( (R1.geodeticMB-R1.modelMB)/ct );
        [R1, ~, ~, ~ ] = MC_Model( GUI_Input,vP_d_th,y,g,mc,ct,fl);     
    elseif R1.modelMB <= 1.1* R1.geodeticMB && R1.modelMB >= 0.9* R1.geodeticMB
        break
    end   
end
while R1.geodeticMB == 0
    if R1.modelMB < -0.0194191 %This is given specifically for the one calibration subset glacier with a 0 geoMB.  If one of the MC glaciers has a 0 geoMB, rethink this.
        ct = ct+1;
        GUI_Input.PC = GUI_Input.PC + abs( (R1.geodeticMB-R1.modelMB)/ct );
        [R1, ~, ~, ~ ] = MC_Model( GUI_Input,vP_d_th,y,g,mc,ct,fl);       
    elseif R1.modelMB > 0.0194191
        ct = ct+1;
        GUI_Input.PC = GUI_Input.PC - abs( (R1.geodeticMB-R1.modelMB)/ct );
        [R1, ~, ~, ~ ] = MC_Model( GUI_Input,vP_d_th,y,g,mc,ct,fl);       
    elseif R1.modelMB >= -0.0194191 && R1.modelMB <= 0.0194191
        break
    end
end

tmp_cal_PC   = GUI_Input.PC;
tmp_cal_MB   = R1.modelMB;
tmp_cal_MELT = R1.modelMELT;

% Randomy select a PC from within 1 std of the calibration distribution to use for a final MC-----------
% load ('/uufs/chpc.utah.edu/common/home/cryosphere/agiese/calibration/calibration_march.mat','cal_PC')
% cal_PC(232) = nan; %(remember to remove glacier 232)
% S = nanstd(cal_PC);
S = 0.5645;
    pd = makedist('Uniform','Lower',-S,'Upper',S);
    s = random(pd);
fl = 'final';
ct = ct+1;
GUI_Input.PC = tmp_cal_PC + s;
[R1, ~, ~, ~] = MC_Model( GUI_Input,vP_d_th,y,g,mc,ct,fl);  

load([GUI_Input.cal_filename,'Glacier_Number_',num2str(iGlacierNumber),'_MonteCarlo.mat'],'cal_PC','MB_corr','MELT_corr','geoMB','varied_PC','varied_MB','varied_MELT')  
    cal_PC(mc)    = tmp_cal_PC;
    MB_corr(mc)   = tmp_cal_MB; 
    MELT_corr(mc) = tmp_cal_MELT;
    geoMB(mc)     = geo_MB;
    varied_PC(mc)   = GUI_Input.PC;
    varied_MB(mc)   = R1.modelMB;
    varied_MELT(mc) = R1.modelMELT;
save([GUI_Input.cal_filename,'Glacier_Number_',num2str(iGlacierNumber),'MonteCarlo.mat'],'cal_PC','MB_corr','MELT_corr','geoMB','varied_PC','varied_MB','varied_MELT') %will be missing ones calculated during other calculations if simultaneous
%^ assess this file to see whether the mean and std are changing much