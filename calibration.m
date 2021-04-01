function[R1, GUI_Input] = calibration( GUI_Input,y,g)
PC = GUI_Input.PC(g)
ct = 1; %CAL counter
[R1, ~, ~ ] = Full_Model( GUI_Input,y,g);
geoMB = R1.geodeticMB;
%         cor = GUI_Input.PC(g)
%         mdlMB = R1.modelMB
        
%  if R1.geodeticMB < 0
%     while R1.modelMB > 0.9* R1.geodeticMB
%         ct = ct+1;
%         GUI_Input.PC(g) = GUI_Input.PC(g) - abs( (R1.geodeticMB-R1.modelMB)/ct );
%         [R1, ~, ~ ] = Full_Model( GUI_Input,y,g);
% %         corG = GUI_Input.PC(g)
% %         mdlMBG = R1.modelMB
%     end
%     while R1.modelMB < 1.1* R1.geodeticMB
%         ct = ct+1;
%         GUI_Input.PC(g) = GUI_Input.PC(g) + abs( (R1.geodeticMB-R1.modelMB)/ct );
%         [R1, ~, ~ ] = Full_Model( GUI_Input,y,g);
% %         corL = GUI_Input.PC(g)
% %         mdlMBL = R1.modelMB
%     end
%     while R1.modelMB > 0.9* R1.geodeticMB
%         ct = ct+1;
%         GUI_Input.PC(g) = GUI_Input.PC(g) - abs( (R1.geodeticMB-R1.modelMB)/ct );
%         [R1, ~, ~ ] = Full_Model( GUI_Input,y,g);
% %         corG = GUI_Input.PC(g)
% %         mdlMBG = R1.modelMB
%     end
% elseif R1.geodeticMB > 0
%     while R1.modelMB < 0.9* R1.geodeticMB
%         ct = ct+1;
%         GUI_Input.PC(g) = GUI_Input.PC(g) + abs( (R1.geodeticMB-R1.modelMB)/ct );
%         [R1, ~, ~ ] = Full_Model( GUI_Input,y,g);
%     end
%     while R1.modelMB > 1.1* R1.geodeticMB
%         ct = ct+1;
%         GUI_Input.PC(g) = GUI_Input.PC(g) - abs( (R1.geodeticMB-R1.modelMB)/(i+1) );
%         [R1, ~, ~ ] = Full_Model( GUI_Input,y,g);
%     end
%     while R1.modelMB < 0.9* R1.geodeticMB
%         ct = ct+1;
%         GUI_Input.PC(g) = GUI_Input.PC(g) + abs( (R1.geodeticMB-R1.modelMB)/ct );
%         [R1, ~, ~ ] = Full_Model( GUI_Input,y,g);
%     end
% % else
%     % put something in for glacier 58! with geoMB == 0
%  end

while R1.geodeticMB < 0
    if R1.modelMB > 0.9* R1.geodeticMB
        ct = ct+1;
        GUI_Input.PC(g) = GUI_Input.PC(g) - abs( (R1.geodeticMB-R1.modelMB)/ct );
        [R1, ~, ~ ] = Full_Model( GUI_Input,y,g);  
    elseif R1.modelMB < 1.1* R1.geodeticMB
        ct = ct+1;
        GUI_Input.PC(g) = GUI_Input.PC(g) + abs( (R1.geodeticMB-R1.modelMB)/ct );
        [R1, ~, ~ ] = Full_Model( GUI_Input,y,g); 
    elseif R1.modelMB >= 1.1* R1.geodeticMB && R1.modelMB <= 0.9* R1.geodeticMB
        break
    end
end
while R1.geodeticMB > 0
    if R1.modelMB < 0.9* R1.geodeticMB
        ct = ct+1;
        GUI_Input.PC(g) = GUI_Input.PC(g) + abs( (R1.geodeticMB-R1.modelMB)/ct );
        [R1, ~, ~ ] = Full_Model( GUI_Input,y,g);
    elseif R1.modelMB > 1.1* R1.geodeticMB
        ct = ct+1;
        GUI_Input.PC(g) = GUI_Input.PC(g) - abs( (R1.geodeticMB-R1.modelMB)/ct );
        [R1, ~, ~ ] = Full_Model( GUI_Input,y,g);     
    elseif R1.modelMB <= 1.1* R1.geodeticMB && R1.modelMB >= 0.9* R1.geodeticMB
        break
    end   
end
while R1.geodeticMB == 0
    if R1.modelMB < -0.0194191
        ct = ct+1;
        GUI_Input.PC(g) = GUI_Input.PC(g) + abs( (R1.geodeticMB-R1.modelMB)/ct );
        [R1, ~, ~ ] = Full_Model( GUI_Input,y,g);       
    elseif R1.modelMB > 0.0194191
        ct = ct+1;
        GUI_Input.PC(g) = GUI_Input.PC(g) - abs( (R1.geodeticMB-R1.modelMB)/ct );
        [R1, ~, ~ ] = Full_Model( GUI_Input,y,g);       
    elseif R1.modelMB >= -0.0194191 && R1.modelMB <= 0.0194191
        break
    end
 end

    tmp_cal_PC = GUI_Input.PC(g); 
    disp(g)
    disp(geoMB)
    tmp_MB_corr = R1.modelMB
load([GUI_Input.cal_filename, 'calibration_march.mat'],'cal_PC','MB_corr')
    cal_PC(g) = tmp_cal_PC;
    MB_corr(g) = tmp_MB_corr;
save([GUI_Input.cal_filename, 'calibration_march.mat'],'cal_PC','MB_corr') %will be missing ones calculated during other calculations if simultaneous
