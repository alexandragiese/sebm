function[R1, GUI_Input] = calibration( GUI_Input,y,g)

ct = 1; %CAL counter
[R1, ~, ~ ] = Full_Model_AG_edits( GUI_Input,y,g);

 if R1.geodeticMB < 0
    ct = ct+1;
    while R1.modelMB > 0.9* R1.geodeticMB
        GUI_Input.PC(g) = GUI_Input.PC(g) - abs( (R1.geodeticMB-R1.modelMB)/ct );
        [R1, ~, ~ ] = Full_Model_AG_edits( GUI_Input,y,g);
    end
    while R1.modelMB < 1.1* R1.geodeticMB
        GUI_Input.PC(g) = GUI_Input.PC(g) + abs( (R1.geodeticMB-R1.modelMB)/ct );
        [R1, ~, ~ ] = Full_Model_AG_edits( GUI_Input,y,g);
    end

elseif R1.geodeticMB > 0
    ct = ct+1;
    while R1.modelMB < 0.9* R1.geodeticMB
        GUI_Input.PC(g) = GUI_Input.PC(g) + abs( (R1.geodeticMB-R1.modelMB)/(i+1) );
        [R1, ~, ~ ] = Full_Model_AG_edits( GUI_Input,y,g);
    end
    while R1.modelMB > 1.1* R1.geodeticMB
        GUI_Input.PC(g) = GUI_Input.PC(g) - abs( (R1.geodeticMB-R1.modelMB)/(i+1) );
        [R1, ~, ~ ] = Full_Model_AG_edits( GUI_Input,y,g);
    end
% else
    % put something in for glacier 58! with geoMB == 0
end
    tmp_cal_PC = GUI_Input.PC(g); 
    tmp_MB_corr = R1.modelMB;
load([GUI_Input.cal_filename, 'calibration.mat'],'cal_PC','MB_corr')
    cal_PC(g) = tmp_cal_PC;
    MB_corr(g) = tmp_MB_corr;
save([GUI_Input.cal_filename, 'calibration.mat'],'cal_PC','MB_corr') %will be missing ones calculated during other calculations if simultaneous
