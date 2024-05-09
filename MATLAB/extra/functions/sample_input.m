function [u_k, new_input_detected] = sample_input(t, patientData)
    u_k.CHO = 0;
    u_k.IIR = 0;
    [new_meal_detected, idx1] = ismember(t, patientData.Meal.time);
    [new_IIR_detected, idx2] = ismember(t, patientData.IIR.time);
    if new_meal_detected
        u_k.CHO = patientData.Meal.values(idx1);
    end
    if new_IIR_detected
        u_k.IIR = patientData.IIR.values(idx2);
    end
    if new_meal_detected || new_IIR_detected
        new_input_detected = true;
    else
        new_input_detected = false;
    end
end
