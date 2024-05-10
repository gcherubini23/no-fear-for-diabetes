function [u_k, new_input_detected] = sample_input(t, patientData)
    u_k = [0,0];
    [new_meal_detected, idx1] = ismember(t, patientData.Meal.time);
    [new_IIR_detected, idx2] = ismember(t, patientData.IIR.time);
    if new_meal_detected && ~isnan(patientData.Meal.values(idx1))
        u_k(1) = patientData.Meal.values(idx1);
    end
    if new_IIR_detected && ~isnan(patientData.IIR.values(idx2))
        u_k(2) = patientData.IIR.values(idx2);
    end
    if new_meal_detected || new_IIR_detected
        new_input_detected = true;
    else
        new_input_detected = false;
    end
end
