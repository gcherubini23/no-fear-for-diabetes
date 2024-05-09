function next_t = next_step(t, patientData)
    idx_CGM = find(patientData.CGM.time > t);
    idx_IIR = find(patientData.IIR.time > t);
    idx_Meal = find(patientData.Meal.time > t);

    next_t = min([min(patientData.CGM.time(idx_CGM)), min(patientData.Meal.time(idx_Meal)), min(patientData.IIR.time(idx_IIR))]);
    
    if isempty(next_t)
        next_t = t + seconds(1);
    end

end

