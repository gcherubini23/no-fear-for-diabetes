function [z_k, new_measurement_detected] = sample_BG(t, patientData)
    z_k = [];
    [new_measurement_detected, idx] = ismember(t, patientData.BG.time);
    if new_measurement_detected
        z_k = patientData.BG.values(idx);
    end
end
