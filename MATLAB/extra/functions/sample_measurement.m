function [z_k, new_measurement_detected] = sample_measurement(t, patientData)
    z_k = [];
    [new_measurement_detected, idx] = ismember(t, patientData.CGM.time);
    if new_measurement_detected
        z_k = patientData.CGM.values(idx);
    end

    % if randi([0, 5]) == 1 && new_measurement_detected
    %     new_measurement_detected = true;
    % else
    %     new_measurement_detected = false;
    % end
end

