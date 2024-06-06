function [z_k, new_measurement_detected, faulty, fault_k] = sample_measurement_with_fault(t, patientData, x, H, how_often, fault_reduce_factor)

    % Initialize output variables
    z_k = [];
    faulty = false;
    fault_k = 0;

    % Check if the current time t is in the patientData CGM time
    [new_measurement_detected, idx] = ismember(t, patientData.CGM.time);

    if new_measurement_detected
        % Get the nominal measurement value
        z_k = patientData.CGM.values(idx);

        z = binornd(1, how_often);
        if z == 1
            Ul = 400;
            % Get the PG value at time t (assuming PG is stored similarly to CGM)
            % PG_k = H * x';
            PG_k = z_k;
    
            % Get the nominal CGM value at time t
            CGM_nominal_k = z_k;
    
            % Calculate the lower and upper bounds for the fault
            lower_bound = PG_k + 110 - CGM_nominal_k;
            upper_bound = Ul - CGM_nominal_k;
    
            % Generate a random fault value from the uniform distribution
            fault_k = unifrnd(lower_bound, upper_bound) / fault_reduce_factor;
    
            % Apply the fault to the measurement
            z_k = z_k + fault_k;
            
            % Indicate that a faulty measurement has been generated
            faulty = true;
        end
    end
end