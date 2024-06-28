close all
clear
clc

rng(42, 'twister');

%% Initialization
run('config.m')

mard_intervals = 0:1:30;
uncertainty_intervals = 0:1:50;

simulate_anomalies = true;
show_confidence_interval = true;
do_chi_sq_test = simulate_anomalies;

% how_often = 0.2;
how_often = 0.05;
alpha = 0.05;
fault_reduce_factor = 4;

if use_tuned_model
    run('param_model.m')
end

uncertainty = 30;
% MARD = 10;

Q = eye(13) * uncertainty;
R = 100;

ekf_dt = 2;

model = non_linear_model(tools);
ekf = ekf(model, tools, params, ekf_dt, Q, R);
ekf.dt = ekf_dt;
% ekf.CGM_MARD = MARD;

anomaly_detector = anomaly_detector(alpha);

[x0, y_minus1] = tools.init_conditions(params);

[X, Y] = meshgrid(mard_intervals, uncertainty_intervals);

results = arrayfun(@(MARD, uncertainty) evaluate_EKF_params(MARD, uncertainty, how_often, fault_reduce_factor, ekf, anomaly_detector, patientData, x0, y_minus1, t_start, t_end, params), X, Y);

Z_sensitivity = reshape([results.sensitivity], size(X));
Z_specificity = reshape([results.specificity], size(X));
Z_youden = reshape([results.youden], size(X));

%% Plotting
close all
clc

font_size = 22;
figureWidth = 5.7;
figureHeight = figureWidth;

% Sensitivity
figure;
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, figureWidth, figureHeight]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figureWidth, figureHeight]);
set(gcf, 'PaperPositionMode', 'auto');
contourf(Z_sensitivity, 'LineColor', 'none');
colormap('jet')
caxis([0 1]);
h = colorbar;
set(h, 'Limits', [0 1]);
xlabel('MARD (R)');
ylabel('Uncertainty (Q)');
title('Sensitivity');
axis square;
set(findall(gcf, '-property', 'FontSize'), 'FontSize', font_size);

% Specificity
figure;
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, figureWidth, figureHeight]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figureWidth, figureHeight]);
set(gcf, 'PaperPositionMode', 'auto');
contourf(Z_specificity, 'LineColor', 'none');
colormap('jet')
caxis([0 1]);
h = colorbar;
set(h, 'Limits', [0 1]);
xlabel('MARD (R)');
ylabel('Uncertainty (Q)');
title('Specificity');
axis square;
set(findall(gcf, '-property', 'FontSize'), 'FontSize', font_size);


% Youden's Index
figure;
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, figureWidth, figureHeight]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figureWidth, figureHeight]);
set(gcf, 'PaperPositionMode', 'auto');
contourf(Z_youden, 'LineColor', 'none');
colormap('jet')
caxis([0 1]);
h = colorbar;
set(h, 'Limits', [0 1]);
xlabel('MARD (R)');
ylabel('Uncertainty (Q)');
title('Youden''s Index');
axis square;
set(findall(gcf, '-property', 'FontSize'), 'FontSize', font_size);


%% Function

function result = evaluate_EKF_params(MARD, uncertainty, how_often, fault_reduce_factor, ekf, anomaly_detector, patientData, x0, y_minus1, t_start, t_end, params)
    
    ekf.CGM_MARD = MARD;
    ekf.Q = eye(13)*uncertainty;

    rng(42, 'twister');

    u0 = [0,0];

    t = t_start;
    x = x0;
    y = y_minus1;
    u = u0;
    last_process_update = t_start;
    P0 = eye(13) * 12000;
    P = P0;

    anomaly_detector.anomalies = [];
    anomaly_detector.time = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');

    used_measurement.z = [];
    used_measurement.t = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');
    
    faulty_measurement.z = [];
    faulty_measurement.t = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');
    
    fault_mean = 0;
    
    while t <= t_end
        
        [z_k, new_measurement_detected, faulty, fault_k] = sample_measurement_with_fault(t, patientData, x, ekf.H, how_often, fault_reduce_factor);
        [u_k, new_input_detected] = sample_input(t, patientData);
        
        if faulty
            fault_mean = fault_mean + fault_k;
            faulty_measurement.z(end+1) = z_k;
            faulty_measurement.t(end+1) = t;
        end
    
        if new_input_detected || new_measurement_detected
    
            dt = convert_to_minutes(t - last_process_update);
            [xp_k, Pp_k, y_kminus1, v_kminus1] = ekf.predict(x, y, u, P, dt, params, true);

            x_current = xp_k;
            P_current = Pp_k;
                        
            if new_measurement_detected
                ekf = ekf.update_sensor_cov(z_k);
                [xm_k, Pm_k, residual_k, innov_cov_k] = ekf.measurement_update(xp_k,Pp_k,z_k);

                anomaly = anomaly_detector.chi_squared_test(residual_k,innov_cov_k);

                if ~anomaly
                    x_current = xm_k;
                    P_current = Pm_k;
                else
                    anomaly_detector.time(end+1, :) = t;
                    anomaly_detector.anomalies(end+1) = z_k;
                end

                used_measurement.z(end+1) = z_k;
                used_measurement.t(end+1) = t;
            end
    
            x = x_current;
            P = P_current;
            y = y_kminus1;
            u = u_k;
                
            last_process_update = t;
    
        end
           
        t = next_step(t, patientData);
            
    end
    fault_mean = fault_mean / length(faulty_measurement.z);

    [sensitivity, specificity] = compute_sensitivity_specificity(faulty_measurement, patientData, anomaly_detector);
    
    youden = sensitivity + specificity - 1;
    
    result.sensitivity = sensitivity;
    result.specificity = specificity;
    result.youden = youden;
end

function [sensitivity, specificity] = compute_sensitivity_specificity(faulty_measurement, patientData, anomaly_detector)
    % Initialize counts
    TP = 0; % True Positives
    FP = 0; % False Positives
    TN = 0; % True Negatives
    FN = 0; % False Negatives

    % Get the unique times for comparison
    all_times = unique([faulty_measurement.t'; patientData.CGM.time; anomaly_detector.time]);

    % Loop through each time and classify
    for i = 1:length(all_times)
        t = all_times(i);

        % Check if time t is in faulty_measurement
        is_faulty = ismember(t, faulty_measurement.t);

        % Check if time t is in anomaly_detector
        is_anomaly = ismember(t, anomaly_detector.time);

        % Check if time t is in original dataset
        is_in_data = ismember(t, patientData.CGM.time);

        if is_in_data
            if is_faulty && is_anomaly
                TP = TP + 1;
            elseif ~is_faulty && is_anomaly
                FP = FP + 1;
            elseif ~is_faulty && ~is_anomaly
                TN = TN + 1;
            elseif is_faulty && ~is_anomaly
                FN = FN + 1;
            end
        end
    end

    % Calculate sensitivity and specificity
    sensitivity = TP / (TP + FN);
    specificity = TN / (TN + FP);

    % Display results
    % disp(['Sensitivity: ', num2str(sensitivity)]);
    % disp(['Specificity: ', num2str(specificity)]);
end
