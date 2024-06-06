close all
clear
clc

% rng default;
rng(42, 'twister');

%% Initialization
run('config.m')

simulate_anomalies = true;
show_confidence_interval = true;
do_chi_sq_test = simulate_anomalies;

how_often = 0.05;
alpha = 0.05;
fault_reduce_factor = 4;

MARD = 7;
uncertainty = 30;

if use_tuned_model
    run('param_model.m')
end


Q = eye(length(state_fields)) * uncertainty;
R = 100;

ekf_dt = 2.3;

model = non_linear_model(tools);
ekf = ekf(model, tools, params, ekf_dt, Q, R);
ekf.dt = ekf_dt;
ekf.CGM_MARD = MARD;

anomaly_detector = anomaly_detector(alpha);


[x0, y_minus1] = tools.init_conditions(params);
u0 = [0,0];

EKF_state_tracking.mean = [];
EKF_state_tracking.variance = [];
EKF_state_tracking.time = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');
residuals = [];

t = t_start;
x = x0;
y = y_minus1;
u = u0;
last_process_update = t_start;
P0 = eye(length(state_fields)) * 12000;
% P0 = Q * 6000;
P = P0;

used_measurement.z = [];
used_measurement.t = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');

faulty_measurement.z = [];
faulty_measurement.t = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');

fault_mean = 0;

%% Start simulation

disp('Starting simulation...')
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
        if plot_complete_history
            [xp_k, Pp_k, y_kminus1, v_kminus1, ekf] = ekf.predict_save_all_variables(x, y, u, P, dt, params, last_process_update, true);
        else
            [xp_k, Pp_k, y_kminus1, v_kminus1] = ekf.predict(x, y, u, P, dt, params, true);
        end
        x_current = xp_k;
        P_current = Pp_k;
                    
        if new_measurement_detected && do_measurment_update
            ekf = ekf.update_sensor_cov(z_k);
            [xm_k, Pm_k, residual_k, innov_cov_k] = ekf.measurement_update(xp_k,Pp_k,z_k);
            residuals(end+1) = residual_k;
            if do_chi_sq_test || do_cusum_test
                if do_chi_sq_test
                    anomaly = anomaly_detector.chi_squared_test(residual_k,innov_cov_k);
                else
                    [anomaly, anomaly_detector] = anomaly_detector.cusum_test(residual_k,innov_cov_k);
                end
                if ~anomaly
                    x_current = xm_k;
                    P_current = Pm_k;
                else
                    anomaly_detector.time(end+1, :) = t;
                    anomaly_detector.anomalies(end+1) = z_k;
                end
            else
                x_current = xm_k;
                P_current = Pm_k;
                anomaly = false;
            end
            used_measurement.z(end+1) = z_k;
            used_measurement.t(end+1) = t;
        end

        x = x_current;
        P = P_current;
        y = y_kminus1;
        u = u_k;
            
        last_process_update = t;

        EKF_state_tracking.mean(:,end+1) = x';
        EKF_state_tracking.variance(end+1) = ekf.H * P * ekf.H';
        EKF_state_tracking.time(end+1, :) = t;

    end
       
    t = next_step(t, patientData);
        
end

fault_mean = fault_mean / length(faulty_measurement.z)
[sensitivity, specificity] = compute_sensitivity_specificity(faulty_measurement, patientData, anomaly_detector);

%% Plot
close all

figure
ms = 6;

name = 'Tracking \pm 2\sigma';

plot(EKF_state_tracking.time, ekf.H * EKF_state_tracking.mean, '-', 'DisplayName', name, 'LineWidth', 2, 'Color', 'm');
hold on

% plot(patientData.CGM.time, patientData.CGM.values, '-o', 'DisplayName', 'CGM', 'Color', 'cyan', 'LineWidth', 1, 'MarkerSize', ms)
plot(used_measurement.t, used_measurement.z, '-o', 'DisplayName', 'CGM', 'Color', 'cyan', 'LineWidth', 1, 'MarkerSize', ms)
hold on

plot(faulty_measurement.t, faulty_measurement.z, 'o', 'DisplayName', 'Anomaly not detected', 'Color', 'blue', 'MarkerSize', ms, 'MarkerFaceColor',[0 0 1])
hold on

if show_confidence_interval
    gamma = 2;
    sigma = sqrt(EKF_state_tracking.variance);
    upper_bound = ekf.H * EKF_state_tracking.mean + gamma * sigma;
    lower_bound = ekf.H * EKF_state_tracking.mean - gamma * sigma;
    x = transpose(EKF_state_tracking.time);
    x2 = [x, fliplr(x)];
    inBetween = [lower_bound, fliplr(upper_bound)];
    fill(x2, inBetween, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor','m','HandleVisibility', 'off');
    hold on
end

if plot_anomalies
    plot(anomaly_detector.time, anomaly_detector.anomalies, 'o','DisplayName', 'Anomaly detected', 'Color', [18, 174, 0]/255, 'MarkerSize', ms, 'MarkerFaceColor', [18, 174, 0]/255)
    hold on

    false_positives.values = [];
    false_positives.time = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');

    for i = 1:length(patientData.CGM.time)
        t_i = patientData.CGM.time(i);

        [is_in_anomalies, fp_idx] = ismember(t_i, anomaly_detector.time);
        [is_in_faults] = ismember(t_i, faulty_measurement.t);

        if is_in_anomalies && ~is_in_faults
            false_positives.values(end+1) = anomaly_detector.anomalies(fp_idx);
            false_positives.time(end+1) = anomaly_detector.time(fp_idx);
        end
    end
    
    plot(false_positives.time, false_positives.values, 'o','DisplayName', 'False Positive', 'Color', [1 0 0], 'MarkerSize', ms, 'MarkerFaceColor', [1 0 0])
    hold on
end


legend show;
xlim([t_start t_end]);
ylim([-0.5, 350]);
set(gcf, 'Position', get(0, 'Screensize'));
xlabel('Time')
ylabel('Plasma Glucose Concentration [mg/dL]')
title('Glucose Estimator Anomaly Detection')
set(gca, 'FontSize', 14)

%% Function

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
    disp(['Sensitivity: ', num2str(sensitivity)]);
    disp(['Specificity: ', num2str(specificity)]);
end
