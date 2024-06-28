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


Q = eye(13) * uncertainty;
R = 100;

ekf_dt = 2;

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
P0 = eye(13) * 12000;
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
figure

ms = 3;
lw = 0.5;

images = [];

font_size = 10;
figureWidth = 5.7;
figureHeight = 4;
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, figureWidth, figureHeight]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figureWidth, figureHeight]);
set(gcf, 'PaperPositionMode', 'auto');

lables = {'CGM', 'Tracking \pm 2\sigma'};
colors = {'m', 'cyan', [1 0 0], [18, 174, 0]/255, 'k'};

subplot(4, 1, [1 3]);

images(end+1) = plot(used_measurement.t, used_measurement.z, '-o', 'Color', colors{2}, 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', colors{2});
hold on

images(end+1) = plot(EKF_state_tracking.time, ekf.H * EKF_state_tracking.mean, '-', 'LineWidth', 2*lw, 'Color', colors{1});
hold on

if ~isempty(faulty_measurement.t)
    lables(end+1) = {'Anomaly not detected'};
    images(end+1) = plot(faulty_measurement.t, faulty_measurement.z, 'o', 'Color', colors{3}, 'MarkerSize', ms, 'MarkerFaceColor',colors{3});
    hold on
end

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
    if ~isempty(anomaly_detector.time)
        lables(end+1) = {'Anomaly detected'};
        images(end+1) = plot(anomaly_detector.time, anomaly_detector.anomalies, 'o', 'Color', colors{4}, 'MarkerSize', ms, 'MarkerFaceColor', colors{4});
        hold on
    end

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
    
    if ~isempty(false_positives.time)
        lables(end+1) = {'False Positive'};
        images(end+1) = plot(false_positives.time, false_positives.values, 'o','DisplayName', 'False Positive', 'Color', colors{5}, 'MarkerSize', ms, 'MarkerFaceColor', colors{5});
        hold on
    end
end


lgd = legend(images, lables);
lgd.Location = 'northoutside';
lgd.Orientation = 'horizontal';
lgd.Box = 'off';
lgd.NumColumns = 2;

xlim([t_start t_end]);
ylim([40, 300]);
ylabel('G [mg/dL]')

title('Glucose Estimator Anomaly Detection')


subplot(4, 1, 4);  % Allocate 1/3 space for the subplot

yyaxis right;
u2 = stem(patientData.IIR.time, patientData.IIR.values, 'filled', 'o', 'Color', [251,142,2]/255,'MarkerSize', 2, 'LineWidth', lw);
ylabel('Insulin [IU]');
if use_true_patient
    ylim([0,7])
else
    ylim([0,15])
end
xlim([t_start, t_end])
hold on
ax = gca;
ax.YColor = 'k'; % Set color of left y-axis to black

yyaxis left;
u1 = stem(patientData.Meal.time, patientData.Meal.values, 'filled', 'square', 'Color', 'g','MarkerSize', 2, 'LineWidth', lw);
ylabel('CHO [g]');
if use_true_patient
    ylim([0,60])
else
    ylim([0,90])
end
xlim([t_start, t_end])
ax = gca;
ax.YColor = 'k';

xlabel('Time');
lgd2 = legend([u1,u2], {'Meal', 'Insulin infused'});
lgd2.Location = 'northoutside';
lgd2.Orientation = 'horizontal';
lgd2.Box = 'off';

set(findall(gcf, '-property', 'FontSize'), 'FontSize', font_size);


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
