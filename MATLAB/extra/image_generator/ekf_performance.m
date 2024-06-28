close all
clear
clc

rng default;

%% Initialization
run('config.m')

show_ceg = false;

if use_tuned_model
    run('param_model.m')
end

horizon = 30;
MARD = 10;
high_uncertainty = 30;

an_value = -5;

Q = eye(13) * high_uncertainty;
R = 100;

ekf_dt = 2; % [min]

model = non_linear_model(tools);
ekf = ekf(model, tools, params, ekf_dt, Q, R);
ekf.dt = ekf_dt;
ekf.CGM_MARD = MARD;

alpha = 0.05;
anomaly_detector = anomaly_detector(alpha);

if use_basal_init_conditions
    [x0, y_minus1] = tools.init_conditions(params);
    u0 = [0,0];
else
    [x0, y_minus1] = tools.set_init_conditions(patientData.CGM.values(1), params);
    u0 = [0,0];
end

EKF_state_tracking.mean = [];
EKF_state_tracking.variance = [];
EKF_state_tracking.time = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');
future_predictions.values = [];
future_predictions.cov = [];
future_predictions.time = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');
residuals = [];

predictions_to_check = [];
CGM_to_check = [];
times_to_check = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');

trajectories = {};

t = t_start;
x = x0;
y = y_minus1;
u = u0;
last_process_update = t_start;
P0 = eye(13) * 12000;
P = P0;

%% Start simulation
% profile on

if simulate_anomalies
    true_CGM.values = [];
    true_CGM.time = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');

    % s = 2/7;
    % e = 5/8;
    s = 2/5;
    e = 2.5/5;
    idxs = floor(length(patientData.CGM.values)*s:length(patientData.CGM.values)*e);
    true_CGM.values = [true_CGM.values, patientData.CGM.values(idxs)'];
    true_CGM.time = [true_CGM.time, patientData.CGM.time(idxs)'];
    patientData.CGM.values(idxs) = an_value;

end

disp('Starting simulation...')
while t <= t_end
    
    [z_k, new_measurement_detected] = sample_measurement(t, patientData);
    [u_k, new_input_detected] = sample_input(t, patientData);
    
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
        end

        x = x_current;
        P = P_current;
        y = y_kminus1;
        u = u_k;
            
        last_process_update = t;

        EKF_state_tracking.mean(:,end+1) = x';
        EKF_state_tracking.variance(end+1) = ekf.H * P * ekf.H';
        EKF_state_tracking.time(end+1, :) = t;

        [temp,temp_cov,~,~] = ekf.predict(x, y, u, P, horizon, params, true);
        future_predictions.values(:,end+1) = temp';
        future_predictions.cov(end+1) = ekf.H * temp_cov * ekf.H';
        future_predictions.time(end+1,:) = t + minutes(horizon);

        if new_measurement_detected && t + minutes(horizon) <= patientData.CGM.time(end)      
            predictions_to_check(end+1) = ekf.H * temp';

            prediction_time = t + minutes(horizon);

            time_differences = abs(patientData.CGM.time - prediction_time);
            [~, nearest_idx] = min(time_differences);
            CGM_to_check(end+1) = patientData.CGM.values(nearest_idx);
            times_to_check(end+1) = prediction_time;
        end

    end
       
    t = next_step(t, patientData);
        
end
% profile off
% profile viewer

if ~show_pred_improvement  
    rmse = mean((predictions_to_check - CGM_to_check).^2)^(1/2)
    T_s = 5;
    tau = delay(CGM_to_check, predictions_to_check);
    horizon
    Delay = abs(tau) * T_s;
    TG = horizon - Delay
end

%% Plot ptred and track

figure

% ms = 2;
ms = 2;
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

lables = {'CGM'};

colors = {'cyan', '#C0C0C0', 'm', [1,0,0], '#77AC30'};


subplot(4, 1, [1 3]);

images(end+1) = plot(patientData.CGM.time, patientData.CGM.values, '-o', 'Color', colors{1}, 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', colors{1});
hold on

if show_confidence_interval
    gamma = 2;
    sigma = sqrt(EKF_state_tracking.variance);
    upper_bound = ekf.H * EKF_state_tracking.mean + gamma * sigma;
    lower_bound = ekf.H * EKF_state_tracking.mean - gamma * sigma;
    x = transpose(EKF_state_tracking.time);
    x2 = [x, fliplr(x)];
    inBetween = [lower_bound, fliplr(upper_bound)];
    fill(x2, inBetween, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor',colors{3},'HandleVisibility', 'off');
    hold on
end

if plot_anomalies
    if simulate_anomalies
        plot(true_CGM.time, true_CGM.values, 'o', 'Color', colors{2}, 'MarkerSize', ms, 'HandleVisibility', 'off', 'MarkerFaceColor', colors{2});
        hold on
    end
    if ~isempty(anomaly_detector.anomalies)
        lables(end+1) = {'Anomalies'};
        images(end+1) = plot(anomaly_detector.time, anomaly_detector.anomalies, 'o', 'Color', colors{4}, 'MarkerSize', ms, 'MarkerFaceColor', colors{4});
        hold on
    end
end

if show_confidence_interval
    lables(end+1) = {'Tracking \pm 2\sigma'};
else
    lables(end+1) = {'Tracking'};
end

images(end+1) = plot(EKF_state_tracking.time, ekf.H * EKF_state_tracking.mean, '-', 'LineWidth', 2*lw, 'Color', colors{3});
hold on

if show_future_predictions
    if show_confidence_interval
        lables(end+1) = {'Prediction \pm 2\sigma'};
    else
        lables(end+1) = {'Prediction'};
    end
    images(end+1) = plot(future_predictions.time, ekf.H * future_predictions.values, '-', 'LineWidth', 2*lw, 'Color', colors{5});
    hold on
    if show_confidence_interval
        gamma = 2;
        sigma = sqrt(future_predictions.cov);
        upper_bound = ekf.H * future_predictions.values + gamma * sigma;
        lower_bound = ekf.H * future_predictions.values - gamma * sigma;
        x = transpose(future_predictions.time);
        x2 = [x, fliplr(x)];
        inBetween = [lower_bound, fliplr(upper_bound)];
        fill(x2, inBetween, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor',colors{5}, 'HandleVisibility', 'off');
        hold on
    end
end

lgd = legend(images, lables);
lgd.Location = 'northoutside';
lgd.Orientation = 'horizontal';
lgd.Box = 'off';

xlim([t_start t_end]);
if ~show_future_predictions
    ylim([an_value-0.5, 300]);
else
    ylim([40, 300]);
end
ylabel('G [mg/dL]')
if do_measurment_update
    if do_chi_sq_test || do_cusum_test
        title('Glucose Estimator Anomaly Detection')
    elseif show_future_predictions
        title('Glucose Estimator Tracking and Prediction')
    else
        title('Glucose Estimator Tracking')
    end
else
    title('Process Model Performance');
end


subplot(4, 1, 4);  % Allocate 1/3 space for the subplot

yyaxis right;
u2 = stem(patientData.IIR.time, patientData.IIR.values, 'filled', 'o', 'Color', [251,142,2]/255,'MarkerSize', ms, 'LineWidth', lw);
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
u1 = stem(patientData.Meal.time, patientData.Meal.values, 'filled', 'square', 'Color', 'g','MarkerSize', ms, 'LineWidth', lw);
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

%% Plot 

if show_ceg
    clarke_for_report(CGM_to_check, predictions_to_check)
end
