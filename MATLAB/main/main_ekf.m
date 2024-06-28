close all
clear
clc

rng default;

%% Initialization
run('config.m')

if use_tuned_model
    run('param_model.m')
end


horizon = 30;

MARD = 8;
high_uncertainty = 30;


Q = eye(length(state_fields)) * high_uncertainty;
R = 100;

ekf_dt = 1; % [min]
% ekf_dt = 2.5;

bounds = [100, 150] * params.VG;

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
P0 = eye(length(state_fields)) * 12000;
% P0 = Q * 6000;
P = P0;

%% Start simulation
% profile on

done = false;

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
    patientData.CGM.values(idxs) = 0;

    % spikes_fractions = [15/16, 2/3, 1/6, 9/16, 11/16, 1/2];
    % spikes_values = [180, 200, 30, 170, 350, 250];
    spikes_fractions = [];
    spikes_values = [];

    for i = 1:length(spikes_values)
        idx = floor(length(patientData.CGM.values)*spikes_fractions(i));
        value = spikes_values(i);
        true_CGM.values(end+1) = patientData.CGM.values(idx);
        true_CGM.time(end+1) = patientData.CGM.time(idx);
        patientData.CGM.values(idx) = value;
    end

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

        delay_t = 0;

        if show_pred_improvement
            [trajectory] = ekf.predict_save_trajectories(x, y, u, P, horizon, params, true, t);
            trajectories{end+1} = trajectory;
        else
            if save_energy && (x(6) >= bounds(1) && x(6) <= bounds(2))
                predict = false;
            else
                predict = true;
            end

            if predict
                [temp,temp_cov,~,~] = ekf.predict(x, y, u, P, horizon, params, true);
                future_predictions.values(:,end+1) = temp';
                future_predictions.cov(end+1) = ekf.H * temp_cov * ekf.H';
                % future_predictions.time(end+1,:) = t + minutes(horizon);
                
                future_predictions.time(end+1,:) = t + minutes(horizon) - minutes(delay_t);
            end
        end
        
        if new_measurement_detected && t + minutes(horizon) <= patientData.CGM.time(end) && ~show_pred_improvement && predict        
            predictions_to_check(end+1) = ekf.H * temp';

            prediction_time = t + minutes(horizon);
            prediction_time = t + minutes(horizon) - minutes(delay_t);

            time_differences = abs(patientData.CGM.time - prediction_time);
            [~, nearest_idx] = min(time_differences);
            % chosen_time = patientData.CGM.time(nearest_idx)
            CGM_to_check(end+1) = patientData.CGM.values(nearest_idx);
            times_to_check(end+1) = prediction_time;
        end

        if loose_track_of_time
            m = floor(length(patientData.CGM.values)/2);
            if t > patientData.CGM.time(m) && ~done
                done = true;
                P = P0;
                x = x0;
            end
        end

    end
       
    t = next_step(t, patientData);
        
end
% profile off
% profile viewer

if ~show_pred_improvement  
    rmse = mean((predictions_to_check - CGM_to_check).^2)^(1/2)
    % T_s = convert_to_minutes(mean(diff(times_to_check)));
    T_s = 5;
    tau = delay(CGM_to_check, predictions_to_check);
    horizon
    Delay = abs(tau) * T_s;
    TG = horizon - Delay

end

%% Plot

if do_plots
    run('plotting.m');
end

disp("Done");
