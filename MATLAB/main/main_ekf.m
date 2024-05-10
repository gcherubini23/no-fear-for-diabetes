close all
clear
clc

rng default;

%% Initialization
run('config.m')

if use_tuned_model
    run('param_model.m')
end


% 25 160 is good for simulated patient

high_uncertainty = 30;
% Q = diag([0,0,0,high_uncertainty,high_uncertainty,high_uncertainty,0,0,0,0,0,0,0]);    % TBD
Q = eye(length(state_fields)) * high_uncertainty;
R = 100;  % TBD
ekf_dt = 1; % [min]

model = non_linear_model(tools);
ekf = ekf(model, tools, params, ekf_dt, Q, R);
ekf.dt = ekf_dt;

alpha = 0.1;
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

t = t_start;
x = x0;
y = y_minus1;
u = u0;
last_process_update = t_start;
P = Q * 600;

%% Start simulation
% profile on

if simulate_anomalies
    patientData.CGM.values(522) = 180;
    patientData.CGM.values(floor(length(patientData.CGM.values)*2/3)) = 200;
    patientData.CGM.values(floor(length(patientData.CGM.values)*1/6)) = 30;
    patientData.CGM.values(floor(length(patientData.CGM.values)*2/6)) = 170;
    patientData.CGM.values(floor(length(patientData.CGM.values)*11/16)) = 350;
    patientData.CGM.values(floor(length(patientData.CGM.values)*2/5:floor(length(patientData.CGM.values)*2.5/5))) = 0;
end

flag = true;
disp('Starting simulation...')
horizon = 45;
i = 0;
while t <= t_end
    
    [z_k, new_measurement_detected] = sample_measurement(t, patientData);
    [u_k, new_input_detected] = sample_input(t, patientData);
    
    % if i < 100 && new_measurement_detected
    %     new_measurement_detected = false;
    % end
    % i = i+1;

    if new_input_detected || new_measurement_detected

        dt = convert_to_minutes(t - last_process_update);
        if plot_complete_history
            [xp_k, Pp_k, y_kminus1, v_kminus1, ekf] = ekf.predict_and_save(x, y, u, P, dt, params, last_process_update, true);
        else
            [xp_k, Pp_k, y_kminus1, v_kminus1] = ekf.predict(x, y, u, P, dt, params, true);
        end
        x_current = xp_k;
        P_current = Pp_k;

        if t == t_start || flag == true
            t_history = t;
            flag = false;
        else
            t_history(end+1,:) = t;
        end
                    
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
            end
        end

        x = x_current;
        P = P_current;
        y = y_kminus1;
        u = u_k;
            
        last_process_update = t;

        EKF_state_tracking.mean(:,end+1) = x';
        EKF_state_tracking.variance(end+1) = P(6,6);
        EKF_state_tracking.time(end+1, :) = t;
        
        [temp,temp_cov,~,~] = ekf.predict(x, y, u, P, horizon, params, true);
        future_predictions.values(:,end+1) = temp';
        future_predictions.cov(end+1) = temp_cov(6,6);
        future_predictions.time(end+1,:) = t + minutes(horizon);
        
        if new_measurement_detected && ~anomaly
            predictions_to_check(end+1) = x(6)/params.VG;
            CGM_to_check(end+1) = z_k;
        end
    end
       
    t = next_step(t, patientData);
        
end
% profile off
% profile viewer

if compute_mse
    % if ~use_true_patient
    %     mse = mean((transpose(patientData.BG.values) - EKF_state_tracking.mean(6, :)/params.VG).^2);
    %     rmse = mse^(1/2)
    % else
    %     mse = 0;
    %     for i = 1:length(patientData.CGM.values)
    %         t = patientData.CGM.time(i);
    %         idx = find(EKF_state_tracking.time == t);
    %         mse = mse + (patientData.CGM.values(i) - EKF_state_tracking.mean(6,idx)/params.VG)^2;
    %     end
    %     mse = mse / length(patientData.CGM.values);
    %     rmse = mse^(1/2)
    % end

    rmse = mean((predictions_to_check - CGM_to_check).^2)^(1/2)
end

%% Plot

if do_plots
    run('plotting.m');
    % [total, percentage] = clarke(patientData.CGM.values',EKF_state_tracking.mean(6,:)/params.VG)
end

disp("Done");
