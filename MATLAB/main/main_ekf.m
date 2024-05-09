close all
clear
clc

rng default;

%% Initialization
run('config.m')

if use_tuned_model
    run('param_model.m')
end

% 1000 1000 is good for true patient
% 25 160 is good for simulated patient

high_uncertainty = 1000;
Q = diag([10,10,10,high_uncertainty,high_uncertainty,high_uncertainty,0,0,0,0,0,0,0]);    % TBD
% Q = eye(length(state_fields)) * high_uncertainty;
R = 1000;  % TBD
ekf_dt = 2; % [min]

% if simulate_anomalies
%     run('error_gen.m')
% end

model = non_linear_model(tools);
lin_model = linearized_model(tools);
ekf = ekf(model, lin_model, tools, params, ekf_dt, Q, R);
ekf.dt = ekf_dt;

alpha = 0.05;
anomaly_detector = anomaly_detector(alpha);

if use_basal_init_conditions
    [x0, y_minus1] = tools.init_conditions(params);
    u0.CHO = 0;
    u0.IIR = 0;
else
    [x0, y_minus1] = tools.set_init_conditions(patientData.CGM.values(1), params);
    u0.CHO = 0;
    u0.IIR = 0;
end

EKF_state_tracking.mean = [];
EKF_state_tracking.variance = [];
EKF_state_tracking.time = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');
future_predictions.values = [];
future_predictions.time = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');
residuals = [];

t = t_start;
x = x0;
y = y_minus1;
u = u0;
last_process_update = t_start;
P = Q * 3;

%% Start simulation
% profile on

flag = true;
disp('Starting simulation...')
horizon = 45;
while t <= t_end

    [z_k, new_measurement_detected] = sample_measurement(t, patientData);
    [u_k, new_input_detected] = sample_input(t, patientData);

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
            % ekf = ekf.update_sensor_cov(z_k);
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

        if new_measurement_detected
            EKF_state_tracking.mean(:,end+1) = tools.convert_to_vector(x);   % (time k)
            % EKF_state_tracking.variance(:,:,end+1) = P;
            EKF_state_tracking.variance(end+1) = P(6,6);
            EKF_state_tracking.time(end+1, :) = t;

            [temp,~,~,~] = ekf.predict(x, y, u, P, horizon, params, false);
            future_predictions.values(end+1) = temp.Gpd;
            future_predictions.time(end+1,:) = t + minutes(horizon);
        end

    end
       
    t = next_step(t, patientData);
        
end
% profile off
% profile viewer

if compute_mse
    if ~use_true_patient
        mse = mean((transpose(patientData.BG.values) - EKF_state_tracking.mean(6, :)/params.VG).^2);
        rmse = mse^(1/2)
    else
        mse = 0;
        for i = 1:length(patientData.CGM.values)
            t = patientData.CGM.time(i);
            idx = find(EKF_state_tracking.time == t);
            mse = mse + (patientData.CGM.values(i) - EKF_state_tracking.mean(6,idx)/params.VG)^2;
        end
        mse = mse / length(patientData.CGM.values);
        rmse = mse^(1/2)
    end
end

%% Plot

if do_plots
    run('plotting.m');
    [total, percentage] = clarke(patientData.CGM.values',EKF_state_tracking.mean(6,:)/params.VG)
end

disp("Done");
