close all
clear
clc

%% Config
state_fields = {'Qsto1','Qsto2','Qgut','Gp','Gt','Gpd','Il','Ip','I1','Id','X','Isc1','Isc2'};
extra_state_fields = {'insulin_to_infuse','last_IIR','CHO_to_eat','D','lastQsto','is_eating'}; 
input_fields = {'CHO', 'IIR'};
true_input_fields = {'CHO_consumed_rate','IIR_dt'};

use_true_patient = true;
use_tuned_model = true;
use_true_model = false;

use_known_init_conditions = true;
do_measurment_update = false;
compute_mse = true;

simulate_anomalies = false;
do_chi_sq_test = true;
do_cusum_test = false;

do_plots = true;
if do_plots
    plot_true_database = false;
    all_states = true;
    only_Gpd = true;
    plot_anomalies = true;
    plot_complete_history = false;
    % show_confidence_interval = false;
end

%% Load data
disp('Loading dataset...')
if ~use_true_patient
    filename = "/Users/giovannicherubini/Desktop/Thesis/Code/data/1minsample/adult#001_6.csv";
    tools = utils(filename, state_fields, extra_state_fields, input_fields, true_input_fields);
    patientData.CGM.values = tools.CGMs;
    patientData.CGM.time = tools.Time;
    patientData.Meal.values = tools.CHOs;
    patientData.Meal.time = tools.Time;
    patientData.IIR.values = tools.IIRs;
    patientData.IIR.time = tools.Time;
    patientData.BG.values = tools.BGs;
    patientData.BG.time = tools.Time;
    basal = tools.IIRs(1);
    if use_true_model
        params = patient_01(basal);
    else
        params = patient_00(basal);
        if use_tuned_model
            run('param_model.m')
        end
    end
else
    filename = "none";
    tools = utils(filename, state_fields, extra_state_fields, input_fields, true_input_fields);
    run('database_preprocessor.m')
    basal = 0;
    dailyBasal = 18;
    params = patient_11(dailyBasal);
    if use_tuned_model
        run('param_model.m')
    end
end

disp('Dataset loaded')

%% Initialization
Q = eye(numel(state_fields)) * 15;    % TBD
R = 20;  % TBD
ekf_dt = 1; % [min]

% if simulate_anomalies
%     run('error_gen.m')
% end

model = non_linear_model(tools);
lin_model = linearized_model(tools);
ekf = ekf(model, lin_model, tools, params, ekf_dt, Q, R);
ekf.dt = ekf_dt;

alpha = 0.01;
anomaly_detector = anomaly_detector(alpha);

if use_known_init_conditions
    [x0, y_minus1] = tools.init_conditions(params);
    u0.CHO = 0;
    u0.IIR = 0;
else
    [x0, y_minus1] = tools.rand_conditions(params);
    u0.CHO = 0;
    u0.IIR = 0;
end

t_start = min([min(patientData.CGM.time), min(patientData.Meal.time), min(patientData.IIR.time)]);
t_end = max([max(patientData.CGM.time), max(patientData.Meal.time), max(patientData.IIR.time)]);
EKF_state_tracking.mean = [];
EKF_state_tracking.variance = [];
EKF_state_tracking.time = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');

t = t_start;
x = x0;
y = y_minus1;
u = u0;
last_process_update = t_start;
P = eye(numel(state_fields),numel(state_fields))*2200;

%% Start simulation
% profile on
disp('Starting simulation...')
while t <= t_end

    [z_k, new_measurement_detected] = sample_measurement(t, patientData);
    [u_k, new_input_detected] = sample_input(t, patientData);

    if new_input_detected || new_measurement_detected

        dt = convert_to_minutes(t - last_process_update);
        if plot_complete_history
            [xp_k, Pp_k, y_kminus1, v_kminus1, ekf] = ekf.predict_and_save(x, y, u, P, dt, params, last_process_update);
        else
            [xp_k, Pp_k, y_kminus1, v_kminus1] = ekf.predict(x, y, u, P, dt, params);
        end
        x_current = xp_k;
        P_current = Pp_k;

        flag = true;
        if t == t_start || flag == true
            t_history = t;
            flag = false;
        else
            t_history(end+1,:) = t;
        end
                    
        if new_measurement_detected && do_measurment_update
            [xm_k, Pm_k, residual_k, innov_cov_k] = ekf.measurement_update(xp_k,Pp_k,z_k);
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

        EKF_state_tracking.mean(:,end+1) = tools.convert_to_vector(x);   % (time k)
        % EKF_state_tracking.variance(:,:,end+1) = P;
        EKF_state_tracking.variance(end+1) = P(6,6);
        EKF_state_tracking.time(end+1, :) = t;

    end
       
    t = next_step(t, patientData);
        
end
% profile off
% profile viewer

if compute_mse
    if ~use_true_patient
        mse = mean((transpose(patientData.BG.values) - EKF_state_tracking.mean(6, :)/params.VG).^2)
        rmse = mse^(1/2)
    else
        mse = 0;
        for i = 1:length(patientData.CGM.values)
            t = patientData.CGM.time(i);
            idx = find(EKF_state_tracking.time == t);
            mse = mse + (patientData.CGM.values(i) - EKF_state_tracking.mean(6,idx)/params.VG)^2;
        end
        mse = mse / i
        rmse = mse^(1/2)
    end
end

disp("Done");

%% Plot

if do_plots
    run('plotting.m');
end

%% Extra functions

function [z_k, new_measurement_detected] = sample_measurement(t, patientData)
    z_k = [];
    [new_measurement_detected, idx] = ismember(t, patientData.CGM.time);
    if new_measurement_detected
        z_k = patientData.CGM.values(idx);
    end
end

function [u_k, new_input_detected] = sample_input(t, patientData)
    u_k.CHO = 0;
    u_k.IIR = 0;
    [new_meal_detected, idx1] = ismember(t, patientData.Meal.time);
    [new_IIR_detected, idx2] = ismember(t, patientData.IIR.time);
    if new_meal_detected
        u_k.CHO = patientData.Meal.values(idx1);
    end
    if new_IIR_detected
        u_k.IIR = patientData.IIR.values(idx2);
    end
    if new_meal_detected || new_IIR_detected
        new_input_detected = true;
    else
        new_input_detected = false;
    end
end

function dt = convert_to_minutes(duration)
    [hours, minutes, seconds] = hms(duration);
    dt = hours * 60 + minutes + seconds / 60;
end

function next_t = next_step(t, patientData)
    idx_CGM = find(patientData.CGM.time > t);
    idx_IIR = find(patientData.IIR.time > t);
    idx_Meal = find(patientData.Meal.time > t);

    next_t = min([min(patientData.CGM.time(idx_CGM)), min(patientData.Meal.time(idx_Meal)), min(patientData.IIR.time(idx_IIR))]);
    
    if isempty(next_t)
        next_t = t + seconds(1);
    end

end