close all
clear
clc

state_fields = {'Qsto1','Qsto2','Qgut','Gp','Gt','Gpd','Il','Ip','I1','Id','X','Isc1','Isc2'};
extra_state_fields = {'insulin_to_infuse','last_IIR','CHO_to_eat','D','lastQsto','is_eating'}; 
input_fields = {'CHO', 'IIR'};
true_input_fields = {'CHO_consumed_rate','IIR_dt'};

%% Load data

simulate_anomalies = false;
use_tuned_model = false;
use_true_model = true;
use_known_init_conditions = true;
do_measurment_update = false;
do_chi_sq_test = false;
do_cusum_test = false;
do_plots = true;
use_true_patient = true;

if ~use_true_patient
    filename = "/Users/giovannicherubini/Desktop/Thesis/Code/data/1minsample/adult#001_6.csv";
    tools = utils(filename, state_fields, extra_state_fields, input_fields, true_input_fields);
    patientData.CGM.values = tools.CGMs;
    patientData.CGM.time = tools.Time;
    patientData.Meal.values = tools.CHOs;
    patientData.Meal.time = tools.Time;
    patientData.IIR.values = tools.IIRs;
    patientData.IIR.time = tools.Time;
    basal = tools.IIRs(1);
else
    filename = "none";
    tools = utils(filename, state_fields, extra_state_fields, input_fields, true_input_fields);
    run('database_preprocessor.m')
    basal = 0;
end

%% Initialization

Q = eye(numel(state_fields)) * 15;    % TBD
R = 100;  % TBD

% ekf_dt_values = [0.1, 0.5, 1];
ekf_dt_values = [1];
ekf_dt = ekf_dt_values(1);

if simulate_anomalies
    run('error_gen.m')
end

if use_true_model
    params = patient_01(basal);
else
    params = patient_00(basal);
    if use_tuned_model
        run('param_model.m')
    end
end

model = non_linear_model(tools);
lin_model = linearized_model(tools);
ekf = ekf(model, lin_model, tools, params, ekf_dt, Q, R);

if use_known_init_conditions
    [x0, y_minus1] = tools.init_conditions(params);
    u0.CHO = 0;
    u0.IIR = 0;
else
    [x0, y_minus1] = tools.rand_conditions(params);
    u0.CHO = 0;
    u0.IIR = 0;
end

%% Start simulation

t_start = min([min(patientData.CGM.time), min(patientData.Meal.time), min(patientData.IIR.time)]);
t_end = max([max(patientData.CGM.time), max(patientData.Meal.time), max(patientData.IIR.time)]);

% Loop over the values of ekf_dt
prediction.first = [];
prediction.second = [];
prediction.third = [];
mse = [];
rmse = [];
i = 1;
for ekf_dt = ekf_dt_values
    
    ekf.dt = ekf_dt;
    residuals = [];
    innovation_cov = [];
    EKF_state_tracking.mean = [];
    EKF_state_tracking.variance = [];
    EKF_state_tracking.time = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');
    z_history = [];
    u_history = [];
    v_history = [];
    y_history.CHO_to_eat = [];
    y_history.D = [];
    y_history.last_Q_sto = [];
    y_history.is_eating = [];
    y_history.last_IIR = [];
    y_history.insulin_to_infuse = [];
    
    t = t_start;
    x = x0;
    y = y_minus1;
    u = u0;
    last_process_update = t_start;
    P = eye(numel(state_fields),numel(state_fields))*2200;

    while t <= t_end

        [z_k, new_measurement_detected] = sample_measurement(t, patientData);
        [u_k, new_input_detected] = sample_input(t, patientData);

        if new_input_detected || new_measurement_detected

            dt = convert_to_minutes(t - last_process_update);
            [xp_k, Pp_k, y_kminus1, v_kminus1] = ekf.predict(x, y, u, P, dt, params);
            x_current = xp_k;
            P_current = Pp_k;
                        
            if new_measurement_detected && do_measurment_update
                [xm_k, Pm_k, residual_k, innov_cov_k] = ekf.measurement_update(xp_k,Pp_k,z_k);
                x_current = xm_k;
                P_current = Pm_k;
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
           
        t = t + seconds(1);
    end
    
    % v_history(end+1) = u.CHO;
    % y_history.CHO_to_eat(end+1) = y.CHO_to_eat;
    % y_history.D(end+1) = y.D;
    % y_history.last_Q_sto(end+1) = y.lastQsto;
    % y_history.is_eating(end+1) = y.is_eating;
    % y_history.insulin_to_infuse(end+1) = y.insulin_to_infuse;
    % y_history.last_IIR(end+1) = y.last_IIR;
    % 
    % if length(ekf_dt_values) > 1
    %     if i == 1
    %         prediction.first = EKF_state_tracking;
    %     elseif i == 2
    %         prediction.second = EKF_state_tracking;    
    %     else
    %         prediction.third = EKF_state_tracking;
    %     end
    %     i = i+1;
    % end
    % 
    % mse(end+1) = mean((transpose(tools.BGs)-(EKF_state_tracking.mean(6,1:1/ekf_dt:end)/params.VG)).^2);
    % rmse(end+1) = mean( ((EKF_state_tracking.mean(6,1:1/ekf_dt:end)/params.VG) - (transpose(tools.BGs))).^2 )^(1/2);    % RMSE

    disp("Done");
end

%% Sensor anomaly detection

if do_chi_sq_test
    run("chi_sq_test.m")
elseif do_cusum_test
    run("cusum_test.m")
end

%% Plot

if do_plots
    run('plotting_new.m');
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