close all
clear
clc

state_fields = {'Qsto1','Qsto2','Qgut','Gp','Gt','Gpd','Il','Ip','I1','Id','X','Isc1','Isc2'};
extra_state_fields = {'insulin_to_infuse','last_IIR','CHO_to_eat','D','lastQsto','is_eating'}; 
input_fields = {'CHO', 'IIR'};
true_input_fields = {'CHO_consumed_rate','IIR_dt'};
% 
% params_to_estimate = {'kp2','k1','k2','kp1','ki','ke1','kmax','kmin','kabs','kp3','Vmx'};
% nvars = length(params_to_estimate);
% ub = [0.02,   0.5,    0.5,    5,   0.01,   0.001,  0.1,  0.01,   0.3,  0.1,    0.1];
% %     kp2,    k1,     k2,     kp1, ki,     ke1,    kmax, kmin,   kabs, kp3,    Vmx  
% lb = [0.0001, 0.0001, 0.0001, 2,   0.0040, 0.0001, 0.01, 0.0001, 0.01, 0.01, 0.001];
% 
params_to_estimate = {'VG','m1','CL','Vmx','k1','Km0','k2','kp2','kmax','kmin','kabs','ki','u2ss'};
nvars = length(params_to_estimate);
ub = [2,   0.4,   1.5,    0.1,   0.1,   300, 0.2, 0.01, 0.1,  0.01,   0.3, 0.01, 10] * 1;
lb = [1.5, 0.1,   0.5,    0.01,  0.01,  200, 0.05, 0.0001, 0.001, 0.0001, 0.05, 0.0040, 0.1] / 1;


use_true_patient = true;
use_anderson = true;
use_tmoore = false;
use_CGM_to_tune = true;


disp('Loading dataset...')
if ~use_true_patient
    filename = "/Users/giovannicherubini/Desktop/Thesis/Code/data/1minsample/adult#001_5.csv";
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
    patient = patient_00(basal);
else
    if use_anderson
        patient_ID = 11;
        date = '11-Feb-2013 06:30:00';
        % date = '26-Jan-2013 06:30:00';
        % date = '28-Jan-2013 06:30:00';

        days_to_examine = 2;
        % days_to_examine = 30;
        % days_to_examine = 'all';
    end

    if use_tmoore
        patient_ID = -1;
        date = '20-Jan-2022 00:00:00';

        days_to_examine = 2;
        % days_to_examine = 30;
        % days_to_examine = 'all';

    end

    start_day = datetime(date,'InputFormat', 'dd-MMM-yyyy HH:mm:ss');

    filename = "none";
    plot_true_database = false;
    tools = utils(filename, state_fields, extra_state_fields, input_fields, true_input_fields);
    run('database_preprocessor.m')
    % basal = 0;
    dailyBasal = 10;
    patient = patient_11(dailyBasal);
end
disp('Dataset loaded...')

t_start = min([min(patientData.CGM.time), min(patientData.Meal.time), min(patientData.IIR.time)]);
t_end = max([max(patientData.CGM.time), max(patientData.Meal.time), max(patientData.IIR.time)]);
experiment_total_time = t_end - t_start;

model = non_linear_model(tools);
window_size = experiment_total_time;

%%
Q = diag([0,0,0,100,100,100,0,0,0,0,0,0,0]);    % TBD
R = 10;  % TBD
ekf_dt = 1; % [min]

lin_model = linearized_model(tools);
ekf = ekf(model, lin_model, tools, patient, ekf_dt, Q, R);
ekf.dt = ekf_dt;

[x0_, ymin1_] = tools.init_conditions(patient);
u0.CHO = 0;
u0.IIR = 0;
x0 = x0_;
ymin1 = ymin1_;

t = t_start;
% profile on;
while t < t_end 
    disp("-------")
    window = set_window(window_size, t, x0, ymin1, u0, t_end);
    
    objective_pso = @(p) objective(p, patient, ekf, patientData, window, params_to_estimate, use_CGM_to_tune);
    
    % options.MaxTime = 3600;
    if t > t_start
        options.InitialPoints = points.X;
    end
    options.Display = 'iter';
    options.MaxIterations = 200;
    options.FunctionTolerance = 0.0001;
    [final_p,fval,~,~,points] = particleswarm(objective_pso, nvars, lb, ub, options);

    % [final_p,fval] = particleswarm(objective_pso, nvars, lb, ub, options);
    
    % To implement recursive approach
    
    t = t + window.size;
end
% profile off;
% profile viewer;

disp('Done');


%% Functions

function f = objective(p, patient, ekf, patientData, window, params_to_estimate, use_CGM_to_tune)
    patient = patient.set_params(params_to_estimate,p);
    t = window.t_start;
    predictions = [];
    % x = window.x0;
    % y = window.ymin1;
    [x, y] = ekf.tools.init_conditions(patient);
    u = window.u0;

    last_process_update = window.t_start;
    P = ekf.Q * 3;
    measurements = [];
    while t <= window.t_end
        [uk, new_input_detected] = sample_input(t, patientData);
        [zk, new_measurement_detected] = sample_measurement(t, patientData);
        if new_input_detected || new_measurement_detected
            dt = convert_to_minutes(t - last_process_update);
            [xk, Pk, ykmin1, ~] = ekf.predict(x, y, u, P, dt, patient);
            x = xk;
            P = Pk;
            y = ykmin1;
            u = uk;
            if new_measurement_detected
                measurements(end+1) = zk;
                predictions(end+1) = x.Gpd;
            end
            last_process_update = t;
        end
        t = next_step(t, patientData);
    end
    
    if use_CGM_to_tune
        gt = measurements;
    else
        idx = find(patientData.BG.time >= window.t_start & ...
                    patientData.BG.time <= window.t_end);
        gt = transpose(patientData.BG.values(idx));
    end
    % f = mean((predictions/patient.VG - gt).^2);
    f = mean(abs(log(gt)-log(predictions/patient.VG)));
    % f = mean(abs((predictions/patient.VG - gt)));

end

function window = set_window(window_size, t, x0, ymin1, u0, t_end)   
    window.t_start = t;
    window.t_end = min(t + window_size, t_end);
    window.size = window.t_end - window.t_start;
    window.x0 = x0;
    window.u0 = u0;
    window.ymin1 = ymin1;
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

function [z_k, new_measurement_detected] = sample_measurement(t, patientData)
    z_k = [];
    [new_measurement_detected, idx] = ismember(t, patientData.CGM.time);
    if new_measurement_detected
        z_k = patientData.CGM.values(idx);
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