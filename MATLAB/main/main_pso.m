close all
clear
clc

rng default;

params_to_estimate = {'kp2','k1','k2','kp1','ki','ke1','kmax','kmin','kabs','kp3','Vmx'};
nvars = length(params_to_estimate);
ub = [0.02,   0.5,    0.5,    5,   0.01,   0.001,  0.1,  0.01,   0.3,  0.1,    0.1] * 2;
%     kp2,    k1,     k2,     kp1, ki,     ke1,    kmax, kmin,   kabs, kp3,    Vmx  
lb = [0.0001, 0.0001, 0.0001, 2,   0.0040, 0.0001, 0.01, 0.0001, 0.01, 0.01, 0.001] / 2;
% 
% params_to_estimate = {'VG','m1','CL','Vmx','k1','Km0','k2','kp2','kmax','kmin','kabs','ki','u2ss'};
% nvars = length(params_to_estimate);
% ub = [2,   0.4,   1.5,    0.1,   0.1,   300, 0.2, 0.01, 0.1,  0.01,   0.3, 0.01, 10] * 1;
% lb = [1.5, 0.1,   0.5,    0.01,  0.01,  200, 0.05, 0.0001, 0.001, 0.0001, 0.05, 0.0040, 0.1] / 1;

% params_to_estimate = {'VG','m1','CL','Vmx','k1','Km0','k2','kp2','kmax','kmin','kabs'};
% nvars = length(params_to_estimate);
% ub = [2,   0.4,   1.5,    0.1,   0.1,   300, 0.2, 0.01,     0.1,  0.01,   0.3] * 3;
% lb = [1.5, 0.1,   0.5,    0.01,  0.01,  200, 0.05, 0.0001, 0.001, 0.0001, 0.05] / 3;

% params_to_estimate = {'VG','m1','CL','Vmx','k1','Km0','k2','kp2', 'kp1', 'ki'};
% nvars = length(params_to_estimate);
% ub = [2,   0.4,   1.5,    0.1,   0.1,   300, 0.2, 0.01,  5, 0.01  ] * 1;
% lb = [1.5, 0.1,   0.5,    0.01,  0.01,  200, 0.05, 0.0001, 2,  0.0040] / 1;


run('config.m')
use_CGM_to_tune = true;

model = non_linear_model(tools);
window_size = experiment_total_time;

%%
high_uncertainty = 1000;
Q = diag([10,10,10,high_uncertainty,high_uncertainty,high_uncertainty,0,0,0,0,0,0,0]);    % TBD
% Q = eye(length(state_fields)) * high_uncertainty;
R = 1000;  % TBD
ekf_dt = 1; % [min]

lin_model = linearized_model(tools);
ekf = ekf(model, lin_model, tools, params, ekf_dt, Q, R);
ekf.dt = ekf_dt;

[x0_, ymin1_] = tools.init_conditions(params);
u0.CHO = 0;
u0.IIR = 0;
x0 = x0_;
ymin1 = ymin1_;

t = t_start;
% profile on;
while t < t_end 
    disp("-------")
    window = set_window(window_size, t, x0, ymin1, u0, t_end);
    
    objective_pso = @(p) objective(p, params, ekf, patientData, window, params_to_estimate, use_CGM_to_tune);
    
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
            [xk, Pk, ykmin1, ~] = ekf.predict(x, y, u, P, dt, patient, false);
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
