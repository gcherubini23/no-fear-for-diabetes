close all
clear
clc

rng default;

params_to_estimate = {'kp2','k1','k2','kp1','ki','ke1','kmax','kmin','kabs','kp3','Vmx'};
nvars = length(params_to_estimate);

prior = [0.0021, 0.065, 0.079, 2.70, 0.0079, 0.0005, 0.0558, 0.008, 0.0570, 0.009, 0.047];

span = 6;
ub = prior*span;
lb = prior/span;

run('config.m')

window_size = experiment_total_time;

%%
high_uncertainty = 1;
Q = eye(13) * high_uncertainty;
R = 20;  % TBD
ekf_dt = 1; % [min]

model = non_linear_model(tools);
ekf = ekf(model, tools, params, ekf_dt, Q, R);
ekf.dt = ekf_dt;
ekf.CGM_MARD = 8;

[x0_, ymin1_] = tools.init_conditions(params);
% [x0_, ymin1_] = tools.set_init_conditions(patientData.CGM.values(1), params);

u0.CHO = 0;
u0.IIR = 0;
x0 = x0_;
ymin1 = ymin1_;

t = t_start;
% profile on;
while t < t_end 
    disp("-------")
    window = set_window(window_size, t, x0, ymin1, u0, t_end);
    
    objective_pso = @(p) objective(p, params, ekf, patientData, window, params_to_estimate);
    
    % options.MaxTime = 3600;
    if t > t_start
        options.InitialPoints = points.X;
    end
    options.Display = 'iter';
    % options.MaxIterations = 50;
    options.FunctionTolerance = 0.0001;
    [final_p,fval,~,~,points] = particleswarm(objective_pso, nvars, lb, ub, options);

    % [final_p,fval] = particleswarm(objective_pso, nvars, lb, ub, options);
    
    t = t + window.size;
end
% profile off;
% profile viewer;

disp('Done');


%% Functions

function f = objective(p, patient, ekf, patientData, window, params_to_estimate)
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
                predictions(end+1) = ekf.H * x';
            end
            last_process_update = t;
        end
        
        % t = next_step(t, patientData);
        t = t + minutes(1);
    end
    
    gt = measurements;

    % f = mean((predictions/patient.VG - gt).^2)^(1/2);
    f = mean(abs(log(gt)-log(predictions)));
    % f = mean(abs((predictions/patient.VG - gt)));

end
