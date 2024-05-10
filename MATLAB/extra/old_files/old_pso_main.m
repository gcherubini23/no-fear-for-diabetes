close all
clear
clc

state_fields = {'Qsto1','Qsto2','Qgut','Gp','Gt','Gpd','Il','Ip','I1','Id','X','Isc1','Isc2'};
extra_state_fields = {'insulin_to_infuse','last_IIR','CHO_to_eat','D','lastQsto','is_eating'}; 
input_fields = {'CHO', 'IIR'};
true_input_fields = {'CHO_consumed_rate','IIR_dt'};
filename = "/Users/giovannicherubini/Desktop/Thesis/Code/data/1minsample/adult#001_5.csv";

params_to_estimate = {'kp2','k1','k2','kp1','ki','ke1','kmax','kmin','kabs','kp3','Vmx','Gb'};
nvars = 12;
ub = [0.5,0.5,0.5,6,0.01,0.001,0.1,0.01,0.1,0.1,0.1,160];
lb = [0.0001,0.0001,0.0001,1,0.0001,0.0001,0.001,0.0001,0.001,0.0001,0.001,50];

% params_to_estimate = {'VG','m1','CL','Vmx','k1','Km0','k2','kp2','kmax','kmin','kabs','ki'};
% nvars = length(params_to_estimate);
% ub = [2,   0.4,   1.5,    0.1,   0.1,   300, 0.2, 0.01, 0.1,  0.01,   0.3, 0.01];
% lb = [1.5, 0.1,   0.5,    0.01,  0.01,  200, 0.05, 0.0001, 0.001, 0.0001, 0.05, 0.0040];

tools = utils(filename);
model = non_linear_model(tools);
basal = tools.IIRs(1);
patient = patient_00(basal);
true_patient = patient_01(basal);

cgm_dt = 1;
experiment_total_time = length(tools.CGMs) * cgm_dt;
% window_size = 25;
window_size = experiment_total_time;
t = 0;

[x0_, ymin1_] = tools.init_conditions(patient);
x0 = x0_;
ymin1 = ymin1_;

while ~stop_simulation(tools, cgm_dt, t)
    disp("-------")
    window = set_window(tools, window_size, t, cgm_dt, x0, ymin1);
    
    objective_pso = @(p) objective_2(p, patient, model, window, tools, params_to_estimate);
    
    % options.ObjectiveLimit = 3.5;
    % options.MaxTime = 1800;
    if t > 0
        options.InitialPoints = points.X;
    end
    options.Display = 'iter';
    options.MaxIterations = 100;
    options.FunctionTolerance = 0.01;
    [final_p,fval,~,~,points] = particleswarm(objective_pso, nvars, lb, ub, options);

    % patient = patient.set_params(params_to_estimate, final_p);
    % k = 1;
    % x = x0_;
    % y = ymin1_;
    % while k < window.t_start
    %     u_k.CHO = tools.CHOs(k);
    %     u_k.IIR = tools.IIRs(k);
    % 
    %     if k > window.t_start
    %         [x_k, y_kminus1, ~] = tools.euler_solve(model,true_patient,x,y,u,window.dt);
    %     else
    %         x_k = x;
    %         y_kminus1 = y;
    %     end
    % 
    %     x = x_k;
    %     y = y_kminus1;
    %     u = u_k;
    % 
    %     k = k + window.dt;
    % end
    % 
    % x0 = x;
    % ymin1 = y;

    t = t + window_size
    fval = fval
    final_p = final_p

end

disp('Done');



%% Functions

function f = objective(p, patient, model, window, tools, params_to_estimate)
    patient = patient.set_params(params_to_estimate,p);
    t = window.t_start;
    dt = window.dt;
    predictions = [];
    x = window.x0;
    y = window.ymin1;
    while t <= window.t_end
        u_k.CHO = tools.CHOs(t);
        u_k.IIR = tools.IIRs(t);

        if t > window.t_start
            [x_k, y_kminus1, ~] = tools.euler_solve(model,patient,x,y,u,dt);
        else
            x_k = x;
            y_kminus1 = y;
        end

        x = x_k;
        y = y_kminus1;
        u = u_k;

        predictions(end+1) = x.Gpd;
        
        t = t + dt;
    end

    gt = transpose(tools.BGs(window.t_start:window.t_end));
    f = mean((predictions/patient.VG - gt).^2);

end

function f = objective_2(p, patient, model, window, tools, params_to_estimate)
    cgm_dt = 1;
    patient = patient.set_params(params_to_estimate,p);
    t = window.t_start;
    dt = window.dt;
    predictions = [];
    % x = window.x0;
    % y = window.ymin1;
    [x, y] = tools.init_conditions(patient);
    while t <= window.t_end
        u_k = [tools.CHOs(t), tools.IIRs(t)];

        if t > window.t_start
            [x_k, y_kminus1, ~] = tools.euler_solve(model,patient,x,y,u,dt);
        else
            x_k = x;
            y_kminus1 = y;
        end

        x = x_k;
        y = y_kminus1;
        u = u_k;

        if mod(t, cgm_dt) == 0
            predictions(end+1) = x(6);
        end
        
        t = t + dt;
    end

    gt = transpose(tools.CGMs);
    % f = mean((predictions/patient.VG - gt).^2);
    % f = mean(abs(predictions/patient.VG - gt).^2)^(1/2);
    f = mean(abs(log(gt)-log(predictions/patient.VG)));
end

function window = set_window(tools, window_size, t, cgm_dt, x0, ymin1)
    window.dt = cgm_dt;
    experiment_total_time = length(tools.CGMs) * cgm_dt;
    
    window.t_start = t+1;
    window.t_end = min(t + window_size, experiment_total_time);

    window.x0 = x0;
    window.ymin1 = ymin1;

end

function stop = stop_simulation(tools, cgm_dt, t)
    experiment_total_time = length(tools.CGMs) * cgm_dt;
    if t >= experiment_total_time
        stop = true;
    else
        stop = false;
    end
end

