close all
clear
clc

state_fields = {'Qsto1','Qsto2','Qgut','Gp','Gt','Gpd','Il','Ip','I1','Id','X','Isc1','Isc2'};
extra_state_fields = {'insulin_to_infuse','last_IIR','CHO_to_eat','D','lastQsto','is_eating'}; 
input_fields = {'CHO', 'IIR'};
true_input_fields = {'CHO_consumed_rate','IIR_dt'};
filename = "/Users/giovannicherubini/Desktop/Thesis/Code/data/1minsample/adult#001.csv";

params_to_estimate = {'kp2','k1','k2','kp1','ki','ke1','kmax','kmin','kabs','kp3','Vmx'};
nvars = 11;
ub = [0.5,0.5,0.5,6,0.01,0.001,0.1,0.01,0.1,0.1,0.1];
lb = [0.0001,0.0001,0.0001,1,0.0001,0.0001,0.001,0.0001,0.001,0.0001,0.001];

tools = utils(filename, state_fields, extra_state_fields, input_fields, true_input_fields);
model = non_linear_model(tools);
basal = tools.IIRs(1);
patient = patient_00(basal);

cgm_dt = 1;
experiment_total_time = (length(tools.CGMs)-1)*cgm_dt;
window_size = experiment_total_time;
num_window = 1;

window = set_window(tools, window_size, num_window, cgm_dt, patient);

objective_pso = @(p) objective(p, patient, model, window, tools, params_to_estimate);

options.ObjectiveLimit = 2;
options.MaxTime = 1800;
options.Display = 'final';
[final_p,fval] = particleswarm(objective_pso, nvars, lb, ub, options);

disp('Done');

function f = objective(p, patient, model, window, tools, params_to_estimate)
    patient = patient.set_params(params_to_estimate,p);
    t = window.t_start;
    dt = window.dt;
    predictions = [];
    x = window.x0;
    y = window.ymin1;
    while t <= window.t_end
        u_k.CHO = window.CHOs(t);
        u_k.IIR = window.IIRs(t);

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

    % f = ( mean( ((predictions/patient.VG) - (transpose(window.BGs))).^2 ) )^(1/2)
    f = mean((predictions/patient.VG - transpose(window.BGs)).^2)

end

function window = set_window(tools, window_size, num_window, cgm_dt, patient)
    window.dt = cgm_dt;
    window.t_start = 1;
    window.t_end = length(tools.CGMs);
    window.CHOs = tools.CHOs(window.t_start:window.t_end);
    window.IIRs = tools.IIRs(window.t_start:window.t_end);
    window.BGs = tools.BGs(window.t_start:window.t_end);
    [window.x0, window.ymin1] = tools.init_conditions(patient);
end



