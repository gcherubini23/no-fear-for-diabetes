close all
clear
clc

state_fields = {'Qsto1','Qsto2','Qgut','Gp','Gt','Gpd','Il','Ip','I1','Id','X','Isc1','Isc2'};
extra_state_fields = {'insulin_to_infuse','last_IIR','CHO_to_eat','D','lastQsto','is_eating'}; 
input_fields = {'CHO', 'IIR'};
true_input_fields = {'CHO_consumed_rate','IIR_dt'};
filename = "/Users/giovannicherubini/Desktop/Thesis/Code/data/1minsample/adult#001.csv";

params_to_estimate = {'k1'};
nvars = 11;
ub = [];
lb = [];

tools = utils(filename, state_fields, extra_state_fields, input_fields, true_input_fields);
model = non_linear_model(tools);
basal = tools.IIRs(1);
patient = patient_00(basal);

cgm_dt = 1;
experiment_total_time = (length(tools.CGMs)-1)*cgm_dt;
time_window = experiment_total_time;

history = 0;

objective_pso = @(p) objective(p, history, patient, model, params_to_estimate);

[final_p,fval] = particleswarm(objective_pso, nvars, lb, ub);


function f = objective(p, patient, model, history, tools, params_to_estimate)
    patient = patient.set_params(params_to_estimate,p);
    t = history.t_start;
    predictions = [];
    x = history.x0;
    y = history.ymin1;
    while t <= history.t_end
        
        
        t = t + history.dt;
    end
end
