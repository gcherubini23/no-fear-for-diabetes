close all
clear
clc

state_fields = {'Qsto1','Qsto2','Qgut','Gp','Gt','Gpd','Il','Ip','I1','Id','X','Isc1','Isc2'};
extra_state_fields = {'insulin_to_infuse','last_IIR','CHO_to_eat','D','lastQsto','is_eating'}; 
input_fields = {'CHO', 'IIR'};
true_input_fields = {'CHO_consumed_rate','IIR_dt'};

simulate_anomalies = false;
use_tuned_model = false;
use_true_model = true;
use_known_init_conditions = true;
do_measurment_update = false;
do_chi_sq_test = false;
do_cusum_test = false;
do_plots = true;

filename = "/Users/giovannicherubini/Desktop/Thesis/Code/data/1minsample/adult#001_4.csv";

Q = eye(numel(state_fields)) * 15;    % TBD
R = 100;  % TBD

% ekf_dt_values = [0.1, 0.5, 1];
ekf_dt_values = [1];
ekf_dt = ekf_dt_values(1);
cgm_dt = 1;
pump_dt = 1;

tools = utils(filename, state_fields, extra_state_fields, input_fields, true_input_fields);
basal = tools.IIRs(1);

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
else
    [x0, y_minus1] = tools.rand_conditions(params);
end

% Loop over the values of ekf_dt
prediction.first = [];
prediction.second = [];
prediction.third = [];
mse = [];
rmse = [];
i = 1;
for ekf_dt = ekf_dt_values
    
    ekf.dt = ekf_dt;
    model_predictions = [];
    residuals = [];
    innovation_cov = [];
    EKF_state_tracking.mean = [];
    EKF_state_tracking.variance = [];
    z_history = [];
    u_history = [];
    v_history = [];
    y_history.CHO_to_eat = [];
    y_history.D = [];
    y_history.last_Q_sto = [];
    y_history.is_eating = [];
    y_history.last_IIR = [];
    y_history.insulin_to_infuse = [];
    
    k = 0;
    t = 0;
    x = x0;
    y = y_minus1;
    P = eye(numel(state_fields),numel(state_fields))*2200;

    while ~stop_simulation(tools, t)
        
        t = k * ekf_dt;
        [z_k, new_measurement_detected] = sample_measurement(t, tools, cgm_dt);
        [u_k, new_input_detected] = sample_input(t, tools, pump_dt);

        % EKF
        if k == 0   % if starting now, set initial conditions
            xp_k = x;
            Pp_k = P;
            y_kminus1 = y;
        else
            [xp_k, Pp_k, y_kminus1, v_kminus1] = ekf.process_update(x, y, u, P, ekf.dt, params); % processupdate(x_k-1, y_k-2, u_k-1)
            v_history(end+1) = v_kminus1.CHO_consumed_rate;
        end

        x_current = xp_k;
        P_current = Pp_k;

        if do_measurment_update && new_measurement_detected
            z_history(end+1) = z_k;
            [xm_k, Pm_k, residual_k, innov_cov_k] = ekf.measurement_update(xp_k,Pp_k,z_k);
            x_current = xm_k;
            P_current = Pm_k;
            residuals(end+1) = residual_k;
            innovation_cov(end+1) = innov_cov_k;
        end
        

        % Update
        x = x_current;
        P = P_current;
        y = y_kminus1;
        u = u_k;
        k = k + 1;
        
        % Store results
        model_predictions(:,end+1) = tools.convert_to_vector(x);   % (time k)
        u_history(end+1) = u.CHO;   % time k

        EKF_state_tracking.variance(end+1) = P_current(6,6);
        
        if k > 1
            y_history.CHO_to_eat(end+1) = y.CHO_to_eat;
            y_history.D(end+1) = y.D;
            y_history.last_Q_sto(end+1) = y.lastQsto;
            y_history.is_eating(end+1) = y.is_eating;
            y_history.insulin_to_infuse(end+1) = y.insulin_to_infuse;
            y_history.last_IIR(end+1) = y.last_IIR;

        end
        
        % pause
    end
    
    v_history(end+1) = u.CHO;
    y_history.CHO_to_eat(end+1) = y.CHO_to_eat;
    y_history.D(end+1) = y.D;
    y_history.last_Q_sto(end+1) = y.lastQsto;
    y_history.is_eating(end+1) = y.is_eating;
    y_history.insulin_to_infuse(end+1) = y.insulin_to_infuse;
    y_history.last_IIR(end+1) = y.last_IIR;

    if length(ekf_dt_values) > 1
        if i == 1
            prediction.first = model_predictions;
        elseif i == 2
            prediction.second = model_predictions;
            
        else
            prediction.third = model_predictions;
        end
        i = i+1;
    end

    mse(end+1) = mean((transpose(tools.BGs)-(model_predictions(6,1:1/ekf_dt:end)/params.VG)).^2);
    rmse(end+1) = mean( ((model_predictions(6,1:1/ekf_dt:end)/params.VG) - (transpose(tools.BGs))).^2 )^(1/2);    % RMSE

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
    run('plotting_old.m');
end

%% Extra functions

function stop = stop_simulation(tools, CGM_read)

    stop = false;
    if CGM_read == length(tools.CGMs)-1
        stop = true;
    end
end

function [z_k, new_measurement_detected] = sample_measurement(t, tools, cgm_dt)
    new_measurement_detected = false;
    z_k = [];
    if mod(t, cgm_dt) == 0
        CGM_read = t + 1;
        z_k = tools.CGMs(CGM_read);
        new_measurement_detected = true;
    end
end

function [u_k, new_input_detected] = sample_input(t, tools, pump_dt)
    new_input_detected = false;
    u_k.CHO = 0;
    u_k.IIR = 0;
    if mod(t, pump_dt) == 0
        U_read = t + 1;
        u_k.CHO = tools.CHOs(U_read);
        u_k.IIR = tools.IIRs(U_read);
        new_input_detected = true;
    end
end