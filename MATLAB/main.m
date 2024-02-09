close
clear
clc

parameters;
state_fields = {'Qsto1','Qsto2','Qgut','Gp','Gt','Gsc','Il','Ip','Id','I1','X','Isc1','Isc2'};
extra_state_fields = {'CHO_to_eat','D','lastQsto','is_eating'}; 
input_fields = {'CHO', 'IIR'};
true_input_fields = {'CHO_consumed_rate','IIR'};

filename = '';  % TBD

ekf_dt = 1;
cgm_dt = 3;
Q = eye(numel(state_fields)) * 0.08;    % TBD
R = 0.001;  % TBD

model = non_linear_model();
tools = utils(filename);
lin_model = linearized_model(params);
ekf = ekf(params, state_fields, R);
ekf.set_process_noise(Q);

[x, y] = tools.init_conditions(params);
P = eye(numel(state_fields));   % TBD

EKF_state_tracking.mean = [];
EKF_state_tracking.variance = [];
t1 = 10;
EKF_10_min_predictions.mean = [];
EKF_10_min_predictions.variance = [];
t2 = 30;
EKF_30_min_predictions.mean = [];
EKF_30_min_predictions.variance = [];

t = 0;
CGM_read = 0;
while ~stop_simulation(tools, CGM_read)
    new_measurement = false;

    % Get control input and eventually measurement
    if mod(t, cgm_dt) == 0
        CGM_read = CGM_read + 1;
        z = tools.CGMs(CGM_read);
        u.CHO = tools.CHOs(CGM_read);
        u.IIR = tools.IIRS(CGM_read);
        new_measurement = true;
    else
        u.CHO = 0;
        u.IIR = 0;
    end
    
    % Preprocess
    [y_current, v] = model.preprocess(x, y, u, params, ekf_dt);
    y = y_current;  % evaluate if updating it here

    % Given current state and input, linearize model
    lin_model = lin_model.linearize(x, y, params);
    ekf = ekf.update_matrices(lin_model.A, lin_model.D);

    % Process update
    [xp, Pp] = ekf.process_update(x,v,P);
    % ensure them to be non-negative    TBD
    
    % If there is measurement, do measurement update
    if new_measurement
        [xm, Pm] = ekf.measurement_update(xp, Pp, z);
        x_current = xm;
        P_current = Pm;
    else
        x_current = xp;
        P_current = Pp;
    end

    % Make predictions and store results

    % Update
    x = x_current;
    P = P_current;
    t = t + ekf_dt;
end

% Plot

function [x1, P1, x2, P2] = predict(x, y, t1, t2, ekf, lin_model)
    

end

function stop = stop_simulation(tools, CGM_read)
    stop = false;
    if CGM_read == length(tools.CGMs)
        stop = true;
    end
end