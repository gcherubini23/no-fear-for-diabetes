close
clear
clc

% parameters;
patient_01;
state_fields = {'Qsto1','Qsto2','Qgut','Gp','Gt','Gsc','Il','Ip','Id','I1','X','Isc1','Isc2'};
extra_state_fields = {'CHO_to_eat','D','lastQsto','is_eating'}; 
input_fields = {'CHO', 'IIR'};
true_input_fields = {'CHO_consumed_rate','IIR'};

filename = "/Users/giovannicherubini/Desktop/Thesis/Code/no-fear-for-diabetes/data/1minsample/adult#001.csv";

ekf_dt = 0.5;
cgm_dt = 1;
Q = eye(numel(state_fields)) * 0.08;    % TBD
R = 0.00;  % TBD

model = non_linear_model();
tools = utils(filename);
lin_model = linearized_model(params);
ekf = ekf(params, state_fields, ekf_dt, Q, R);

[x0, y0] = tools.init_conditions(params);
P = zeros(numel(state_fields),numel(state_fields));   % TBD

EKF_state_tracking.mean = [];
EKF_state_tracking.variance = [];
model_predictions = [];
u_history = [];
v_history = [];
y_history.CHO_to_eat = [];
y_history.D = [];
y_history.last_Q_sto = [];
y_history.is_eating = [];


t = 0;
CGM_read = 0;
x = x0;
y = y0;
x_true = x0;
IIRb = 0.021125;
while ~stop_simulation(tools, CGM_read)
    new_measurement = false;

    % Get control input and eventually measurement
    if mod(t, cgm_dt) == 0
        CGM_read = CGM_read + 1;
        z = tools.CGMs(CGM_read);
        u.CHO = tools.CHOs(CGM_read);
        u.IIR = tools.IIRs(CGM_read);
        new_measurement = true;
    else
        u.CHO = 0;
        u.IIR = 0;
    end
    
    % Preprocess
    [y_current, v] = model.preprocess(x, y, u, params, ekf_dt);
    y = y_current;

    u_history(end+1) = u.CHO;
    v_history(end+1) = v.CHO_consumed_rate;
    y_history.CHO_to_eat(end+1) = y.CHO_to_eat;
    y_history.D(end+1) = y.D;
    y_history.last_Q_sto(end+1) = y.lastQsto;
    y_history.is_eating(end+1) = y.is_eating;

    % Non-linear model
    dx_dt_new = model.step(x_true, y, v, params);
    x_true_new = ekf.euler_solve(x_true, dx_dt_new, ekf_dt);

    % EKF
    % Given current state and input, linearize model
    lin_model = lin_model.linearize(x, y, params);
    ekf = ekf.update_matrices(lin_model.A, lin_model.D);

    % Process update
    [xp, Pp] = ekf.process_update(x,v,P);

    % If there is measurement, do measurement update
    if new_measurement
        [xm, Pm] = ekf.measurement_update(xp, Pp, z);
        x_current = xm;
        P_current = Pm;
    else
        x_current = xp;
        P_current = Pp;
    end
    
    % % Test just linear model
    % x_current = xp;
    % P_current = Pp;


    % Store results
    EKF_state_tracking.mean(end+1) = x_current.Gp / params.VG;
    EKF_state_tracking.variance(end+1) = 1 / params.VG * sqrt(P_current(4,4));
    model_predictions(end+1) = x_true_new.Gp / params.VG;

    % Update
    x_true = x_true_new;
    x = x_current;
    P = P_current;
    t = t + ekf_dt;
end

disp("Done");

%% Plot

% Plot one figure with:

%   - 1 subplot that shows EKF_state_tracking.mean and model_predictions (timestep is 1 min)
%   - 1 subplot that shows tools.BGs (timestep is 3 min)
%   - 1 subplot that shows u_history and v_history (timestep is 1 min)

% Define time steps (dt) in minutes
dt_EKF_Model = ekf_dt; % Time step for EKF mean and model predictions
dt_BGs = cgm_dt; % Time step for BGs
dt_UV = ekf_dt; % Time step for u_history and v_history

% Assuming the total simulation time is known
totalMinutes = length(EKF_state_tracking.mean) * dt_EKF_Model;

% Time vectors based on respective dt and total simulation time
timeVecEKF_Model = 0:dt_EKF_Model:totalMinutes-dt_EKF_Model;
timeVecBGs = 0:dt_BGs:totalMinutes-dt_BGs;
timeVecUV = 0:dt_UV:totalMinutes-dt_UV;

% Number of points for each vector might differ, adjust according to the longest one
numPointsEKFModel = length(timeVecEKF_Model);
numPointsBGs = length(timeVecBGs);
numPointsUV = length(timeVecUV);

% Create figure
figure;

% Subplot 1: EKF_state_tracking.mean and model_predictions
subplot(3,1,1);
% plot(timeVecEKF_Model, EKF_state_tracking.mean(1:numPointsEKFModel), 'b-', 'DisplayName', 'EKF Mean');
% hold on;
plot(timeVecEKF_Model, model_predictions(1:numPointsEKFModel), 'g-', 'DisplayName', 'Model Predictions');
hold off;
legend show;
grid on;
% title('EKF Mean and Model Predictions');
title('Model Predictions');
xlabel('Time (min)');
ylabel('Value');

% Subplot 2: tools.BGs
subplot(3,1,2);
plot(timeVecBGs, tools.BGs(1:numPointsBGs), 'r-', 'DisplayName', 'BGs');
legend show;
grid on;
title('BGs');
xlabel('Time (min)');
ylabel('BG Value');

% Subplot 3: u_history and v_history
subplot(3,1,3);
plot(timeVecUV, u_history(1:numPointsUV), 'c-', 'DisplayName', 'u history');
hold on;
plot(timeVecUV, v_history(1:numPointsUV), 'm-', 'DisplayName', 'v history');
hold off;
legend show;
grid on;
title('Control Input History');
xlabel('Time (min)');
ylabel('Input Value');

% Adjust figure properties for better visibility
set(gcf, 'Position', get(0, 'Screensize')); % Maximize figure window for clarity


%% Extra functions

% function [x1, P1, x2, P2] = predict(x, y, t1, t2, ekf, lin_model)
% 
% 
% end

function stop = stop_simulation(tools, CGM_read)
    stop = false;
    if CGM_read == length(tools.CGMs)
        stop = true;
    end
end