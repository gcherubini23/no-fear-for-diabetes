close all
clear
clc

rng default;

%% Initialization

do_glucose = false;


% filename = "/Users/giovannicherubini/Desktop/Thesis/Code/data/1minsample/adult#001_3.csv";
filename = "/Users/giovannicherubini/Desktop/Thesis/Code/data/10minsample/1.csv";

tools = utils(filename);
patientData.CGM.values = tools.CGMs;
patientData.CGM.time = tools.Time;
patientData.Meal.values = tools.CHOs;
patientData.Meal.time = tools.Time;
patientData.IIR.values = tools.IIRs;
patientData.IIR.time = tools.Time;
patientData.BG.values = tools.BGs;
patientData.BG.time = tools.Time;

basal = tools.IIRs(1);
params = patient_01(basal);

figure

t_start = min([min(patientData.BG.time), min(patientData.Meal.time), min(patientData.IIR.time)]);
t_end = max([max(patientData.BG.time), max(patientData.Meal.time), max(patientData.IIR.time)]);
experiment_total_time = t_end - t_start;

%%
high_uncertainty = 0.1;
Q = eye(13) * high_uncertainty;
R = 100;

ekf_dt = [0.001, 0.01, 0.1, 0.5, 1, 2, 3, 3.7, 4, 4.5, 4.7, 5]; % [min]

% ekf_dt = [0.01, 0.1, 0.5, 1];

model = non_linear_model(tools);
ekf = ekf(model, tools, params, ekf_dt(1), Q, R);

alpha = 0.05;
anomaly_detector = anomaly_detector(alpha);

[x0, y_minus1] = tools.init_conditions(params);
u0 = [0,0];

RMSE = [];

for dt_idx = 1:length(ekf_dt)
    dte = ekf_dt(dt_idx);
    ekf.dt = dte;
    
    EKF_state_tracking.mean = [];
    EKF_state_tracking.variance = [];
    EKF_state_tracking.time = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');

    predictions_to_check = [];
    BG_to_check = [];
    times_to_check = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');

    t = t_start;
    x = x0;
    y = y_minus1;
    u = u0;
    last_process_update = t_start;
    P0 = eye(13) * 12000;
    P = P0;

    flag = true;
    disp('Starting simulation...')

    while t <= t_end
        [z_k, new_measurement_detected] = sample_BG(t, patientData);
        [u_k, new_input_detected] = sample_input(t, patientData);
        
        if new_measurement_detected || new_input_detected
    
            dt = convert_to_minutes(t - last_process_update);
            
            [xp_k, Pp_k, y_kminus1, v_kminus1] = ekf.predict(x, y, u, P, dt, params, true);
            
            x_current = xp_k;
            P_current = Pp_k;
                       
            x = x_current;
            P = P_current;
            y = y_kminus1;
            u = u_k * dt;

            predictions_to_check(end+1) = ekf.H * x';
            BG_to_check(end+1) = z_k;
                
            last_process_update = t;
    
            if do_glucose
                EKF_state_tracking.mean(:,end+1) = ekf.H * x';
            else
                EKF_state_tracking.mean(:,end+1) = x(8) / params.VI;
            end
            EKF_state_tracking.time(end+1, :) = t;
    
        end 
        t = next_step(t, patientData);
    end
    disp('Done');

    plot(EKF_state_tracking.time, EKF_state_tracking.mean, '-o', 'LineWidth', 1.5);
    hold on


    if dt_idx == 1
        first = predictions_to_check;
    else
        RMSE(end+1) = mean((predictions_to_check - first).^2)^(1/2);
    end

    % RMSE(end+1) = mean((predictions_to_check - BG_to_check).^2)^(1/2);


end

% plot(patientData.BG.time, patientData.BG.values, 'DisplayName', 'BG');
legend '\Deltat = 0.001' '\Deltat = 0.01' '\Deltat = 0.1' '\Deltat = 0.5' '\Deltat = 1' '\Deltat = 2' '\Deltat = 3' '\Deltat = 3.7' '\Deltat = 4' '\Deltat = 4.5' '\Deltat = 4.7' '\Deltat = 5' '\Deltat = 10'

if do_glucose
    ylim([0, 220])
else
    ylim([80, 400])
end

xlabel('Time')
if do_glucose
    ylabel('Plasma Glucose Concentration [mg/dL]')
    title('Time Step Effect on Glucose')
else
    ylabel('Plasma Insulin Concentration [pmol/L]')
    title('Time Step Effect on Insulin')
end

set(gca, 'FontSize', 14)

RMSE
