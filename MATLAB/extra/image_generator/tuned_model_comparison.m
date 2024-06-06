close all
clear
clc
rng default;

%% INIT
use_true_patient = false;
use_CGM_to_compare = true;

disp('Loading dataset...')
if use_true_patient
    use_anderson = true;
    use_tmoore = false;
    use_shanghai = false;
    plot_true_database = false;
    patient_ID = 11;
    dailyBasal = 18;
    date = '11-Feb-2013 06:30:00';
    % date = '26-Jan-2013 06:30:00';
    days_to_examine = 2;
    % days_to_examine = 30;
    % days_to_examine = 'all';
    start_day = datetime(date,'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
    tools = utils("none");
    run('database_preprocessor.m')
    params = patient_11(dailyBasal);
else
    filename = "/Users/giovannicherubini/Desktop/Thesis/Code/data/1minsample/adult#001_5.csv";
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
    params = patient_00(basal);
end
disp('Dataset loaded')


t_start = min([min(patientData.CGM.time), min(patientData.Meal.time), min(patientData.IIR.time)]);
t_end = max([max(patientData.CGM.time), max(patientData.Meal.time), max(patientData.IIR.time)]);
experiment_total_time = t_end - t_start;


% if use_tuned_model
%     if use_true_patient
%         run('true_p.m')
%     else
%         run('virtual_p.m')
%     end
% end


high_uncertainty = 10;
Q = eye(13) * high_uncertainty;
R = 100;

ekf_dt = 1;

model = non_linear_model(tools);
ekf = ekf(model, tools, params, ekf_dt, Q, R);

RMSE = [];
figure

for i = 1:3

    
    if i == 1

    else
        if use_true_patient
            run('true_p.m')
        else
            run('virtual_p.m')
        end
    end
    
    [x0, y_minus1] = tools.init_conditions(params);
    u0 = [0,0];
    
    EKF_state_tracking.mean = [];
    EKF_state_tracking.variance = [];
    EKF_state_tracking.time = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');
    
    gt = [];
    times_to_check = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');
    
    t = t_start;
    x = x0;
    y = y_minus1;
    u = u0;
    last_process_update = t_start;
    P0 = eye(13) * 12000;
    P = P0;
    
    %% Start simulation
    
    disp('Starting simulation...')
    while t <= t_end
        
        if use_CGM_to_compare
            [z_k, new_measurement_detected] = sample_measurement(t, patientData);
        else
            [z_k, new_measurement_detected] = sample_BG(t, patientData);
        end
    
        [u_k, new_input_detected] = sample_input(t, patientData);
        
        if new_input_detected || new_measurement_detected
    
            dt = convert_to_minutes(t - last_process_update);
            [xp_k, Pp_k, y_kminus1, v_kminus1] = ekf.predict(x, y, u, P, dt, params, true);
            x_current = xp_k;
            P_current = Pp_k;
    
            x = x_current;
            P = P_current;
            y = y_kminus1;
            u = u_k;
                
            last_process_update = t;
    
            if new_measurement_detected
                EKF_state_tracking.mean(:,end+1) = ekf.H * x';
                EKF_state_tracking.variance(end+1) = ekf.H * P * ekf.H';
                EKF_state_tracking.time(end+1, :) = t;
                gt(end+1) = z_k;
            end
        end
           
        t = next_step(t, patientData);
            
    end
    
    plot(EKF_state_tracking.time, EKF_state_tracking.mean, 'LineWidth', 1.5)
    hold on
    RMSE(end+1) = mean((EKF_state_tracking.mean - gt).^2)^(1/2);

end

RMSE
plot(EKF_state_tracking.time, gt, '-o', 'Color', 'cyan', 'LineWidth', 0.8)

legend 'Untuned' 'MAP-Tuned' 'PSO-Tuned' 'CGM'

xlim([t_start, t_end])
if use_true_patient
    title('True Patient - Model Comparison')
else
    title('In-Silico Patient - Model Comparison')
    ylim([40,250])
end
xlabel('Time')
ylabel('Plasma Glucose Concentration [mg/dL]')
set(gca, 'FontSize', 14)

disp("Done");
