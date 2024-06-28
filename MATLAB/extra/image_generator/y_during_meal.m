close all
clear
clc

rng default;

%% Initialization

do_meal = false;

run('config.m')

if use_tuned_model
    run('param_model.m')
end

horizon = 30;
MARD = 10;
high_uncertainty = 30;

Q = eye(13) * high_uncertainty;
R = 100;

if do_meal
    ekf_dt = 0.5; % insulin
else
    ekf_dt = 0.05; % insulin
end


model = non_linear_model(tools);
ekf = ekf(model, tools, params, ekf_dt, Q, R);
ekf.dt = ekf_dt;
ekf.CGM_MARD = MARD;

alpha = 0.05;
anomaly_detector = anomaly_detector(alpha);

if use_basal_init_conditions
    [x0, y_minus1] = tools.init_conditions(params);
    u0 = [0,0];
else
    [x0, y_minus1] = tools.set_init_conditions(patientData.CGM.values(1), params);
    u0 = [0,0];
end

EKF_state_tracking.mean = [];
EKF_state_tracking.variance = [];
EKF_state_tracking.time = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');
future_predictions.values = [];
future_predictions.cov = [];
future_predictions.time = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');
residuals = [];

predictions_to_check = [];
CGM_to_check = [];
times_to_check = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');

trajectories = {};

t = t_start;
x = x0;
y = y_minus1;
u = u0;
last_process_update = t_start;
P0 = eye(13) * 12000;
P = P0;

%% Start simulation
% profile on

disp('Starting simulation...')
while t <= t_end
    
    [z_k, new_measurement_detected] = sample_measurement(t, patientData);
    [u_k, new_input_detected] = sample_input(t, patientData);
    
    if new_input_detected || new_measurement_detected

        dt = convert_to_minutes(t - last_process_update);
        [xp_k, Pp_k, y_kminus1, v_kminus1, ekf] = ekf.predict_save_all_variables(x, y, u, P, dt, params, last_process_update, true);
        
        x_current = xp_k;
        P_current = Pp_k;
                    
        if new_measurement_detected && do_measurment_update
            ekf = ekf.update_sensor_cov(z_k);
            [xm_k, Pm_k, residual_k, innov_cov_k] = ekf.measurement_update(xp_k,Pp_k,z_k);
            residuals(end+1) = residual_k;
            if do_chi_sq_test || do_cusum_test
                if do_chi_sq_test
                    anomaly = anomaly_detector.chi_squared_test(residual_k,innov_cov_k);
                else
                    [anomaly, anomaly_detector] = anomaly_detector.cusum_test(residual_k,innov_cov_k);
                end
                if ~anomaly
                    x_current = xm_k;
                    P_current = Pm_k;
                else
                    anomaly_detector.time(end+1, :) = t;
                    anomaly_detector.anomalies(end+1) = z_k;
                end
            else
                x_current = xm_k;
                P_current = Pm_k;
                anomaly = false;
            end
        end

        x = x_current;
        P = P_current;
        y = y_kminus1;
        u = u_k;
            
        last_process_update = t;

        EKF_state_tracking.mean(:,end+1) = x';
        EKF_state_tracking.variance(end+1) = ekf.H * P * ekf.H';
        EKF_state_tracking.time(end+1, :) = t;
    end
       
    t = next_step(t, patientData);
        
end
% profile off
% profile viewer


%% Plot ptred and track

ms = 2;
lw = 0.5;
font_size = 10;

if do_meal
    figure

    images = [];

    figureWidth = 5.7;
    figureHeight = 3.5;
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [1, 1, figureWidth, figureHeight]);
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [figureWidth, figureHeight]);
    set(gcf, 'PaperPositionMode', 'auto');


    lables1 = {'Meal size', 'CHO to eat', 'D', 'Model input (CHO/min)'};
    lables3 = {'Is eating'};

    ts = datetime('11-Mar-2024 08:35:30');
    te = datetime('11-Mar-2024 08:49:30');
    subplot(3,1,[1 2]);
    images = [];
    i = find(ekf.u_history(1,:)~=0);
    images(end+1) = stem(ekf.t_history(i), ekf.u_history(1,i), 'filled', 'Color', 'g','MarkerSize', 4*ms, 'LineWidth', 2*lw);
    hold on

    f = find(ekf.u_history(1,:)>0);
    ekf.y_history(3,f) = ekf.u_history(1,f);
    iii = find(ekf.y_history(3,:)~=0);
    ii = find(ekf.t_history <= te+minutes(3));
    i = intersect(iii,ii);
    images(end+1) = stem(ekf.t_history(i),ekf.y_history(3,i), 'LineWidth', 2*lw, 'Color', 'm');
    hold on
    i = find(ekf.y_history(4,:)~=0);
    images(end+1) = stem(ekf.t_history(i),ekf.y_history(4,i), 'LineWidth', 2*lw, 'Color', [255, 184, 0]/255, 'Marker', 'x');
    hold on
    i = find(ekf.v_history(1,:)~=0);
    images(end+1) = stem(ekf.t_history(i)-minutes(ekf_dt),ekf.v_history(1,i), 'filled', 'o', 'Color', 'b','MarkerSize', 1.5*ms, 'LineWidth', lw);
    hold on
    xlim([ts te])
    ylabel('CHO [g]')
    lgd = legend(images, lables1);
    lgd.Location = 'northoutside';
    lgd.Orientation = 'horizontal';
    lgd.Box = 'off';
    lgd.NumColumns = 2;
    
    subplot(3,1,3);
    images = [];
    images(end+1) = plot(ekf.t_history-minutes(ekf_dt),ekf.y_history(6,:), 'o', 'MarkerFaceColor','auto', 'MarkerSize', 2*ms, 'Color', [118,0,200]/255, 'MarkerFaceColor', [118,0,200]/255);
    xlim([ts te])
    lgd = legend(images, lables3);
    lgd.Location = 'northoutside';
    lgd.Orientation = 'horizontal';
    lgd.Box = 'off';
    xlabel('Time');

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', font_size);
else
    figure

    lables2 = {'Injected insulin', 'Model input (IU/min)', 'Insulin to infuse'};

    images = [];

    % ts = datetime('11-Mar-2024 08:37:45');
    % te = datetime('11-Mar-2024 08:39:15');
    
    ts = datetime('11-Mar-2024 08:37:30');
    te = datetime('11-Mar-2024 08:39:30');


    figureWidth = 5.7;
    figureHeight = 2.2;
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [1, 1, figureWidth, figureHeight]);
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [figureWidth, figureHeight]);
    set(gcf, 'PaperPositionMode', 'auto');
    
    i = find(ekf.u_history(2,:)~=0);
    images(end+1) = stem(ekf.t_history(i), ekf.u_history(2,i), 'filled', 'o', 'Color', 'g','MarkerSize', 4*ms, 'LineWidth', 2*lw);
    hold on
    i = find(ekf.v_history(2,:)~=0);
    images(end+1) = stem(ekf.t_history(i)-minutes(ekf_dt),ekf.v_history(2,i), 'filled', 'o', 'Color', 'b','MarkerSize', 1.5*ms, 'LineWidth', lw);
    hold on
    % i = find(ekf.u_history(2,:) > 0);
    % ekf.y_history(1,i) = ekf.u_history(2,i);
    i = find(ekf.y_history(1,:)~=0);
    images(end+1) = stem(ekf.t_history(i)-minutes(ekf_dt),ekf.y_history(1,i), 'LineWidth', 2*lw, 'Color', 'm');
    hold on

    xlim([ts te])
    ylabel('Insulin [IU]')
    lgd = legend(images, lables2);
    lgd.Location = 'northoutside';
    lgd.Orientation = 'horizontal';
    lgd.Box = 'off';
    xlabel('Time');

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', font_size);

end
