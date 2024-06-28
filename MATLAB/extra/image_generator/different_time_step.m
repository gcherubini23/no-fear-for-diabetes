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

nonZeroIndex = tools.CHOs ~= 0;
patientData.Meal.values = tools.CHOs(nonZeroIndex);
patientData.Meal.time = tools.Time(nonZeroIndex);

patientData.IIR.values = tools.IIRs;
patientData.IIR.time = tools.Time;
patientData.BG.values = tools.BGs;
patientData.BG.time = tools.Time;

basal = tools.IIRs(1);
params = patient_01(basal);

t_start = min([min(patientData.BG.time), min(patientData.Meal.time), min(patientData.IIR.time)]);
t_end = max([max(patientData.BG.time), max(patientData.Meal.time), max(patientData.IIR.time)]);
experiment_total_time = t_end - t_start;


%%

figure

ms = 2;
lw = 0.5;

font_size = 10;
figureWidth = 5.7;
figureHeight = 8;
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, figureWidth, figureHeight]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figureWidth, figureHeight]);
set(gcf, 'PaperPositionMode', 'auto');

lables = {'\Deltat = 0.001 min' '\Deltat = 0.01 min' '\Deltat = 0.1 min' '\Deltat = 0.5 min' '\Deltat = 1 min' '\Deltat = 2 min' '\Deltat = 3 min' '\Deltat = 3.7 min' '\Deltat = 4 min' '\Deltat = 4.5 min' '\Deltat = 4.7 min' '\Deltat = 5 min' '\Deltat = 10 min'};

rainbowColors = [
    1.0, 0.0, 0.0;     % Red
    1.0, 0.2, 0.0;     % Red-Orange
    1.0, 0.4, 0.0;     % Orange
    1.0, 0.6, 0.0;     % Orange-Yellow
    1.0, 0.8, 0.0;     % Yellow
    0.8, 1.0, 0.0;     % Yellow-Green
    0.6, 1.0, 0.0;     % Green
    0.0, 1.0, 0.6;     % Green-Cyan
    0.0, 1.0, 1.0;     % Cyan
    0.0, 0.6, 1.0;     % Cyan-Blue
    0.0, 0.0, 1.0;     % Blue
    0.6, 0.0, 1.0;     % Blue-Violet
    117/255, 0, 214/255 
];


% subfigures = 7;
subfigures = 7 * 4;

%%
high_uncertainty = 0.1;
Q = eye(13) * high_uncertainty;
R = 100;

ekf_dt = [0.001, 0.01, 0.1, 0.5, 1, 2, 3, 3.7, 4, 4.5, 4.7, 5]; % [min]
% ekf_dt = [0.001, 0.01, 0.1, 0.5, 1, 2, 3, 3.7, 4, 4.5, 4.7, 5, 10]; % [min]


model = non_linear_model(tools);
ekf = ekf(model, tools, params, ekf_dt(1), Q, R);

alpha = 0.05;
anomaly_detector = anomaly_detector(alpha);

[x0, y_minus1] = tools.init_conditions(params);
u0 = [0,0];

RMSE = [];

imagesG = [];
imagesI = [];

EKF_insulin = [];

subplot(subfigures, 1, [1 13]);

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
    
            
            EKF_state_tracking.mean(:,end+1) = x';
            
                % EKF_state_tracking.mean(:,end+1) = x(8) / params.VI;
            
            EKF_state_tracking.time(end+1, :) = t;
    
        end 
        t = next_step(t, patientData);
    end
    disp('Done');

    imagesG(end+1) = plot(EKF_state_tracking.time, EKF_state_tracking.mean(6,:)/params.VG , '-o', 'LineWidth', lw, 'MarkerSize', ms, 'Color', rainbowColors(dt_idx,:), 'MarkerFaceColor', rainbowColors(dt_idx,:));
    hold on
 
    % EKF_insulin(end+1,:) = EKF_state_tracking.mean(8,:)/params.VI;
    EKF_insulin(end+1,:) = EKF_state_tracking.mean(7,:);

    if dt_idx == 1
        first = predictions_to_check;
    else
        RMSE(end+1) = mean((predictions_to_check - first).^2)^(1/2);
    end

    % RMSE(end+1) = mean((predictions_to_check - BG_to_check).^2)^(1/2);


end

% plot(patientData.BG.time, patientData.BG.values, 'DisplayName', 'BG');


% if do_glucose
%     ylim([0, 220])
% else
%     ylim([80, 400])
% end

lgd = legend(imagesG, lables);
lgd.Location = 'northoutside';
lgd.Orientation = 'horizontal';
lgd.Box = 'off';
lgd.NumColumns = 3;
title('Discretization Numerical Instability')
ylim([0, 220])
ylabel('G [mg/dL]')

% if do_glucose
%     ylabel('G [mg/dL]')
%     title('Time Step Effect on Glucose')
% else
%     ylabel('I [pmol/L]')
%     title('Time Step Effect on Insulin')
% end

subplot(subfigures, 1, [16 23]);

for ig = 1:length(ekf_dt)
    imagesI(end+1) = plot(EKF_state_tracking.time, EKF_insulin(ig,:), '-o', 'LineWidth', lw, 'MarkerSize', ms, 'Color', rainbowColors(ig,:), 'MarkerFaceColor', rainbowColors(ig,:));
    hold on
end
% ylim([80, 400])
ylim([0, 11])
ylabel('I_l [pmol/kg]')

RMSE

subplot(subfigures, 1, [25 28]);  % Allocate 1/3 space for the subplot
yyaxis right;
u2 = stem(patientData.IIR.time, patientData.IIR.values, 'filled', 'o', 'Color', [251,142,2]/255,'MarkerSize', ms, 'LineWidth', lw);
ylabel('Insulin [IU]');
% ylim([0,7])
xlim([t_start, t_end])
hold on
ax = gca;
ax.YColor = 'k'; % Set color of left y-axis to black

yyaxis left;
u1 = stem(patientData.Meal.time, patientData.Meal.values, 'filled', 'square', 'Color', 'g','MarkerSize', ms, 'LineWidth', lw);
ylabel('CHO [g]');
% ylim([0,60])
xlim([t_start, t_end])
ax = gca;
ax.YColor = 'k';

xlabel('Time');
lgd2 = legend([u1,u2], {'Meal', 'Insulin infused'});
lgd2.Location = 'northoutside';
lgd2.Orientation = 'horizontal';
lgd2.Box = 'off';


set(findall(gcf, '-property', 'FontSize'), 'FontSize', font_size);


