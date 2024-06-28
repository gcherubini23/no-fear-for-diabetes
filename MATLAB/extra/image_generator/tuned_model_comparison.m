close all
clear
clc
rng default;

%% INIT
use_true_patient = false;
use_CGM_to_compare = true;
plot_model = true;

ms = 2;
lw = 0.5;

figure

if plot_model
    if use_true_patient
        lables = {'CGM' 'Untuned' 'MAP-Tuned' 'PSO-Tuned'};
    else
        lables = {'CGM' 'True' 'Untuned' 'MAP-Tuned' 'PSO-Tuned'};
    end
else
    lables = {'CGM'};
end


font_size = 10;
figureWidth = 5.7;
figureHeight = 4;
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, figureWidth, figureHeight]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figureWidth, figureHeight]);
set(gcf, 'PaperPositionMode', 'auto');


disp('Loading dataset...')
if use_true_patient
    use_anderson = true;
    use_tmoore = false;
    use_shanghai = false;
    plot_true_database = false;
    patient_ID = 11;
    dailyBasal = 16.9;
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


high_uncertainty = 10;
Q = eye(13) * high_uncertainty;
R = 100;

ekf_dt = 1;

model = non_linear_model(tools);
ekf = ekf(model, tools, params, ekf_dt, Q, R);

RMSE = [];
images = [];

colors = [76, 8, 243;
          199, 0, 57;
          15, 169, 46]/255;

subplot(4, 1, [1 3]);

images(end+1) = plot(patientData.CGM.time, patientData.CGM.values, '-o', 'Color', 'cyan', 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', 'cyan');
hold on

if ~use_true_patient
    images(end+1) = plot(patientData.BG.time, patientData.BG.values, '-', 'Color', [255, 213, 0]/255, 'LineWidth', 2*lw);
end

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
    
    if plot_model
        images(end+1) = plot(EKF_state_tracking.time, EKF_state_tracking.mean, 'LineWidth', 2*lw, 'Color', colors(i,:));
        hold on
    end
    RMSE(end+1) = mean((EKF_state_tracking.mean - gt).^2)^(1/2);

end

RMSE
   
lgd = legend(images, lables);
lgd.Location = 'northoutside';
lgd.Orientation = 'horizontal';
lgd.Box = 'off';

if ~use_true_patient
    lgd.NumColumns = 3;
end


xlim([t_start, t_end])
if plot_model
    if use_true_patient
        ylim([40,300])
    else
        ylim([40,250])
    end
end
ylabel('G [mg/dL]')


if use_true_patient && plot_model
    title('True Patient - Model Comparison')
elseif ~use_true_patient && plot_model
    title('Simulated Patient - Model Comparison')
else
    title('True Patient one-day data')
end

%%

subplot(4, 1, 4);  % Allocate 1/3 space for the subplot
yyaxis right;
u2 = stem(patientData.IIR.time, patientData.IIR.values, 'filled', 'o', 'Color', [251,142,2]/255,'MarkerSize', 2, 'LineWidth', lw);
ylabel('Insulin [IU]');
if use_true_patient
    ylim([0,7])
else
    ylim([0,14])
end
xlim([t_start, t_end])
hold on
ax = gca;
ax.YColor = 'k'; % Set color of left y-axis to black

yyaxis left;
u1 = stem(patientData.Meal.time, patientData.Meal.values, 'filled', 'square', 'Color', 'g','MarkerSize', 2, 'LineWidth', lw);
ylabel('CHO [g]');
if use_true_patient
    ylim([0,60])
else
    ylim([0,90])
end
% xlim([t_start + seconds(1), t_end-seconds(1)])
xlim([t_start, t_end])
ax = gca;
ax.YColor = 'k';

xlabel('Time');
lgd2 = legend([u1,u2], {'Meal', 'Insulin infused'});
lgd2.Location = 'northoutside';
lgd2.Orientation = 'horizontal';
lgd2.Box = 'off';


set(findall(gcf, '-property', 'FontSize'), 'FontSize', font_size);

