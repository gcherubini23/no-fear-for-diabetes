close all
clear
clc

state_fields = {'Qsto1','Qsto2','Qgut','Gp','Gt','Gpd','Il','Ip','I1','Id','X','Isc1','Isc2'};
extra_state_fields = {'insulin_to_infuse','last_IIR','CHO_to_eat','D','lastQsto','is_eating'}; 
input_fields = {'CHO', 'IIR'};
true_input_fields = {'CHO_consumed_rate','IIR_dt'};

filename = "/Users/giovannicherubini/Desktop/Thesis/Code/data/1minsample/adult#001_3.csv";

Q = eye(numel(state_fields)) * 1;    % TBD
R = 20;  % TBD

% ekf_dt_values = [0.1, 0.5, 1];
ekf_dt_values = [1];
ekf_dt = ekf_dt_values(1);
cgm_dt = 1;

tools = utils(filename, state_fields, extra_state_fields, input_fields, true_input_fields);
model = non_linear_model(tools);
lin_model = linearized_model(tools);
ekf = ekf(model, lin_model, tools, ekf_dt, Q, R);

basal = tools.IIRs(1);

params = patient_01(basal);

% params = patient_00(basal);
% % Best using mse
% p_pso_00 = [0.00613673067477461,	0.400287164971832,	0.494617050929902,	4.92072844013362,	0.00382498896706389,	0.000826933211308203,	0.0998386175768517,	0.00854364392451055,	0.0997034961579507,	0.0309196769142360,	0.0389878346812570];
% % Best using rmse
% p_pso_00 = [0.00726798039653122,	0.288211953663130,	0.357087874454711,	5.50484148758251,	0.00351366020364207,	0.000786244943534256,	0.0999827530364573,	0.00890161075502704,	0.0999958739540014,	0.0365651255883235,	0.0360709498914858];
% params = params.set_params({'kp2','k1','k2','kp1','ki','ke1','kmax','kmin','kabs','kp3','Vmx'}, p_pso_00);

[x0, y_minus1] = tools.init_conditions(params);
% [x0, y_minus1] = tools.rand_conditions(params);

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
    ground_truth = [];
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
        [z_k, new_measurement_detected] = sample_measurement(t, tools);
        [u_k, new_input_detected] = sample_input(t, tools);

        % EKF

        if i == 1
            what = true;
        else
            what = false;
        end
        
        if length(ekf_dt_values) == 3
            what = true;
        end

        if k == 0   % if starting now, set initial conditions
            xp_k = x;
            Pp_k = P;
            y_kminus1 = y;
        else
            [xp_k, Pp_k, y_kminus1, v_kminus1] = ekf.process_update(x, y, u, P, params, what, t); % processupdate(x_k-1, y_k-2, u_k-1)
            v_history(end+1) = v_kminus1.CHO_consumed_rate;
        end

        x_current = xp_k;
        P_current = Pp_k;

        if false && new_measurement_detected
            z_history(end+1) = z_k;
            [xm_k, Pm_k] = ekf.measurement_update(xp_k,Pp_k,z_k,params);
            x_current = xm_k;
            P_current = Pm_k;
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
        
        
        if k > 1
            y_history.CHO_to_eat(end+1) = y.CHO_to_eat;
            y_history.D(end+1) = y.D;
            y_history.last_Q_sto(end+1) = y.lastQsto;
            y_history.is_eating(end+1) = y.is_eating;
            y_history.insulin_to_infuse(end+1) = y.insulin_to_infuse;
            y_history.last_IIR(end+1) = y.last_IIR;

        end
        
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

    if ekf_dt == 1
        mse(end+1) = mean((transpose(tools.BGs)-(model_predictions(6,:)/params.VG)).^2);
        rmse(end+1) = mean( ((model_predictions(6,:)/params.VG) - (transpose(tools.BGs))).^2 )^(1/2);    % RMSE
    end

    disp("Done");
end

%% Plot

if length(ekf_dt_values) == 1 && true
    if ekf_dt <= 1
        dt_EKF_Model = ekf_dt; % Time step for EKF mean and model predictions
        dt_BGs = cgm_dt; % Time step for BGs
        dt_UV = ekf_dt; % Time step for u_history and v_history
    else
        dt_EKF_Model = ekf_dt; % Time step for EKF mean and model predictions
        dt_BGs = ekf_dt; % Time step for BGs
        dt_UV = ekf_dt; % Time step for u_history and v_history
    end
    
    % Assuming the total simulation time is known
    totalMinutes = length(model_predictions) * ekf_dt;
    
    % Time vectors based on respective dt and total simulation time
    timeVecEKF_Model = 0:dt_EKF_Model:totalMinutes-dt_EKF_Model;
    timeVecBGs = 0:dt_BGs:totalMinutes-dt_BGs;
    timeVecUV = 0:dt_UV:totalMinutes-dt_UV;
    
    % Number of points for each vector might differ, adjust according to the longest one
    numPointsEKFModel = length(timeVecEKF_Model);
    numPointsBGs = length(timeVecBGs);
    numPointsUV = length(timeVecUV);
    
    if (true)
        
        % Create figure
        figure;
        
        % Subplot 1: EKF_state_tracking.mean and model_predictions
        subplot(1,1,1);
        plot(timeVecBGs, tools.BGs(1:numPointsBGs), 'g-', 'DisplayName', 'True Gpd');
        hold on
        plot(timeVecEKF_Model, model_predictions(6,1:numPointsEKFModel)/params.VG, 'b-', 'DisplayName', 'EKF Gpd');
        hold on
        stem(timeVecBGs, tools.CGMs(1:numPointsBGs), 'x', 'MarkerSize', 6, 'LineStyle', 'none', 'DisplayName', 'CGM')
        hold off;
        legend show;
        grid on;
        % title('EKF Mean and Model Predictions');
        % title('Model Predictions');
        xlabel('Time (min)');
        ylabel('Value');
        
        set(gcf, 'Position', get(0, 'Screensize')); % Maximize figure window for clarity
    end
   
    % if (false)
    %     % Create figure
    %     figure;
    % 
    %     % Subplot 1
    %     subplot(3,1,1);
    %     plot(timeVecEKF_Model, model_predictions(1,1:numPointsEKFModel), 'g-', 'DisplayName', 'Qsto1');
    %     legend show;
    %     grid on;
    %     xlabel('Time (min)');
    %     ylabel('Value');
    % 
    %     % Subplot 2
    %     subplot(3,1,2);
    %     plot(timeVecEKF_Model, model_predictions(2,1:numPointsEKFModel), 'r-', 'DisplayName', 'Qsto2');
    %     legend show;
    %     grid on;
    %     xlabel('Time (min)');
    %     ylabel('Value');
    % 
    %     % Subplot 3
    %     subplot(3,1,3);
    %     plot(timeVecEKF_Model, model_predictions(3,1:numPointsEKFModel), 'b-', 'DisplayName', 'Qgut');
    %     legend show;
    %     grid on;
    %     xlabel('Time (min)');
    %     ylabel('Value');
    % 
    % 
    % end
    
    
    % if (false)
    % 
    %     % Create figure
    %     figure;
    % 
    %     % Subplot 1
    %     subplot(2,1,1);
    %     plot(timeVecEKF_Model, y_history.D(1:numPointsUV), 'g-', 'DisplayName', 'D');
    %     hold on;
    %     plot(timeVecEKF_Model, u_history(1:numPointsUV), 'b-', 'DisplayName', 'u');
    %     hold on
    %     plot(timeVecEKF_Model, v_history(1:numPointsUV), 'r-', 'DisplayName', 'v');
    %     hold off;
    %     legend show;
    %     grid on;
    %     xlabel('Time (min)');
    %     ylabel('Value');
    % 
    %     % xlim([418,426]);
    % 
    %     % Subplot 2
    %     subplot(2,1,2);
    %     plot(timeVecEKF_Model, y_history.CHO_to_eat(1:numPointsUV), 'm-', 'DisplayName', 'CHO_to_eat');
    %     legend show;
    %     grid on;
    %     xlabel('Time (min)');
    %     ylabel('Value');
    % 
    %      % xlim([418,426])
    % 
    % 
    % end
end

if (false && length(ekf_dt_values) > 1)

    timeVec_first = 0:ekf_dt_values(1):(length(prediction.first(1,:)) - 1) * ekf_dt_values(1);
    timeVec_second = 0:ekf_dt_values(2):(length(prediction.second(1,:)) - 1) * ekf_dt_values(2);

    if length(ekf_dt_values) == 3
        timeVec_third = 0:ekf_dt_values(3):(length(prediction.third(1,:)) - 1) * ekf_dt_values(3);
    end

    timeVec_BGs = 0:cgm_dt:(length(tools.BGs) - 1) * cgm_dt;

    % Plotting the predictions and BGs
    figure; % Create a new figure
    hold on; % Hold on to plot multiple lines

    % Plot the predictions
    what = 6;

    plot(timeVec_first, prediction.first(what,:)/params.VG,'b-','DisplayName', sprintf('Prediction dt=%.2f', ekf_dt_values(1)));
    plot(timeVec_second, prediction.second(what,:)/params.VG,'r-', 'DisplayName', sprintf('Prediction dt=%.2f', ekf_dt_values(2)));
    if length(ekf_dt_values) == 3
        plot(timeVec_third, prediction.third(what,:)/params.VG,'g-', 'DisplayName', sprintf('Prediction dt=%.2f', ekf_dt_values(3)));
    end

    % Plot the BGs
    plot(timeVec_BGs, tools.BGs, 'k', 'DisplayName', 'BGs');

    % Customize the plot
    xlabel('Time (min)');
    ylabel('Glucose Level');
    title('Glucose Predictions and BGs');
    legend('show'); % Show legend
    hold off; % Release the hold on the current figure

end

if false
    totalMinutes = length(model_predictions) * ekf_dt;
    time_vec = 0:ekf_dt:totalMinutes-ekf_dt;

    figure;

    subplot(7,2,1);
    plot(time_vec, model_predictions(1,:) + model_predictions(2,:), 'b-');
    xlabel('Time [min]');
    ylabel('Qsto [mg]');
    subplot(7,2,3);
    plot(time_vec, model_predictions(3,:), 'b-');
    xlabel('Time [min]');
    ylabel('Qgut [mg]');


    subplot(7,2,5);
    plot(time_vec, model_predictions(4,:), 'r-');
    xlabel('Time [min]');
    ylabel('Gp [mg/kg]');
    subplot(7,2,7);
    plot(time_vec, model_predictions(5,:), 'r-');
    xlabel('Time [min]');
    ylabel('Gt [mg/kg]');
    subplot(7,2,9);
    plot(time_vec, model_predictions(6,:), 'r-');
    xlabel('Time [min]');
    ylabel('Gsc [mg/kg]');

    subplot(7,2,11);
    plot(time_vec, model_predictions(7,:), 'g-');
    xlabel('Time [min]');
    ylabel('Il [pmol/kg]');
    subplot(7,2,13);
    plot(time_vec, model_predictions(8,:), 'g-');
    xlabel('Time [min]');
    ylabel('Ip [pmol/kg]');
    subplot(7,2,2);
    plot(time_vec, model_predictions(9,:), 'g-');
    xlabel('Time [min]');
    ylabel('I1 [pmol/l]');
    subplot(7,2,4);
    plot(time_vec, model_predictions(10,:), 'g-');
    xlabel('Time [min]');
    ylabel('Id [pmol/l]');
    subplot(7,2,6);
    plot(time_vec, model_predictions(11,:), 'g-');
    xlabel('Time [min]');
    ylabel('X [pmol/l]');

    subplot(7,2,8);
    plot(time_vec, model_predictions(12,:), 'k-');
    xlabel('Time [min]');
    ylabel('Isc1 [pmol/kg]');
    subplot(7,2,10);
    plot(time_vec, model_predictions(13,:), 'k-');
    xlabel('Time [min]');
    ylabel('Isc2 [pmol/kg]');

    subplot(7,2,12);
    plot(time_vec, tools.CHOs, 'm-');
    xlabel('Time [min]');
    ylabel('CHO[g]');
    subplot(7,2,14);
    plot(time_vec, tools.IIRs, 'm-');
    xlabel('Time [min]');
    ylabel('IIR[U]');

    set(gcf, 'Position', get(0, 'Screensize')); % Maximize figure window for clarity

end

%% Extra functions

function stop = stop_simulation(tools, CGM_read)

    stop = false;
    if CGM_read == length(tools.CGMs)-1
        stop = true;
    end
end

function [z_k, new_measurement_detected] = sample_measurement(t, tools)
    new_measurement_detected = false;
    z_k = [];
    if mod(t, 1) == 0
        CGM_read = t + 1;
        z_k = tools.CGMs(CGM_read);
        new_measurement_detected = true;
    end
end

function [u_k, new_input_detected] = sample_input(t, tools)
    new_input_detected = false;
    u_k.CHO = 0;
    u_k.IIR = 0;
    if mod(t, 1) == 0
        U_read = t + 1;
        u_k.CHO = tools.CHOs(U_read);
        u_k.IIR = tools.IIRs(U_read);
        new_input_detected = true;
    end
end