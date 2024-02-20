close
clear
clc

state_fields = {'Qsto1','Qsto2','Qgut','Gp','Gt','Gsc','Il','Ip','I1','Id','X','Isc1','Isc2'};
extra_state_fields = {'insulin_to_infuse','last_IIR','CHO_to_eat','D','lastQsto','is_eating'}; 
input_fields = {'CHO', 'IIR'};
true_input_fields = {'CHO_consumed_rate','IIR'};

filename = "/Users/giovannicherubini/Desktop/Thesis/Code/no-fear-for-diabetes/data/1minsample/adult#001_3.csv";

Q = eye(numel(state_fields)) * 0.08;    % TBD
R = 0.00;  % TBD

% ekf_dt_values = [0.01, 0.1, 1];
ekf_dt_values = [1, 1];
ekf_dt = ekf_dt_values(1);
cgm_dt = 1;

tools = utils(filename, state_fields, extra_state_fields, input_fields, true_input_fields);
model = non_linear_model(tools);
lin_model = linearized_model(tools);
ekf = ekf(model, lin_model, tools, ekf_dt, Q, R);

basal = tools.IIRs(1);
params = patient_01(basal);

[x0, y0] = tools.init_conditions(params);
P = zeros(numel(state_fields),numel(state_fields));   % TBD

% Loop over the values of ekf_dt
prediction.first = [];
prediction.second = [];
prediction.third = [];
i = 1;
for ekf_dt = ekf_dt_values
    
    ekf.dt = ekf_dt;
    model_predictions = [];
    EKF_state_tracking.mean = [];
    EKF_state_tracking.variance = [];
    u_history = [];
    v_history = [];
    y_history.CHO_to_eat = [];
    y_history.D = [];
    y_history.last_Q_sto = [];
    y_history.is_eating = [];
    y_history.last_IIR = [];
    y_history.insulin_to_infuse = [];
    
    t = 0;
    CGM_read = 0;
    x = x0;
    y = y0;
    while ~stop_simulation(tools, CGM_read)
        new_measurement = false;
    
    
        % Get control input and eventually measurement
        if ekf_dt <= 1
            if t == 1 / ekf_dt || t == 0
                check = true;
                t = 0;
                increment = CGM_read + 1;
            else
                check = false;
            end
        else
            if mod(t,ekf_dt) == 0
                increment = t+1;
                check = true;
            else
                check = false;
            end
        end
    
        if check
            CGM_read = increment
            z = tools.CGMs(CGM_read);
            u.CHO = tools.CHOs(CGM_read);
            u.IIR = tools.IIRs(CGM_read);
            new_measurement = true;
        else
            u.CHO = 0;
            u.IIR = 0;
        end
        
        %% EKF
        % % Preprocess
        % [y_current, v] = model.preprocess(x, y, u, params, ekf_dt);
        % y = y_current;
    
        % % EKF
        % % Process update
        % [xp, Pp] = ekf.process_update(x,y,v,P,params);
        % 
        % % If there is measurement, do measurement update
        % if new_measurement
        %     [xm, Pm] = ekf.measurement_update(xp, Pp, z);
        %     x_current = xm;
        %     P_current = Pm;
        % else
        %     x_current = xp;
        %     P_current = Pp;
        % end
    
        %% Model
        % [y_true_current, v_true] = model.preprocess(x_true, y_true, u, params, ekf_dt);
        % y_true = y_true_current;
        % 
        % x_true_new = ekf.process_update(x_true, y_true, v_true, P, params);

        if i == 1
            what = true;
        else
            what = false;
        end

        [x_current, P_current, y_current, v] = ekf.process_update(x, y, u, P, params, what);




       
        %% Store results
        % EKF_state_tracking.mean(end+1) = x_current.Gp / params.VG;
        % EKF_state_tracking.variance(end+1) = 1 / params.VG * sqrt(P_current(4,4));
    
        model_predictions(:,end+1) = tools.convert_to_vector(x_current);
    
        u_history(end+1) = u.IIR;
        v_history(end+1) = v.IIR_dt;
        y_history.CHO_to_eat(end+1) = y_current.CHO_to_eat;
        y_history.D(end+1) = y_current.D;
        y_history.last_Q_sto(end+1) = y_current.lastQsto;
        y_history.is_eating(end+1) = y_current.is_eating;
        y_history.insulin_to_infuse(end+1) = y_current.insulin_to_infuse;
        y_history.last_IIR(end+1) = y_current.last_IIR;
        
        % Update
        x = x_current;
        y = y_current;
        % P = P_current;
        t = t + 1;
    end
    
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

    disp("Done");
end

%% Plot

if length(ekf_dt_values) == 1
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
        plot(timeVecBGs, tools.BGs(1:numPointsBGs), 'r-', 'DisplayName', 'BGs');
        hold on
        plot(timeVecEKF_Model, model_predictions(4,1:numPointsEKFModel)/params.VG, 'g-', 'DisplayName', 'Model Predictions');
        hold off;
        legend show;
        grid on;
        % title('EKF Mean and Model Predictions');
        % title('Model Predictions');
        xlabel('Time (min)');
        ylabel('Value');
        
        % % Subplot 2: tools.BGs
        % subplot(2,1,2);
        % plot(timeVecBGs, tools.BGs(1:numPointsBGs), 'r-', 'DisplayName', 'BGs');
        % legend show;
        % grid on;
        % title('BGs');
        % xlabel('Time (min)');
        % ylabel('BG Value');
        
        % % Subplot 3: u_history and v_history
        % subplot(3,1,3);
        % % plot(timeVecUV, u_history(1:numPointsUV), 'c-', 'DisplayName', 'u history');
        % % hold on;
        % plot(timeVecUV, y_history.last_Q_sto(1:numPointsUV), 'm-', 'DisplayName', 'Last Qsto');
        % % hold off;
        % legend show;
        % grid on;
        % title('Control Input History');
        % xlabel('Time (min)');
        % ylabel('Input Value');
        
        % Adjust figure properties for better visibility
        set(gcf, 'Position', get(0, 'Screensize')); % Maximize figure window for clarity
    end
    
    
    if (false)
        % Create figure
        figure;
        
        % Subplot 1
        subplot(3,1,1);
        plot(timeVecEKF_Model, model_predictions(10,1:numPointsEKFModel), 'g-', 'DisplayName', 'Id');
        legend show;
        grid on;
        xlabel('Time (min)');
        ylabel('Value');
    
        % Subplot 2
        subplot(3,1,2);
        plot(timeVecEKF_Model, model_predictions(9,1:numPointsEKFModel), 'g-', 'DisplayName', 'I1');
        hold on
        plot(timeVecEKF_Model, model_predictions(8,1:numPointsEKFModel) / params.VI, 'r-', 'DisplayName', 'It');
        hold off
        legend show;
        grid on;
        xlabel('Time (min)');
        ylabel('Value');
    
        % Subplot 3
        subplot(3,1,3);
        plot(timeVecEKF_Model, model_predictions(11,1:numPointsEKFModel), 'g-', 'DisplayName', 'X');
        legend show;
        grid on;
        xlabel('Time (min)');
        ylabel('Value');
    
    end
    
    
    if (false)
        % Create figure
        figure;
        
        % Subplot 1
        subplot(3,1,1);
        plot(timeVecEKF_Model, model_predictions(12,1:numPointsEKFModel)-params.Isc1ss, 'g-', 'DisplayName', 'Isc1');
        legend show;
        grid on;
        xlabel('Time (min)');
        ylabel('Value');
    
        % Subplot 2
        subplot(3,1,2);
        plot(timeVecEKF_Model, params.ka1 * model_predictions(12,1:numPointsEKFModel) + params.ka2 * model_predictions(13,1:numPointsEKFModel), 'g-', 'DisplayName', 'Rit');
        legend show;
        grid on;
        xlabel('Time (min)');
        ylabel('Value');
    
        % Subplot 3
        subplot(3,1,3);
        plot(timeVecEKF_Model, u_history(1:numPointsEKFModel)* 6000 / params.BW, 'g-', 'DisplayName', 'u');
        legend show;
        grid on;
        xlabel('Time (min)');
        ylabel('Value');
    
    
    end
    
    
    if (true)
        
        % Create figure
        figure;
        
        % Subplot 1
        subplot(2,1,1);
        plot(timeVecEKF_Model, y_history.insulin_to_infuse(1:numPointsUV), 'g-', 'DisplayName', 'Insulin to infuse');
        hold on;
        plot(timeVecEKF_Model, u_history(1:numPointsUV), 'b-', 'DisplayName', 'u');
        hold on
        plot(timeVecEKF_Model, v_history(1:numPointsUV), 'r-', 'DisplayName', 'v');
        hold off;
        legend show;
        grid on;
        xlabel('Time (min)');
        ylabel('Value');
    
        % Subplot 2
        subplot(2,1,2);
        plot(timeVecEKF_Model, y_history.last_IIR(1:numPointsUV), 'm-', 'DisplayName', 'last_IIR');
        legend show;
        grid on;
        xlabel('Time (min)');
        ylabel('Value');
    
    
    end
end

if (true && length(ekf_dt_values) > 1)

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
    what = 4;

    plot(timeVec_first, prediction.first(what,:)/params.VG, 'DisplayName', sprintf('Prediction dt=%.2f , euler', ekf_dt_values(1)));
    plot(timeVec_second, prediction.second(what,:)/params.VG, 'DisplayName', sprintf('Prediction dt=%.2f , rk4', ekf_dt_values(2)));
    if length(ekf_dt_values) == 3
        plot(timeVec_third, prediction.third(what,:)/params.VG, 'DisplayName', sprintf('Prediction dt=%.2f', ekf_dt_values(3)));
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