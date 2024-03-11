
if ekf_dt <= 1
    dt_EKF_Model = ekf_dt; % Time step for EKF mean and model predictions
    dt_BGs = cgm_dt; % Time step for BGs
    dt_UV = pump_dt; % Time step for u_history and v_history
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
    plot(timeVecBGs, tools.BGs(1:numPointsBGs), 'g-', 'DisplayName', 'True model');
    hold on
    plot(timeVecEKF_Model, model_predictions(6,1:numPointsEKFModel)/params.VG, 'b-', 'DisplayName', 'Nominal model + PSO optimization');
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


% if (true && length(ekf_dt_values) > 1)
% 
%     timeVec_first = 0:ekf_dt_values(1):(length(prediction.first(1,:)) - 1) * ekf_dt_values(1);
%     timeVec_second = 0:ekf_dt_values(2):(length(prediction.second(1,:)) - 1) * ekf_dt_values(2);
% 
%     if length(ekf_dt_values) == 3
%         timeVec_third = 0:ekf_dt_values(3):(length(prediction.third(1,:)) - 1) * ekf_dt_values(3);
%     end
% 
%     timeVec_BGs = 0:cgm_dt:(length(tools.BGs) - 1) * cgm_dt;
% 
%     % Plotting the predictions and BGs
%     figure; % Create a new figure
%     hold on; % Hold on to plot multiple lines
% 
%     % Plot the predictions
%     what = 6;
% 
%     plot(timeVec_first, prediction.first(what,:)/params.VG,'b-','DisplayName', sprintf('Prediction dt=%.2f', ekf_dt_values(1)));
%     plot(timeVec_second, prediction.second(what,:)/params.VG,'r-', 'DisplayName', sprintf('Prediction dt=%.2f', ekf_dt_values(2)));
%     if length(ekf_dt_values) == 3
%         plot(timeVec_third, prediction.third(what,:)/params.VG,'g-', 'DisplayName', sprintf('Prediction dt=%.2f', ekf_dt_values(3)));
%     end
% 
%     % Plot the BGs
%     plot(timeVec_BGs, tools.BGs, 'k', 'DisplayName', 'BGs');
% 
%     % Customize the plot
%     xlabel('Time (min)');
%     ylabel('Glucose Level');
%     title('Glucose Predictions and BGs');
%     legend('show'); % Show legend
%     hold off; % Release the hold on the current figure
% 
% end

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