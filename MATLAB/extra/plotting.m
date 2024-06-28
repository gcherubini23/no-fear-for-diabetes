if all_states

    figure;

    subplot(7,2,1);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(1,:) + EKF_state_tracking.mean(2,:), 'b-');
    ylabel('Qsto [mg]');
    xlim([t_start t_end]);
    subplot(7,2,3);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(3,:), 'b-');
    ylabel('Qgut [mg]');
    xlim([t_start t_end]);


    subplot(7,2,5);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(4,:), 'r-');
    ylabel('Gp [mg/kg]');
    xlim([t_start t_end]);
    subplot(7,2,7);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(5,:), 'r-');
    ylabel('Gt [mg/kg]');
    xlim([t_start t_end]);
    subplot(7,2,9);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(6,:), 'r-');
    ylabel('Gsc [mg/dl]');
    xlim([t_start t_end]);

    subplot(7,2,11);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(7,:), 'g-');
    ylabel('Il [pmol/kg]');
    xlim([t_start t_end]);
    subplot(7,2,13);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(8,:), 'g-');
    ylabel('Ip [pmol/kg]');
    xlim([t_start t_end]);
    subplot(7,2,2);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(9,:), 'g-');
    ylabel('I1 [pmol/l]');
    xlim([t_start t_end]);
    subplot(7,2,4);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(10,:), 'g-');
    ylabel('Id [pmol/l]');
    xlim([t_start t_end]);
    subplot(7,2,6);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(11,:), 'g-');
    ylabel('X [pmol/l]');
    xlim([t_start t_end]);

    subplot(7,2,8);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(12,:), 'k-');
    ylabel('Isc1 [pmol/kg]');
    xlim([t_start t_end]);
    subplot(7,2,10);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(13,:), 'k-');
    ylabel('Isc2 [pmol/kg]');
    xlim([t_start t_end]);

    subplot(7,2,12);
    stem(patientData.Meal.time, patientData.Meal.values, 'm-', 'Marker', '.');
    ylabel('CHO[g]');
    xlim([t_start t_end]);
    subplot(7,2,14);
    stem(patientData.IIR.time, patientData.IIR.values, 'm-', 'Marker', '.');
    ylabel('IIR[U]');
    xlim([t_start t_end]);

    set(gcf, 'Position', get(0, 'Screensize')); % Maximize figure window for clarity
end

if only_Gpd

    figure
    ms = 6;
    
    if show_confidence_interval
        name = 'Tracking \pm 2\sigma';
    else
        name = 'Tracking';
    end

    plot(EKF_state_tracking.time, ekf.H * EKF_state_tracking.mean, '-', 'DisplayName', name, 'LineWidth', 2, 'Color', 'm');
    hold on
    if ~use_true_patient
        plot(tools.Time, tools.BGs, 'b-', 'DisplayName', 'Ground truth', 'LineWidth', 1);
        hold on
    end
    % plot(patientData.CGM.time, patientData.CGM.values, '-o', 'DisplayName', 'CGM', 'Color', 'cyan', 'MarkerSize', 8)
    plot(patientData.CGM.time, patientData.CGM.values, '-o', 'DisplayName', 'CGM', 'Color', 'cyan', 'LineWidth', 1, 'MarkerSize', ms)
    hold on

    if show_confidence_interval
        gamma = 2;
        sigma = sqrt(EKF_state_tracking.variance);
        upper_bound = ekf.H * EKF_state_tracking.mean + gamma * sigma;
        lower_bound = ekf.H * EKF_state_tracking.mean - gamma * sigma;
        x = transpose(EKF_state_tracking.time);
        x2 = [x, fliplr(x)];
        inBetween = [lower_bound, fliplr(upper_bound)];
        fill(x2, inBetween, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor','m','HandleVisibility', 'off');
        hold on
    end
    
    if plot_anomalies
        if simulate_anomalies
            plot(true_CGM.time, true_CGM.values, 'o', 'Color', "#C0C0C0", 'MarkerSize', ms, 'HandleVisibility', 'off')
            hold on
        end
        plot(anomaly_detector.time, anomaly_detector.anomalies, 'o','DisplayName', 'Anomalies', 'Color', [1 0 0], 'MarkerSize', ms, 'MarkerFaceColor', [1 0 0])
        hold on
    end


    if show_future_predictions
        if show_pred_improvement
            for i = 1:length(trajectories)
                traj = trajectories{i};
                plot(traj.time, traj.values, '-', 'LineWidth', 1, 'Color', "#C0C0C0", 'HandleVisibility', 'off')
                hold on
            end
            
        else
            if show_confidence_interval
                name = 'Prediction \pm 2\sigma';
            else
                name = 'Prediction';
            end
            plot(future_predictions.time, ekf.H * future_predictions.values, '-', 'DisplayName', name, 'LineWidth', 2, 'Color', "#77AC30")
            hold on
            if show_confidence_interval
                gamma = 2;
                sigma = sqrt(future_predictions.cov);
                upper_bound = ekf.H * future_predictions.values + gamma * sigma;
                lower_bound = ekf.H * future_predictions.values - gamma * sigma;
                x = transpose(future_predictions.time);
                x2 = [x, fliplr(x)];
                inBetween = [lower_bound, fliplr(upper_bound)];
                fill(x2, inBetween, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor',"#77AC30", 'HandleVisibility', 'off');
                hold on
            end
        end
    end


    
    legend_position = [1-0.18, 0.85];
    

    x_range = xlim;
    y_range = ylim;
    
    % Define the normalized position for the horizon segment
    x1 = x_range(1) + (x_range(2) - x_range(1)) * legend_position(1);
    x2 = x1 + minutes(horizon);
    y = y_range(1) + (y_range(2) - y_range(1)) * legend_position(2);
    
    % Add the segment
    line([x1 x2], [y y], 'Color', 'r', 'LineWidth', 2, 'HandleVisibility', 'off');
    
    % Add the text label for the horizon
    text(x2+minutes(100) , y, sprintf('%d minutes', horizon), 'Color', 'r', 'HorizontalAlignment', 'center','FontSize',14);


    
    legend show;
    xlim([t_start t_end]);
    ylim([-0.5, 350]);
    set(gcf, 'Position', get(0, 'Screensize'));
    xlabel('Time')
    ylabel('Plasma Glucose Concentration [mg/dL]')
    if do_measurment_update
        if do_chi_sq_test || do_cusum_test
            title('Glucose Estimator Anomaly Detection')
        elseif show_future_predictions
            title('Glucose Estimator Tracking and Prediction')
        else
            title('Glucose Estimator Tracking')
        end
    else
        title('Process Model Performance');
    end
    set(gca, 'FontSize', 14)
end

if plot_complete_history

    figure

    subplot(2,1,1)
    stem([ekf.t_history], [ekf.u_history(2,:)], 'square', 'Color', 'g', 'DisplayName', 'u')
    hold on
    stem([ekf.t_history], [ekf.v_history(2,:)], 'x', 'Color', 'b', 'DisplayName', 'v (IIR dt)');
    hold on
    stem([ekf.t_history], [ekf.y_history(1,:)], '.', 'Color', 'r', 'DisplayName', 'Insulin to infuse');
    hold off
    grid on
    ylabel('IIR')
    xlim([t_start t_end]);
    legend show

    subplot(2,1,2)
    stem([ekf.t_history], [ekf.u_history(1,:)], 'square', 'Color', 'g', 'DisplayName', 'u');
    hold on
    stem([ekf.t_history], [ekf.v_history(1,:)], 'x', 'Color', 'b', 'DisplayName', 'v (CHO consumed rate)');
    hold on
    stem([ekf.t_history], [ekf.y_history(3,:)], '.', 'Color', 'r', 'DisplayName', 'CHO to eat');
    hold on
    plot([ekf.t_history], [ekf.y_history(4,:)], 'Color', 'm', 'DisplayName', 'CHO eaten');
    grid on
    ylabel('CHO')
    xlim([t_start t_end]);
    legend show


    set(gcf, 'Position', get(0, 'Screensize'));

end

% 
% if false
% 
%     figure
%     plot(EKF_state_tracking.time, EKF_state_tracking.mean(8,:), 'g-');
%     ylabel('Ip [pmol/kg]');
%     xlim([t_start t_end]);
% 
% end

%%

if ~show_pred_improvement  
    [total, percentage] = colored_clarke(CGM_to_check,predictions_to_check, horizon);

    figure
    areas = {'A','B','C','D','E'};
    colors = [53,255,22;
              22,161,0;
              240,255,6;
              250,168,3;
              255,0,0]/255;

    b = bar(areas, percentage, 'facecolor', 'flat');
    b.CData = colors;
    xlabel('Areas');
    ylabel('Percentage');
    title('Clarke''s Error Grid Percentage')
    % annotation('textbox', [0.65, 0.795, 0.1, 0.1], 'String', 'Horizon: 180 min', 'FitBoxToText', 'on', 'BackgroundColor', 'white', 'FontSize', 14);
    set(gca, 'FontSize', 14);
    % set(gca,'XLim',[1, 4]);
    set(gca,'YLim',[0 100]);
    axis square

    % [0.75, 0.895, 0.1, 0.1]
end


