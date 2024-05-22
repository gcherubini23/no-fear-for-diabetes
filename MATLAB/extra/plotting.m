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
    
    if show_confidence_interval
        name = 'EKF Track \pm 2\sigma';
    else
        name = 'EKF Track';
    end

    plot(EKF_state_tracking.time, ekf.H * EKF_state_tracking.mean, '-', 'DisplayName', name, 'LineWidth', 1, 'Color', 'm');
    hold on
    if ~use_true_patient
        plot(tools.Time, tools.BGs, 'b-', 'DisplayName', 'Ground truth', 'LineWidth', 1);
        hold on
    end
    plot(patientData.CGM.time, patientData.CGM.values, 'o', 'DisplayName', 'CGM', 'Color', 'cyan', 'MarkerSize', 4)
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
        plot(anomaly_detector.time, anomaly_detector.anomalies, 'o','DisplayName', 'Anomalies', 'Color', [1 0 0], 'MarkerSize', 4)
        hold on
        if simulate_anomalies
            plot(true_CGM.time, true_CGM.values, 'o', 'Color', "#C0C0C0", 'MarkerSize', 4, 'HandleVisibility', 'off')
            hold on
        end
    end


    if show_future_predictions
        if show_pred_improvement
            for i = 1:length(trajectories)
                traj = trajectories{i};
                plot(traj.time, traj.values/params.VG, '-', 'LineWidth', 1, 'Color', "#C0C0C0", 'HandleVisibility', 'off')
                hold on
            end
            
        else
            if show_confidence_interval
                name = 'EKF Pred \pm 2\sigma';
            else
                name = 'EKF Pred';
            end
            plot(future_predictions.time, ekf.H * future_predictions.values, 'x', 'DisplayName', name, 'LineWidth', 1, 'Color', "#77AC30")
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
    
    legend show;
    xlim([t_start t_end]);
    ylim([-0.5, 400]);
    set(gcf, 'Position', get(0, 'Screensize'));
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