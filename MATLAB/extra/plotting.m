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
    ylabel('Gsc [mg/kg]');
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

    plot(EKF_state_tracking.time, EKF_state_tracking.mean(6,:)/params.VG, 'm-', 'DisplayName', 'EKF', 'LineWidth', 1);
    hold on
    if ~use_true_patient
        plot(tools.Time, tools.BGs, 'b-', 'DisplayName', 'Ground truth', 'LineWidth', 1);
        hold on
    end
    plot(patientData.CGM.time, patientData.CGM.values, 'o', 'DisplayName', 'CGM', 'Color', 'cyan', 'MarkerSize', 4)
    hold on

    if show_confidence_interval
        gamma = 1.96;
        sigma = sqrt(EKF_state_tracking.variance);
        upper_bound = (EKF_state_tracking.mean(6,:) + gamma * sigma)/params.VG;
        lower_bound = (EKF_state_tracking.mean(6,:) - gamma * sigma)/params.VG;
        x = transpose(EKF_state_tracking.time);
        x2 = [x, fliplr(x)];
        inBetween = [lower_bound, fliplr(upper_bound)];
        fill(x2, inBetween, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor','m','DisplayName','95% CI');
        hold on
    end
    
    if plot_anomalies
        plot(anomaly_detector.time, anomaly_detector.anomalies, 'o','DisplayName', 'Anomalies', 'Color', [1 0 0], 'MarkerSize', 4)
        hold on 
    end
    
    legend show;
    xlim([t_start t_end]);
    set(gcf, 'Position', get(0, 'Screensize'));
end

if plot_complete_history

    figure

    subplot(2,1,1)
    stem([ekf.t_history], [ekf.u_history.IIR], 'square', 'Color', 'g', 'DisplayName', 'u')
    hold on
    stem([ekf.t_history], [ekf.v_history.IIR_dt], 'x', 'Color', 'b', 'DisplayName', 'v (IIR dt)');
    hold on
    stem([ekf.t_history], [ekf.y_history.insulin_to_infuse], '.', 'Color', 'r', 'DisplayName', 'Insulin to infuse');
    hold off
    grid on
    ylabel('IIR')
    xlim([t_start t_end]);
    legend show

    subplot(2,1,2)
    stem([ekf.t_history], [ekf.u_history.CHO], 'square', 'Color', 'g', 'DisplayName', 'u');
    hold on
    stem([ekf.t_history], [ekf.v_history.CHO_consumed_rate], 'x', 'Color', 'b', 'DisplayName', 'v (CHO consumed rate)');
    hold on
    stem([ekf.t_history], [ekf.y_history.CHO_to_eat], '.', 'Color', 'r', 'DisplayName', 'CHO to eat');
    hold on
    plot([ekf.t_history], [ekf.y_history.D], 'Color', 'm', 'DisplayName', 'CHO eaten');
    grid on
    ylabel('CHO')
    xlim([t_start t_end]);
    legend show


    set(gcf, 'Position', get(0, 'Screensize'));

end