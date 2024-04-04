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

    plot(EKF_state_tracking.time, EKF_state_tracking.mean(6,:)/params.VG, 'm-', 'DisplayName', 'EKF');
    hold on
    if ~use_true_patient
        plot(tools.Time, tools.BGs, 'b-', 'DisplayName', 'Ground truth');
        hold on
    end
    plot(patientData.CGM.time, patientData.CGM.values, 'o', 'DisplayName', 'CGM', 'Color', 'cyan', 'MarkerSize', 4)
    hold on
    
    if plot_anomalies        
        plot(anomaly_detector.time, anomaly_detector.anomalies, 'o','DisplayName', 'Anomalies', 'Color', [1 0 0], 'MarkerSize', 4)
        hold on 
    end
    
    legend show;
    set(gcf, 'Position', get(0, 'Screensize'));
end

if plot_real_data && use_true_patient

    figure;

    subplot(3,1,1)
    plot(patientData.CGM.time, patientData.CGM.values, '-o', 'DisplayName', 'CGM', 'Color', 'cyan', 'MarkerSize', 4);
    xlabel('Time');
    xlim([t_start t_end]);
    legend show
    grid on

    subplot(3,1,2)
    % stem(patientData.Insulin.time, patientData.Insulin.values, 'o', 'Color', 'r', 'DisplayName', 'Insulin');
    stem(patientData.IIR.time, patientData.IIR.values, 'o', 'Color', 'r', 'DisplayName', 'Insulin');
    hold on
    stem(patientData.Meal.time, patientData.Meal.values, 'square', 'Color', 'g', 'DisplayName', 'Meal');
    xlabel('Time');
    xlim([t_start t_end]);
    legend show
    grid on

    subplot(3,1,3)
    plot(patientData.System.time, patientData.System.DiAsState, '.', 'Color', 'm', 'DisplayName', 'Operational mode');
    xlabel('Time');
    xlim([t_start t_end]);
    legend show
    grid on
    
    set(gcf, 'Position', get(0, 'Screensize'));
    
end

if plot_complete_history

    figure

    subplot(3,1,1)
    stem([ekf.t_history], [ekf.u_history.IIR], 'square', 'Color', 'g', 'DisplayName', 'u')
    hold on
    stem([ekf.t_history], [ekf.v_history.IIR_dt], 'o', 'Color', 'b', 'DisplayName', 'v (IIR dt)');
    hold on
    stem([ekf.t_history], [ekf.y_history.insulin_to_infuse], '.', 'Color', 'r', 'DisplayName', 'Insulin to infuse');
    hold off
    grid on
    ylabel('IIR')
    xlim([t_start t_end]);
    legend show

    subplot(3,1,2)
    stem([ekf.t_history], [ekf.y_history.last_IIR], '.', 'Color', 'r', 'DisplayName', 'last IIR')
    ylabel('IIR')
    xlim([t_start t_end]);
    grid on
    legend show

    subplot(3,1,3)
    stem([ekf.t_history], [ekf.u_history.CHO], 'square', 'Color', 'g', 'DisplayName', 'u');
    hold on
    stem([ekf.t_history], [ekf.v_history.CHO_consumed_rate], 'o', 'Color', 'b', 'DisplayName', 'v (CHO consumed rate)');
    hold on
    stem([ekf.t_history], [ekf.y_history.CHO_to_eat], '.', 'Color', 'r', 'DisplayName', 'CHO to eat');
    hold off
    grid on
    ylabel('CHO')
    xlim([t_start t_end]);
    legend show

    set(gcf, 'Position', get(0, 'Screensize'));

end