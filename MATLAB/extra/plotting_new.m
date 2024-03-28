if all_states

    figure;

    subplot(7,2,1);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(1,:) + EKF_state_tracking.mean(2,:), 'b-');
    ylabel('Qsto [mg]');
    subplot(7,2,3);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(3,:), 'b-');
    ylabel('Qgut [mg]');


    subplot(7,2,5);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(4,:), 'r-');
    ylabel('Gp [mg/kg]');
    subplot(7,2,7);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(5,:), 'r-');
    ylabel('Gt [mg/kg]');
    subplot(7,2,9);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(6,:), 'r-');
    ylabel('Gsc [mg/kg]');

    subplot(7,2,11);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(7,:), 'g-');
    ylabel('Il [pmol/kg]');
    subplot(7,2,13);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(8,:), 'g-');
    ylabel('Ip [pmol/kg]');
    subplot(7,2,2);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(9,:), 'g-');
    ylabel('I1 [pmol/l]');
    subplot(7,2,4);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(10,:), 'g-');
    ylabel('Id [pmol/l]');
    subplot(7,2,6);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(11,:), 'g-');
    ylabel('X [pmol/l]');

    subplot(7,2,8);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(12,:), 'k-');
    ylabel('Isc1 [pmol/kg]');
    subplot(7,2,10);
    plot(EKF_state_tracking.time, EKF_state_tracking.mean(13,:), 'k-');
    ylabel('Isc2 [pmol/kg]');

    subplot(7,2,12);
    stem(patientData.Meal.time, patientData.Meal.values, 'm-', 'Marker', '.');
    ylabel('CHO[g]');
    subplot(7,2,14);
    stem(patientData.IIR.time, patientData.IIR.values, 'm-', 'Marker', '.');
    ylabel('IIR[U]');

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