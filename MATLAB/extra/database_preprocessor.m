if use_anderson
    folder_path = '/Users/giovannicherubini/Desktop/Thesis/Code/data/Anderson2016/processed_data';
    filenames = {'MonitorCGM_processed.csv', 'MonitorTotalBolus_processed.csv', 'MonitorMeal_processed.csv', 'MonitorSystem_processed.csv', 'CGM_processed.csv'};
    
    for filename = filenames
        file = strcat(folder_path, '/', string(filename));
    
        if strcmp(string(filename), 'MonitorCGM_processed.csv')
            databaseCGM = readtable(file);
        elseif strcmp(string(filename), 'MonitorTotalBolus_processed.csv')
            databaseInsulin = readtable(file);
        elseif strcmp(string(filename), 'MonitorMeal_processed.csv')
            databaseMeal = readtable(file);
        elseif strcmp(string(filename), 'MonitorSystem_processed.csv')
            databaseSystem = readtable(file);
        else
            database = readtable(file);
        end
    
    end
    
    
    idx = extract_idx(databaseCGM, patient_ID, start_day, days_to_examine);
    patientData.CGM.time = databaseCGM.Time(idx);
    patientData.CGM.values = databaseCGM.CGM(idx);
    
    idx = extract_idx(databaseInsulin, patient_ID, start_day, days_to_examine);
    patientData.IIR.time = databaseInsulin.Time(idx);
    patientData.IIR.values = databaseInsulin.DeliveredValue(idx);
    
    idx = extract_idx(databaseMeal, patient_ID, start_day, days_to_examine);
    patientData.Meal.time = databaseMeal.Time(idx);
    patientData.Meal.values = databaseMeal.MealSize(idx);
    
    idx = extract_idx(databaseSystem, patient_ID, start_day, days_to_examine);
    patientData.System.time = databaseSystem.Time(idx);
    patientData.System.DiAsState = databaseSystem.DiAsState(idx);

end

if use_tmoore
    patient_ID = -1;
    folder_path = '/Users/giovannicherubini/Desktop/Thesis/Code/data/TMoore/processed_data';
    filenames = {'BloodGlucose_processed.csv','DietaryCarbohydrates_processed.csv','InsulinDelivery_processed.csv'};
    
    for filename = filenames
        file = strcat(folder_path, '/', string(filename));
        if strcmp(string(filename), 'BloodGlucose_processed.csv')
            databaseCGM = readtable(file);
        elseif strcmp(string(filename), 'InsulinDelivery_processed.csv')
            databaseInsulin = readtable(file);
        elseif strcmp(string(filename), 'DietaryCarbohydrates_processed.csv')
            databaseMeal = readtable(file);
        end
    end
    
    idx = extract_idx(databaseCGM, patient_ID, start_day, days_to_examine);
    patientData.CGM.time = databaseCGM.Time(idx);
    patientData.CGM.values = databaseCGM.value(idx);
    
    idx = extract_idx(databaseInsulin, patient_ID, start_day, days_to_examine);
    patientData.IIR.time = databaseInsulin.Time(idx);
    patientData.IIR.values = databaseInsulin.value(idx);
    
    idx = extract_idx(databaseMeal, patient_ID, start_day, days_to_examine);
    patientData.Meal.time = databaseMeal.Time(idx);
    patientData.Meal.values = databaseMeal.value(idx);

end
    


if plot_true_database

    t_start = min([min(patientData.CGM.time), min(patientData.Meal.time), min(patientData.IIR.time)]);
    t_end = max([max(patientData.CGM.time), max(patientData.Meal.time), max(patientData.IIR.time)]);

    figure;

    subplot(2,1,1)
    plot(patientData.CGM.time, patientData.CGM.values, '-o', 'DisplayName', 'CGM', 'Color', 'cyan', 'MarkerSize', 4);
    xlabel('Time');
    xlim([t_start t_end]);
    legend show
    grid on

    subplot(2,1,2)
    % stem(patientData.Insulin.time, patientData.Insulin.values, 'o', 'Color', 'r', 'DisplayName', 'Insulin');
    stem(patientData.IIR.time, patientData.IIR.values, 'o', 'Color', 'r', 'DisplayName', 'Insulin');
    hold on
    stem(patientData.Meal.time, patientData.Meal.values, 'square', 'Color', 'g', 'DisplayName', 'Meal');
    xlabel('Time');
    xlim([t_start t_end]);
    legend show
    grid on

    % subplot(3,1,3)
    % plot(patientData.System.time, patientData.System.DiAsState, '.', 'Color', 'm', 'DisplayName', 'Operational mode');
    % xlabel('Time');
    % xlim([t_start t_end]);
    % legend show
    % grid on
    
    set(gcf, 'Position', get(0, 'Screensize'));
    
end

function idx = extract_idx(database, patient_ID, start_day, days_to_examine)
    
    if patient_ID > 0
        if ~strcmp(string(days_to_examine), 'all')
            idx = database.DeidentID == patient_ID & ...
                  database.Time >= start_day & ...
                  database.Time < start_day + days(days_to_examine);
        else
            idx = database.DeidentID == patient_ID;
        end
    else
        if ~strcmp(string(days_to_examine), 'all')
            idx = database.Time >= start_day & ...
                  database.Time < start_day + days(days_to_examine);
        else
            num_rows = size(database,1);
            idx = 1:num_rows;
        end
    end
end
