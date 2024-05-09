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
    
if use_shanghai
    folder_path = '/Users/giovannicherubini/Desktop/Thesis/Code/data/ShanghaiDataset/processed_data';

    file = strcat(folder_path, '/', string(filename));

    database = readtable(file);
    
    idx = extract_idx(database, patient_ID, start_day, days_to_examine);
    patientData.CGM.time = database.Time(idx);
    patientData.CGM.values = database.CGM(idx);
    
    patientData.IIR.time = database.Time(idx);
    patientData.IIR.values = database.Insulin(idx);
    
    patientData.Meal.time = database.Time(idx);
    patientData.Meal.values = database.CHO(idx);
    
end


if plot_true_database

    t_start = min([min(patientData.CGM.time), min(patientData.Meal.time), min(patientData.IIR.time)]);
    t_end = max([max(patientData.CGM.time), max(patientData.Meal.time), max(patientData.IIR.time)]);

    figure;

    if use_anderson
        tot = 3;
    else
        tot = 2;
    end

    subplot(tot,1,1)
    plot(patientData.CGM.time, patientData.CGM.values, '-o', 'DisplayName', 'CGM', 'Color', 'cyan', 'MarkerSize', 4);
    xlabel('Time');
    xlim([t_start t_end]);
    legend show
    grid on

    subplot(tot,1,2)
    % stem(patientData.Insulin.time, patientData.Insulin.values, 'o', 'Color', 'r', 'DisplayName', 'Insulin');
    stem(patientData.IIR.time, patientData.IIR.values, 'o', 'Color', 'r', 'DisplayName', 'Insulin');
    hold on
    stem(patientData.Meal.time, patientData.Meal.values, 'square', 'Color', 'g', 'DisplayName', 'Meal');
    xlabel('Time');
    xlim([t_start t_end]);
    legend show
    grid on

    if use_anderson
        subplot(tot,1,tot)
        plot(patientData.System.time, patientData.System.DiAsState, '.', 'Color', 'm', 'DisplayName', 'Operational mode');
        xlabel('Time');
        xlim([t_start t_end]);
        legend show
        grid on
    end
    
    set(gcf, 'Position', get(0, 'Screensize'));
    
end

function idx = extract_idx(database, patient_ID, start_day, days_to_examine)
    
    if patient_ID > 0 && patient_ID < 700
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
