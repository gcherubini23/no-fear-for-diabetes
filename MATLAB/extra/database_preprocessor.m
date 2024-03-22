folder_path = '/Users/giovannicherubini/Desktop/Thesis/Code/data/Anderson2016/processed_data';
filenames = {'MonitorCGM_processed.csv', 'MonitorTotalBolus_processed.csv', 'MonitorMeal_processed.csv', 'MonitorSystem_processed.csv'};

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

%% Select days and patient and plot
plot_real_data_in_advance = true;
% patient_ID = 3;
% date = '25-Apr-2015';
% date = '24-Apr-2015';
patient_ID = 11;
date = '10-Feb-2013';
start_day = datetime(date,'InputFormat', 'dd-MMM-yyyy');
days_to_examine = 2;

idx = extract_idx(databaseCGM, patient_ID, start_day, days_to_examine);
patientData.CGM.time = databaseCGM.Time(idx);
patientData.CGM.values = databaseCGM.CGM(idx);

idx = extract_idx(databaseInsulin, patient_ID, start_day, days_to_examine);
patientData.Insulin.time = databaseInsulin.Time(idx);
patientData.Insulin.values = databaseInsulin.DeliveredValue(idx);

idx = extract_idx(databaseMeal, patient_ID, start_day, days_to_examine);
patientData.Meal.time = databaseMeal.Time(idx);
patientData.Meal.values = databaseMeal.MealSize(idx);

idx = extract_idx(databaseSystem, patient_ID, start_day, days_to_examine);
patientData.System.time = databaseSystem.Time(idx);
patientData.System.DiAsState = databaseSystem.DiAsState(idx);

close all
clc

if plot_real_data_in_advance

    figure;

    subplot(3,1,1)
    plot(patientData.CGM.time, patientData.CGM.values, '-o', 'DisplayName', 'CGM');
    xlabel('Time');
    legend show
    grid on

    subplot(3,1,2)
    stem(patientData.Insulin.time, patientData.Insulin.values, 'o', 'Color', 'r', 'DisplayName', 'Insulin');
    hold on
    stem(patientData.Meal.time, patientData.Meal.values, 'square', 'Color', 'g', 'DisplayName', 'Meal');
    xlabel('Time');
    legend show
    grid on

    subplot(3,1,3)
    plot(patientData.System.time, patientData.System.DiAsState, '.', 'Color', 'm', 'DisplayName', 'Operational mode');
    xlabel('Time');
    legend show
    grid on
    
    set(gcf, 'Position', get(0, 'Screensize'));

end

function idx = extract_idx(database, patient_ID, start_day, days_to_examine)
    idx = database.DeidentID == patient_ID & ...
          database.Time >= start_day & ...
          database.Time < start_day + days(days_to_examine);

    % idx = database.DeidentID == patient_ID;
end
