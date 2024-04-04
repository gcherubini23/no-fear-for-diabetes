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

%% Select days and patient
patient_ID = 11;
date = '14-Feb-2013 06:30:00';
start_day = datetime(date,'InputFormat', 'dd-MMM-yyyy HH:mm:ss');

days_to_examine = 2;
% days_to_examine = 'all';

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

function idx = extract_idx(database, patient_ID, start_day, days_to_examine)
    if ~strcmp(string(days_to_examine), 'all')
        idx = database.DeidentID == patient_ID & ...
              database.Time >= start_day & ...
              database.Time < start_day + days(days_to_examine);
    else
        idx = database.DeidentID == patient_ID;
    end
end
