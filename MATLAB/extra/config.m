% state_fields = {'Qsto1','Qsto2','Qgut','Gp','Gt','Gpd','Il','Ip','I1','Id','X','Isc1','Isc2'};
% extra_state_fields = {'insulin_to_infuse','last_IIR','CHO_to_eat','D','lastQsto','is_eating'}; 
% input_fields = {'CHO', 'IIR'};
% true_input_fields = {'CHO_consumed_rate','IIR_dt'};

indices;

use_true_patient = false;
use_tuned_model = false;

if use_true_patient
    use_anderson = true;
    use_tmoore = false;
    use_shanghai = false;

    if use_anderson
        patient_ID = 11;
        dailyBasal = 16.9;
        date = '11-Feb-2013 06:30:00';
        % date = '26-Jan-2013 06:30:00';
        days_to_examine = 2;
        % days_to_examine = 30;
        % days_to_examine = 'all';

        % patient_ID = 17;
        % % date = '20-Feb-2014 06:30:00';
        % date = '28-Jul-2013 12:18:22';
        % days_to_examine = 2;
        % % days_to_examine = 'all';
        % dailyBasal = 21.5;
    end

    if use_tmoore
        patient_ID = -1;
        dailyBasal = 10;
        date = '20-Jan-2022 00:00:00';

        days_to_examine = 2;
        % days_to_examine = 30;
        % days_to_examine = 'all';

    end

    if use_shanghai

        patient_ID = 1007;
        date = '27-Jul-2021 06:00:00';
        % date = '03-Aug-2021 09:00:00';
        days_to_examine = 2;
        % days_to_examine = 'all';
         
        % patient_ID = 1002;
        % date = '09-Sep-2021 16:00:00';
        % days_to_examine = 2;
        % % days_to_examine = 'all';


        % patient_ID = 1010;
        % date = '15-Sep-2021 06:30:00';
        % days_to_examine = 3;
        % % days_to_examine = 'all';


        if patient_ID == 1007
            filename = '1007_0_20210726_processed.csv';
        elseif patient_ID == 1002
            filename = '1002_2_20210909_processed.csv';
        elseif patient_ID == 1010
            filename = '1010_0_20210915_processed.csv';
        end

    end

    start_day = datetime(date,'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
end

use_basal_init_conditions = true;
do_measurment_update = true;

simulate_anomalies = true;
do_chi_sq_test = simulate_anomalies;
do_cusum_test = false;
loose_track_of_time = false;

do_plots = true;
if do_plots
    plot_true_database = false;
    all_states = true;
    only_Gpd = true;
    plot_anomalies = true;
    plot_complete_history = false;
    show_confidence_interval = true;
    show_future_predictions = false;
    show_pred_improvement = false;
    save_energy = false;
else
    plot_true_database = false;
    all_states = false;
    only_Gpd = false;
    plot_anomalies = false;
    plot_complete_history = false;
    show_confidence_interval = false;
    show_future_predictions = false;
    show_pred_improvement = false;
    save_energy = false;
end

%% Load data
disp('Loading dataset...')
if ~use_true_patient
    filename = "/Users/giovannicherubini/Desktop/Thesis/Code/data/1minsample/adult#001_5.csv";
    tools = utils(filename);
    patientData.CGM.values = tools.CGMs;
    patientData.CGM.time = tools.Time;
    patientData.Meal.values = tools.CHOs;
    patientData.Meal.time = tools.Time;
    patientData.IIR.values = tools.IIRs;
    patientData.IIR.time = tools.Time;
    patientData.BG.values = tools.BGs;
    patientData.BG.time = tools.Time;

    basal = tools.IIRs(1);

    use_true_model = false;
    if use_true_model
        params = patient_01(basal);
    else
        params = patient_00(basal);
    end
else
    
    tools = utils("none");
    run('database_preprocessor.m')

    if patient_ID == 17
        params = patient_17(dailyBasal);
    elseif patient_ID == 11
        params = patient_11(dailyBasal);
    elseif patient_ID == 1007
        params = patient_1007();
    elseif patient_ID == 1002
        params = patient_1002(); 
    elseif patient_ID == 1010
        params = patient_1010(); 
    end
end

disp('Dataset loaded')

t_start = min([min(patientData.CGM.time), min(patientData.Meal.time), min(patientData.IIR.time)]);
t_end = max([max(patientData.CGM.time), max(patientData.Meal.time), max(patientData.IIR.time)]);
experiment_total_time = t_end - t_start;