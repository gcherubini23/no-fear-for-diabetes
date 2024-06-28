clc
clear all
close all


% Sample data arrays (replace these with your actual data)
horizon = 15:15:180;  % Example horizon values [min]
rmse_tuned = [11.25, 16.22, 20.08, 23.15, 25.79, 28.22, 31.05, 34.58, 38.93, 44.13, 49.55, 54.82];  % Example RMSE values for Tuned model [mg/dL]
time_gain_tuned = [0, 5, 15, 25, 30, 40, 50, 55, 65, 70, 70, 75];  % Example Time Gain values for Tuned model [min]

time_gain_untuned = [0, 5, 10, 15, 20, 25, 30, 30, 35, 35, 35, 35];  % Example Time Gain values for Untuned model [min]
rmse_untuned = [13.91, 19.37, 24.26, 28.76, 32.87, 36.37, 39.36, 41.95, 44.10, 46.16, 48.22, 50.17];  % Example RMSE values for Untuned model [mg/dL]

% Set font size
fontSize = 14;

% Plot RMSE
figure;
plot(horizon, rmse_tuned, '-o', 'LineWidth', 1.5, 'DisplayName', 'Tuned');
hold on;
plot(horizon, rmse_untuned, '-x', 'LineWidth', 1.5, 'DisplayName', 'Untuned');
xlabel('Horizon [min]', 'FontSize', fontSize);
ylabel('RMSE [mg/dL]', 'FontSize', fontSize);
title('RMSE', 'FontSize', fontSize);
legend('Location', 'best', 'FontSize', fontSize);
grid on;
xlim([0 180])
axis square;
set(gca, 'FontSize', fontSize);

% Plot Time Gain
figure;
plot(horizon, time_gain_tuned, '-o', 'LineWidth', 1.5, 'DisplayName', 'Tuned');
hold on;
plot(horizon, time_gain_untuned, '-x', 'LineWidth', 1.5, 'DisplayName', 'Untuned');
xlabel('Horizon [min]', 'FontSize', fontSize);
ylabel('Time Gain [min]', 'FontSize', fontSize);
title('Time Gain', 'FontSize', fontSize);
legend('Location', 'best', 'FontSize', fontSize);
grid on;
xlim([0 180])
axis square;
set(gca, 'FontSize', fontSize);
