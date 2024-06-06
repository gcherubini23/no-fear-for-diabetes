clc
clear all
close all



% Define the file paths (adjust these as necessary)
file1 = "/Users/giovannicherubini/Desktop/Thesis/Code/data/energy/am.csv";
file2 = "/Users/giovannicherubini/Desktop/Thesis/Code/data/energy/lpm3.csv";

% Read the data from the first CSV file
data1 = readtable(file1);
time1 = data1.Time_ns_ * 1e-9; % Convert ns to s
current1 = data1.Current_nA_ * 1e-9; % Convert nA to A
voltage1 = data1.Voltage_mV_ * 1e-3; % Convert mV to V
energy1 = data1.Energy_uJ_ * 1e-6; % Convert µJ to J

% Read the data from the second CSV file
data2 = readtable(file2);
time2 = data2.Time_ns_ * 1e-9; % Convert ns to s
current2 = data2.Current_nA_ * 1e-9; % Convert nA to A
voltage2 = data2.Voltage_mV_ * 1e-3; % Convert mV to V
energy2 = data2.Energy_uJ_ * 1e-6; % Convert µJ to J

% Plot comparing current of the two files
figure;
plot(time1, current1, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Current AM');
hold on;
plot(time2, current2, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Current LPM3');
xlabel('Time (s)');
ylabel('Current (A)');
title('Current Comparison');
legend;
xlim([1 300])
set(gca, 'FontSize', 14)
grid on;

% Plot comparing voltage of the two files
figure;
plot(time1, voltage1, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Voltage AM');
hold on;
plot(time2, voltage2, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Voltage LPM3');
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Voltage Comparison');
legend;
xlim([1 300])
set(gca, 'FontSize', 14)
grid on;

% Plot comparing energy consumption of the two files
figure;
plot(time1, energy1, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Energy AM');
hold on;
plot(time2, energy2, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Energy LPM3');
xlabel('Time (s)');
ylabel('Energy (J)');
title('Energy Consumption Comparison');
legend;
xlim([1 300])
set(gca, 'FontSize', 14)
grid on;