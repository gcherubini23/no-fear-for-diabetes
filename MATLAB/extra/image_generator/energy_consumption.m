clc
clear all
close all


% Define the file paths (adjust these as necessary)
file1 = "/Users/giovannicherubini/Desktop/Thesis/Code/data/energy/am.csv";
file2 = "/Users/giovannicherubini/Desktop/Thesis/Code/data/energy/lpm3.csv";

% Read the data from the first CSV file
data1 = readtable(file1);
time1 = data1.Time_ns_ * 1e-9 / 60; % Convert ns to min
current1 = data1.Current_nA_ * 1e-9; % Convert nA to A
voltage1 = data1.Voltage_mV_ * 1e-3; % Convert mV to V
energy1 = data1.Energy_uJ_ * 1e-6; % Convert µJ to J

% Read the data from the second CSV file
data2 = readtable(file2);
time2 = data2.Time_ns_ * 1e-9 / 60; % Convert ns to min
current2 = data2.Current_nA_ * 1e-9; % Convert nA to A
voltage2 = data2.Voltage_mV_ * 1e-3; % Convert mV to V
energy2 = data2.Energy_uJ_ * 1e-6; % Convert µJ to J


font_size = 10;
figureWidth = 5.7;
figureHeight = 4;
lw = 0.5;
ms = 2;
labels = {'Active Mode', 'Active Mode + LPM3'};

% Plot comparing current of the two files
figure;
images = [];
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, figureWidth, figureHeight]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figureWidth, figureHeight]);
set(gcf, 'PaperPositionMode', 'auto');
images(end+1) = plot(time1, current1, 'b-', 'LineWidth', 2*lw);
hold on;
images(end+1) = plot(time2, current2, 'r-', 'LineWidth', 2*lw);
xlabel('Time [min]');
ylabel('Current [A]');
title('Current Consumption');
lgd = legend(images, labels);
lgd.Location = 'northoutside';
lgd.Orientation = 'horizontal';
lgd.Box = 'off';
xlim([0.01 5])
ylim([1.65 4.5]*1e-4)
% axis square;
grid on;
set(findall(gcf, '-property', 'FontSize'), 'FontSize', font_size);

% % Plot comparing voltage of the two files
% figure;
% plot(time1, voltage1, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Voltage AM');
% hold on;
% plot(time2, voltage2, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Voltage LPM3');
% xlabel('Time [min]');
% ylabel('Voltage [V]');
% title('Voltage Comparison');
% legend;
% xlim([0.01 5])
% set(gca, 'FontSize', 14)
% axis square;
% grid on;

% Plot comparing energy consumption of the two files
figure;
images = [];
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, figureWidth, figureHeight]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figureWidth, figureHeight]);
set(gcf, 'PaperPositionMode', 'auto');
images(end+1) = plot(time1, energy1, 'b-', 'LineWidth', 2*lw);
hold on;
images(end+1) = plot(time2, energy2, 'r-', 'LineWidth', 2*lw);
xlabel('Time [min]');
ylabel('Energy [J]');
title('Energy Consumption');
lgd = legend(images, labels);
lgd.Location = 'northoutside';
lgd.Orientation = 'horizontal';
lgd.Box = 'off';
xlim([0.01 5])
% axis square;
grid on;
set(findall(gcf, '-property', 'FontSize'), 'FontSize', font_size);