clc
close all
clear all


h = figure;

font_size = 10;
figureWidth = 5.7;
figureHeight = 4;
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, figureWidth, figureHeight]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figureWidth, figureHeight]);
set(gcf, 'PaperPositionMode', 'auto');


hold on;

% Define colors
zone_colors = [[53,255,22];
              [39,128,0];
              [240,255,6];
              [250,168,3];
              [255,0,0]]/255;

% Zone A
fill([0 70 70 400 400 1000/3 175/3 0], [0 0 56 320 400 400 70 70], zone_colors(1,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Zone B
fill([70 130 180 240 240 400 400 70], [0 0 70 70 180 180 320 56], zone_colors(2,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([70 1000/3 290 70], [84 400 400 180], zone_colors(2,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% % Zone C
fill([130 180 180], [0 0 70], zone_colors(3,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([70 290 70], [180 400 400], zone_colors(3,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% 
% % Zone D
fill([240 400 400 240], [70 70 180 180], zone_colors(4,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([0 175/3 70 70 0], [70 70 84 180 180], zone_colors(4,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% 
% % Zone E
fill([0 70 70 0], [180 180 400 400], zone_colors(5,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([180 400 400 180], [0 0 70 70], zone_colors(5,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plotting lines
plot([0 400], [0 400], 'k:'); % Theoretical 45º regression line
plot([0 175/3], [70 70], 'k-');
plot([175/3 400/1.2], [70 400], 'k-'); % Replace 320 with 400/1.2 because 100*(400 - 400/1.2)/(400/1.2) = 20% error
plot([70 70], [84 400], 'k-');
plot([0 70], [180 180], 'k-');
plot([70 290], [180 400], 'k-'); % Corrected upper B-C boundary
plot([70 70], [0 56], 'k-'); % Replace 175.3 with 56 because 100*abs(56-70)/70) = 20% error
plot([70 400], [56 320], 'k-');
plot([180 180], [0 70], 'k-');
plot([180 400], [70 70], 'k-');
plot([240 240], [70 180], 'k-');
plot([240 400], [180 180], 'k-');
plot([130 180], [0 70], 'k-'); % Lower B-C boundary slope OK

% Adding labels
text(30, 20, 'A', 'FontSize', 14);
text(30, 150, 'D', 'FontSize', 14);
text(30, 380, 'E', 'FontSize', 14);
text(150, 380, 'C', 'FontSize', 14);
text(160, 20, 'C', 'FontSize', 14);
text(380, 20, 'E', 'FontSize', 14);
text(380, 120, 'D', 'FontSize', 14);
text(380, 260, 'B', 'FontSize', 14);
text(280, 380, 'B', 'FontSize', 14);

xlabel('Reference Concentration [mg/dL]');
ylabel('Predicted Concentration [mg/dL]');
title('Clarke''s Error Grid');
set(gca, 'XLim', [0 400]);
set(gca, 'YLim', [0 400]);
axis square;
set(h, 'color', 'white');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', font_size);
