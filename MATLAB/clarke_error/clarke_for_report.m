function [total, percentage] = clarke_for_report(y, yp)
    % CLARKE    Performs Clarke Error Grid Analysis
    %
    % The Clarke error grid approach is used to assess the clinical
    % significance of differences between the glucose measurement technique
    % under test and the venous blood glucose reference measurements. The
    % method uses a Cartesian diagram, in which the values predicted by the
    % technique under test are displayed on the y-axis, whereas the values
    % received from the reference method are displayed on the x-axis. The
    % diagonal represents the perfect agreement between the two, whereas the
    % points below and above the line indicate, respectively, overestimation
    % and underestimation of the actual values. Zone A (acceptable) represents
    % the glucose values that deviate from the reference values by ±20% or are
    % in the hypoglycemic range (<70 mg/dl), when the reference is also within
    % the hypoglycemic range. The values within this range are clinically exact
    % and are thus characterized by correct clinical treatment. Zone B (benign
    % errors) is located above and below zone A; this zone represents those
    % values that deviate from the reference values, which are incremented by
    % 20%. The values that fall within zones A and B are clinically acceptable,
    % whereas the values included in areas C-E are potentially dangerous, and
    % there is a possibility of making clinically significant mistakes. [1-4]
    %
    % SYNTAX:
    % 
    % [total, percentage] = clarke(y, yp)
    % 
    % INPUTS: 
    % y             Reference values (mg/dl) 
    % yp            Predicted/estimated values (mg/dl)
    % 
    % OUTPUTS: 
    % total         Total points per zone: 
    %               total(1) = zone A, 
    %               total(2) = zone B, and so on
    % percentage    Percentage of data which fell in certain region:
    %               percentage(1) = zone A, 
    %               percentage(2) = zone B, and so on.
    % 
    % EXAMPLE:      load example_data.mat 
    %               [tot, per] = clarke(y, yp)
    % 
    % © Edgar Guevara Codina 
    % codina@REMOVETHIScactus.iico.uaslp.mx 
    % File Version 1.2 
    % March 29 2013 

    % Error checking
    if nargin == 0
        error('clarke:Inputs', 'There are no inputs.')
    end
    if length(yp) ~= length(y)
        error('clarke:Inputs', 'Vectors y and yp must be the same length.')
    end
    if (max(y) > 400) || (max(yp) > 400) || (min(y) < 0) || (min(yp) < 0)
        error('clarke:Inputs', 'Vectors y and yp are not in the physiological range of glucose (<400mg/dl).')
    end
    
    % Print figure flag
    PRINT_FIGURE = true;
    
    % Determine data length
    n = length(y);
    
    % Plot Clarke's Error Grid
    h = figure;

    font_size = 9;
    ms = 1;
    figureWidth = 5.7/2;
    figureHeight = 5.7/2;
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


    plot(y,yp,'ko','MarkerSize',ms,'MarkerFaceColor','k','MarkerEdgeColor','k');
    xlabel('CGM measurement [mg/dL]');
    ylabel('Predicted Blood Glucose [mg/dL]');
    title('Clarke''s Error Grid');
    set(gca, 'XLim', [0 400]);
    set(gca, 'YLim', [0 400]);
    axis square;
    set(h, 'color', 'white');
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', font_size);

    
    total = zeros(5, 1); % Initializes output
    
    % Statistics
    for i = 1:n
        if (yp(i) <= 70 && y(i) <= 70) || (yp(i) <= 1.2 * y(i) && yp(i) >= 0.8 * y(i))
            total(1) = total(1) + 1; % Zone A
        else
            if ((y(i) >= 180) && (yp(i) <= 70)) || ((y(i) <= 70) && yp(i) >= 180)
                total(5) = total(5) + 1; % Zone E
            else
                if ((y(i) >= 70 && y(i) <= 290) && (yp(i) >= y(i) + 110)) || ((y(i) >= 130 && y(i) <= 180) && (yp(i) <= (7/5) * y(i) - 182))
                    total(3) = total(3) + 1; % Zone C
                else
                    if ((y(i) >= 240) && ((yp(i) >= 70) && (yp(i) <= 180))) || (y(i) <= 175/3 && (yp(i) <= 180) && (yp(i) >= 70)) || ((y(i) >= 175/3 && y(i) <= 70) && (yp(i) >= (6/5) * y(i)))
                        total(4) = total(4) + 1; % Zone D
                    else
                        total(2) = total(2) + 1; % Zone B
                    end % End of 4th if
                end % End of 3rd if
            end % End of 2nd if
        end % End of 1st if
    end % End of for loop
    
    percentage = (total./n) * 100;
end


