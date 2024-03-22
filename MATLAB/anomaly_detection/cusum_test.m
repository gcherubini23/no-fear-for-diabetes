anomalies = [];
std_residuals = [];
Ss = [];
S = 0;
alpha = 0.01;
tau = 1/(alpha);

b_tilda = 8;
b = 1 * b_tilda;

time = length(residuals);

for i = 1:length(residuals)

    % std_residual = abs(residuals(i) / innovation_cov(i)^(1/2));
    std_residual = residuals(i)^2 / innovation_cov(i);
    
    increment = b * (std_residual - b * 0.5);

    S = max(0, S + increment);
    Ss(end+1) = S;
    if S <= tau
        anomaly = false;
    else
        anomaly = true;
        i = i
        % S = 0;
    end
    anomalies(end+1) = anomaly;
end
