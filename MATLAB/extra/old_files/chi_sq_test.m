anomalies = [];
alpha = 0.01;
tau = chi2inv(1-alpha,1);
for i = 1:length(residuals)
    residual = residuals(i);
    inn_cov = innovation_cov(i);
    anomaly = chi_squared_test(residual, inn_cov, tau);
    anomalies(end+1) = anomaly;
    if anomaly
        i = i
    end
end

function anomaly = chi_squared_test(residual, innovation_covariance, tau)
    chi_sq_err = residual^2 / innovation_covariance;
    if chi_sq_err > tau
        anomaly = true;
    else
        anomaly = false;
    end
end