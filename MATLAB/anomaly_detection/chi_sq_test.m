anomalies = [];
alpha = 0.99;
for i = 1:length(residuals)
    residual = residuals(i);
    inn_cov = innovation_cov(i);
    anomalies(end+1) = chi_squared_test(residual, inn_cov, alpha);
end

function anomaly = chi_squared_test(residual, innovation_covariance, alpha)
    chi_sq_err = residual^2 / innovation_covariance;
    if chi_sq_err > chi2inv(alpha,1)
        anomaly = true;
    else
        anomaly = false;
    end
end