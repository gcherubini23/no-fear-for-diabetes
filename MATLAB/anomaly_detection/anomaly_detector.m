classdef anomaly_detector
    properties
        alpha;
        time = datetime([], 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss');
        anomalies = [];
    end

    methods
        function obj = anomaly_detector(alpha)
            obj.alpha = alpha;
        end

        function anomaly_detected = chi_squared_test(obj, residual, innovation_covariance)
            tau = chi2inv(1-obj.alpha,1);
            err = residual^2 / innovation_covariance;
            if err > tau
                anomaly_detected = true;
            else
                anomaly_detected = false;
            end
        end

        function [anomaly_detected, obj] = cusum_test(obj)

        end
    end
end