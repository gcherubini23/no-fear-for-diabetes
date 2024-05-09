close all
clear
clc

rng default;

run('config.m')
eps = 1e-5;


params_to_estimate = {'kp2','k1','k2','kp1','ki','ke1','kmax','kmin','kabs','kp3','Vmx'};
nvars = length(params_to_estimate);
prior = [0.0021, 0.065, 0.079, 2.70, 0.0079, 0.0005, 0.0558, 0.008, 0.0570, 0.009, 0.047];
cov_p = diag([0.1, 0.1, 0.1, 3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]);
ub = [1, 1, 1, 6, 1, 1, 1, 1, 1, 1, 1];  % Upper bounds
lb = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] + eps;  % Lower bounds


window_size = experiment_total_time;
high_uncertainty = 1000;
Q = diag([10,10,10,high_uncertainty,high_uncertainty,high_uncertainty,0,0,0,0,0,0,0]);    % TBD
% Q = eye(length(state_fields)) * high_uncertainty;
R = 20;  % TBD
ekf_dt = 1; % [min]
model = non_linear_model(tools);
lin_model = linearized_model(tools);
ekf = ekf(model, lin_model, tools, params, ekf_dt, Q, R);
ekf.dt = ekf_dt;

if use_basal_init_conditions
    [x0_, ymin1_] = tools.init_conditions(params);
else
    [x0_, ymin1_] = tools.set_init_conditions(patientData.CGM.values(1), params);
end

u0.CHO = 0;
u0.IIR = 0;
x0 = x0_;
ymin1 = ymin1_;
window = set_window(window_size, t_start, x0, ymin1, u0, t_end);


initial_params = prior;  
proposal_dist = @(p) mvnrnd(p, cov_p);
logpdf = @(p) log_posterior(p, cov_p, prior, params, ekf, patientData, window, params_to_estimate);

% Number of samples and burn-in period
nsamples = 600;  % Adjust based on your convergence diagnostics
burnin = 100;

% Set up the MCMC sampler using mhsample
[samples, accept_rate] = mhsample(initial_params, nsamples, ...
    'logpdf', logpdf, ...
    'proprnd', proposal_dist, ...
    'symmetric', true, ...
    'burnin', burnin);

% Output results
fprintf('Acceptance Rate: %.2f%%\n', accept_rate * 100);


function mapLoss = log_posterior(p, cov_p, prior, patient, ekf, patientData, window, params_to_estimate)
    patient = patient.set_params(params_to_estimate,p);
    t = window.t_start;
    x = window.x0;
    y = window.ymin1;
    u = window.u0;
    last_process_update = window.t_start;
    P = ekf.Q * 3;
    log_lik = 0;
    while t <= window.t_end
        [uk, new_input_detected] = sample_input(t, patientData);
        [zk, new_measurement_detected] = sample_measurement(t, patientData);
        if new_input_detected || new_measurement_detected
            dt = convert_to_minutes(t - last_process_update);
            [xk, Pk, ykmin1, ~] = ekf.predict(x, y, u, P, dt, patient);
            x = xk;
            P = Pk;
            y = ykmin1;
            u = uk;
            if new_measurement_detected
                % ekf = ekf.update_sensor_cov(zk);
                r = zk - ekf.H(6) * x.Gpd;
                log_lik = log_lik + r * r / ekf.R;
            end
            last_process_update = t;
        end
        t = next_step(t, patientData);
    end

    regularization_term = (p - prior) * cov_p * (p - prior)';
    
    alpha = 1;
    mapLoss = log_lik + alpha * regularization_term;
end
