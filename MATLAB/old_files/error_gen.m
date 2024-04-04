add_bias = false;
turnoff_sensor = true;

t_start = 1;
t_end = length(tools.CGMs);

rng("default")

if add_bias
    t_bias_start = 1000;
    drift = 0;
    % Standard deviation of the Gaussian random walk step
    sigma = 1; % This can be changed based on how strong you want the drift to be
    
    % Iterate from the bias start time to the end and generate the drift
    for k = t_bias_start:t_end
        % Generate a random step from a Gaussian distribution
        step = sigma * randn;
        
        % Accumulate the drift
        drift = drift + step;
        
        % Add the drift to the original reading
        tools.CGMs(k) = tools.CGMs(k) + drift;
    end
end

if turnoff_sensor
    start_turnoff_t = 2000;
    turnoff_dt = 250;
    end_turnoff_t = min(t_end, start_turnoff_t + turnoff_dt);

    for k = start_turnoff_t:end_turnoff_t
        tools.CGMs(k) = 0;
    end
end


