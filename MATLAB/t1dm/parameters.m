% Patient features
params.BW = 60; % body weight [kg]
params.eat_rate = 5; % [g/min]

% To be set by physician
params.u2ss = 0.001; % steady state (basal) insulin rate (IIRb)

% To be measured at basal state
params.Gb = 200;

% Sensor features
params.Td = 10; % glucose sensor delay (tbd)

% Glucose Kinetics
params.VG = 1.88;
params.k1 = 0.065;
params.k2 = 0.079;
params.Gpb = params.Gb * params.VG; % basal glucose in plasma
% Insulin Kinetics
params.VI = 0.05;
params.HEb = 0.6;
params.m1 = 0.190;
params.m2 = 0.484;
params.m30 = params.m1 * params.HEb / (1 - params.HEb);
params.m4 = 0.194;
params.m5 = 0.0304;
params.m6 = 0.6471;
params.Ipb = params.u2ss / (params.m2 + params.m4 - params.m1 * params.m2 / (params.m1 + params.m30));  % basal insulin in plasma
params.Ilb = params.m2 / (params.m1 + params.m30) * params.Ipb;
params.Ib = params.Ipb / params.VI;
% Rate of Appearance
params.kmax = 0.0558;
params.kmin = 0.0080;
params.kabs = 0.0570;
params.kgri = params.kmax;
params.f = 0.9;
params.b = 0.82;
params.d = 0.010;
% Endogenous Glucose Production
params.kp1 = 2.70;
params.kp2 = 0.0021;
params.kp3 = 0.009;
params.kp4 = 0.0618;
params.ki = 0.0079;
params.EGPb = params.kp1 - params.kp2 * params.Gpb - params.kp3 * params.Ib;
% Utilization
params.Fcns = 1;
params.Gtb = 1 / params.k2 * (params.Fcns - params.EGPb + params.k1 * params.Gpb);  % basal glucose in slowly equilibrating tissues
params.Km0 = 225.59;
params.Vm0 = (params.EGPb - params.Fcns) * (params.Km0 + params.Gtb) / params.Gtb;
params.Vmx = 0.047;
params.p2U = 0.0331;
% Insulin Infusion
params.kd = 0.0164;
params.ka1 = 0.0018;
params.ka2 = 0.0182;
params.Isc1ss = params.u2ss / (params.kd + params.ka1);
params.Isc2ss = params.Isc1ss * params.kd / params.ka2;
% Renal Excretion
params.ke1 = 0.0005;
params.ke2 = 339;
