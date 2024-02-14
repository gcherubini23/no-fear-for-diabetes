%% patient_01 true parameters

% Patient features
params.BW = 102.32; % body weight [kg]
params.eat_rate = 5; % [g/min]

% To be set by physician
params.u2ss = 1.2386244136; % steady state (basal) insulin rate (IIRb)

% To be measured at basal state
params.Gb = 138.56;

% Sensor features
params.Td = 10; % glucose sensor delay (tbd)

% Glucose Kinetics
params.VG = 1.9152;
params.k1 = 0.058138;
params.k2 = 0.087114;
params.Gpb = params.Gb * params.VG; % basal glucose in plasma
% Insulin Kinetics
params.VI = 0.054906;
params.HEb = 0.6;
params.m1 = 0.15446;
params.m2 = 0.225;
params.m30 = params.m1 * params.HEb / (1 - params.HEb);
params.m4 = 0.09;
params.m5 = 0.027345;
params.Ipb = params.u2ss / (params.m2 + params.m4 - params.m1 * params.m2 / (params.m1 + params.m30));  % basal insulin in plasma
params.Ilb = params.m2 / (params.m1 + params.m30) * params.Ipb;
params.Ib = params.Ipb / params.VI;
% Rate of Appearance
params.kmax = 0.046122;
params.kmin = 0.0037927;
params.kabs = 0.08906;
params.kgri = params.kmax;
params.f = 0.9;
params.b = 0.70391;
params.d = 0.21057;
% Endogenous Glucose Production
params.kp1 = 4.7314;
params.kp2 = 0.00469;
params.kp3 = 0.01208;
params.ki = 0.0046374;
params.EGPb = params.kp1 - params.kp2 * params.Gpb - params.kp3 * params.Ib;
% Utilization
params.Fcns = 1;
params.Gtb = 1 / params.k2 * (params.Fcns - params.EGPb + params.k1 * params.Gpb);  % basal glucose in slowly equilibrating tissues
params.Km0 = 253.52;
params.Vm0 = (params.EGPb - params.Fcns) * (params.Km0 + params.Gtb) / params.Gtb;
params.Vmx = 0.031319;
params.p2U = 0.0278;
% Insulin Infusion
params.kd = 0.0152;
params.ka1 = 0.0019;
params.ka2 = 0.0078;
params.Isc1ss = params.u2ss / (params.kd + params.ka1);
params.Isc2ss = params.Isc1ss * params.kd / params.ka2;
% Renal Excretion
params.ke1 = 0.0005;
params.ke2 = 339;
