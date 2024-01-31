clear
close
clc


%% Parameters
% Glucose Kinetics
params.VG = 1.88;
params.k1 = 0.065;
params.k2 = 0.079;
% Insulin Kinetics
params.VI = 0.05;
params.HEb = 0.6;
params.m1 = 0.190;
params.m2 = 0.484;
params.m30 = params.m1 * params.HEb / (1 - params.HEb);
params.m4 = 0.194;
params.m5 = 0.0304;
params.m6 = 0.6471;
% Rate of Appearance
params.kmax = 0.0558;
params.kmin = 0.0080;
params.kabs = 0.0570;
params.kgri = 0.0558;
params.f = 0.9;
params.a = 0.00013;
params.b = 0.82;
params.c = 0.00236;
params.d = 0.010;
% Endogenous Glucose Production
params.kp1 = 2.70;
params.kp2 = 0.0021;
params.kp3 = 0.009;
params.kp4 = 0.0618;
params.ki = 0.0079;
% Utilization
params.Fcns = 1;
params.Vm0 = 2.5;
params.Vmx = 0.047;
params.Km0 = 225.59;
params.p2U = 0.0331;
% Insulin Secretion
params.K = 2.30;
params.alpha = 0.05;
params.beta = 0.11;
params.gamma = 0.5;
% Insulin Infusion
params.kd = 0.0164;
params.ka1 = 0.0018;
params.ka2 = 0.0182;
% params.Td = 10; % glucose sensor delay tbd
% Renal Excretion
params.ke1 = 0.0005;
params.ke2 = 339;

% Extra
params.BW = 60; % body weight
params.u2ss = 1.5; % steady state insulin rate
params.Ib = 115;
params.ksc = 0.09;
params.kgut = params.kmax;

% Extra for insulin secretion
params.h = 0; % to be set to Gb (basal glucose)
params.Sb = 0;
