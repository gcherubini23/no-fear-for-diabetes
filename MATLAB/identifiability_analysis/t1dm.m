function model = t1dm()

% Define symbolic variables for states, parameters, and inputs
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 y4
syms kgri kmin kmax b c kabs kp1 kp2 kp3 f Fcns ke1 ke2 k1 k2 Vm0 Vmx Km0 m1 m2 m6 m4 ka1 ka2 ki p2u kd VI VG

BW = 102.32;
Ts = 10;
Ib = 51.2;
Gb = 138.56;
% VG = 1.88; 
% VI = 0.05;

Gpb = Gb*VG;
EGPb = kp1 - kp2 * Gpb - kp3 * Ib;
Gtb = 1 / (k2 * (Fcns - EGPb + k1 * Gpb));
basal = 0.0211;
u2ss = basal * 6000 / BW;

eps = 1e-6;

Ipb = Ib * VI;
m3 = m1*m6/(1-m6);
Ilb = m2 / (m1 + m3) * Ipb;
Isc1ss = u2ss / (kd + ka1);
Isc2ss = Isc1ss * kd / ka2;

% States vector
model.sym.x = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12; x13; y4];

% Parameters vector
model.sym.p = [kgri; kmin; kmax; b; c; kabs; kp1; kp2; kp3; f; Fcns; ke1; ke2; k1; k2; Vm0; Vmx; Km0; m1; m2; m6; m4; ka1; ka2; ki; p2u; kd; VI; VG];

% Inputs vector
model.sym.g = [1.0,0;
               0,0;
               0,0;
               0,0;
               0,0;
               0,0;
               0,0;
               0,0;
               0,0;
               0,0;
               0,0;
               0,1.0;
               0,0;
               1.0,0;];

% ODEs
model.sym.xdot = [
    -kgri * x1;
    -x2 * (kmin + (kmax - kmin)/2 * (tanh(5 / (2 * (1-b) * y4 + eps) * (x1 + x2 - b * y4)) - tanh(5 / (2 * c * y4 + eps) * (x1 + x2 - c * y4)) + 2)) + kgri * x1;
    -kabs * x3 + x2 * (kmin + (kmax - kmin)/2 * (tanh(5 / (2 * (1-b) * y4 + eps) * (x1 + x2 - b * y4)) - tanh(5 / (2 * c * y4 + eps) * (x1 + x2 - c * y4)) + 2));
    kp1 - kp2 * x4 - kp3 * x10 + f * kabs * x3 / BW - Fcns - ke1 * (x4 - ke2) - k1 * x4 + k2 * x5;
    -(Vm0 + Vmx * x11) * x5 / (Km0 + x5) + k1 * x4 - k2 * x5;
    -x6 / Ts + x4 / (Ts * VG);
    -(m1 + (m1*m6/(1-m6))) * x7 + m2 * x8;
    -(m2 + m4) * x8 + m1 * x7 + ka1 * x12 + ka2 * x13;
    -ki * (x9 - x8 / (VI));
    -ki * (x10 - x9);
    -p2u * x11 + p2u * (x8 / VI - Ib);
    -(kd + ka1) * x12;
    kd * x12 - ka2 * x13;
    0;
];

% Outputs
model.sym.y = x6/VG;

% Init condition
model.sym.x0 = [0,0,0,Gpb,Gtb,Gb,Ilb,Ipb,Ib,Ib,0,Isc1ss,Isc2ss,0];
% model.sym.x0 = [0,0,0,264,202,139,0,0,0,0,0,0,0,0];

end


