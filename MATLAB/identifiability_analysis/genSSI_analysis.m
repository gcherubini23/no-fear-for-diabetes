clc
clear all
close all

syms kgri kmin kmax b c kabs kp1 kp2 kp3 f Fcns ke1 ke2 k1 k2 Vm0 Vmx Km0 m1 m2 m6 m4 ka1 ka2 ki p2u kd VI VG


% p = [kgri; kmin; kmax; kabs; kp1; kp2; kp3; ke1; ke2; k1; k2; Vm0; Vmx; Km0; m1; m2; m6; m4; ka1; ka2; ki; p2u; kd];
p = [kp2; k1; k2; kp1; ki; ke1; kmax; kmin; kabs; kp3; Vmx];
options.verbose = true;
options.noRank = false;

n_derivative = 8;
genssiMain('t1dm',n_derivative,p,options);