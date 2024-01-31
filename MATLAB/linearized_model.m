clear
close
clc

%% Parameters definition

state_fields = {'Qsto1','Qsto2','Qgut','Gp','Gt','Gsc','Il','Ip','Id','I1','X','Isc1','Isc2'};
% c1 = cell(length(state_fields),1);
% x = cell2struct(c1,state_fields);
c1 = num2cell(zeros(size(state_fields)));
x = cell2struct(c1,state_fields,2);

input_fields = {'CHO', 'IIR'};
c2 = cell(length(input_fields),1);
u = cell2struct(c2,input_fields);

dt = 0.1;

%% Linearized model (for T1DM)
% x[k+1] = A[k] * x[k] + B[k] * u[k] + D[k]
% z[k] = C * x[k]

% x(t+dt) = x(t) + dt * (A(t) * x(t) + B(t) * u(t) + D(t))
% z(t) = C * x(t)

% Partial derivative for Endogenous glucose production (EGPt)
if (params.kp1 - params.kp2 * x.Gp - params.kp3 * x.Id > 0)
    dEGPt_dGp = -params.kp2;
    dEGPt_dId = -params.kp3;
    bias_Gp_EGPt = params.kp1;
else
    dEGPt_dGp = 0;
    dEGPt_dId = 0;
    bias_Gp_EGPt = 0;
end

% Partial derivative for renal excretion (Et)
if (x.Gp - params.ke2) > 0
    dEt_dGp = params.ke1;
    bias_Gp_Et = -params.ke1*params.ke2;
else
    dEt_dGp = 0;
    bias_Gp_Et = 0;
end

% Non-linear term Gt
dUidt_dGt = (params.Vm0 + params.Vmx * x.X) * params.Km0 / ((params.Km0 + x.Gt)^2);
dUidt_dX = params.Vmx * x.Gt / (params.Km0 + x.Gt);

% Matrices
A = [-params.kmax, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     params.kmax, -params.kgut, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, params.kgut, -params.kabs, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

     0, 0, params.f*params.kabs/params.BW, -params.k1+dEGPt_dGp-dEt_dGp, params.k2, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, params.k1, -params.k2-dUidt_dGt, 0, 0, 0, 0, 0, -dUidt_dX, 0, 0;
     0, 0, 0, params.ksc, 0, -params.ksc, 0, 0, 0, 0, 0, 0, 0;

     0, 0, 0, 0, 0, 0, -params.m1-params.m30, -params.m2, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, params.m1, -params.m2-params.m4, 0, 0, 0, params.ka1, params.ka2;
     0, 0, 0, 0, 0, 0, 0, 0, -params.ki, params.ki, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, params.ki/params.VI, 0, -params.ki, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, params.p2U/params.VI, 0, 0, -params.p2U, 0, 0;

     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -params.kd-params.ka1, 0;
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, params.kd, -params.ka2];

B = [1000, 0;
     0, 0;
     0, 0;
     0, 0;
     0, 0;
     0, 0;
     0, 0;
     0, 0;
     0, 0;
     0, 0;
     0, 0;
     0, 6000/params.BW;
     0, 0;];

C = [0,0,0,0,0,1,0,0,0,0,0,0,0];

D = [0;
     0;
     0;
     bias_Gp_EGPt-params.Fcns-bias_Gp_Et;
     0;
     0;
     0;
     0;
     0;
     0;
     -params.Ib*params.p2U;
     0;
     0];


