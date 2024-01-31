clear
close
clc

%% Parameters definition

state_fields = {'Qsto1','Qsto2','Qgut','Gp','Gt','Gsc','Il','Ip','Id','I1','X','Isc1','Isc2'};
c1 = cell(length(state_fields),1);
x = cell2struct(c1,state_fields);
x0 = cell2struct(c1,state_fields);

input_fields = {'CHO', 'IIR'};
c2 = cell(length(input_fields),1);
u = cell2struct(c2,input_fields);

dt = 0.1;

%% Non-linear model (for T1DM)

function x_new = model(x,u,lastQsto,lastfoodtaken,params,dt)
    
    [new_Qsto1, new_Qsto2, new_Qgut] = gastro_intestinal_tract(x,u,lastQsto,lastfoodtaken,params,dt);
    [new_Gp, new_Gt, new_Gsc] = glucose_subystem(x,params,dt);
    [new_Il, new_Ip, new_Id, new_I1, new_X, new_Isc1, new_Isc2] = insulin_infusion_subsystem(x,u,params,dt);

    x_new = struct('Qsto1',new_Qsto1,'Qsto2',new_Qsto2,'Qgut',new_Qgut,'Gp',new_Gp,'Gt',new_Gt,'Gsc',new_Gsc,'Il',new_Il,'Ip',new_Ip,'Id',new_Id,'I1',new_I1,'X',new_X,'Isc1',new_Isc1,'Isc2',new_Isc2);

end

function [new_Qsto1, new_Qsto2, new_Qgut] = gastro_intestinal_tract(x,u,lastQsto,lastfoodtaken,params,dt)
    d = u.CHO * 1000; % g->mg

    % Stomach
    Qsto = x.Qsto1 + x.Qsto2;
    Dbar = lastQsto + lastfoodtaken * 1000;
    
    new_Qsto1 = x.Qsto1 + dt * (-params.kmax * x.Qsto1 + d);
    
    if Dbar > 0
        aa = 5 / 2 / (1 - params.b) / Dbar;
        cc = 5 / 2 / params.d / Dbar;
        params.kgut = params.kmin + (params.kmax - params.kmin) / 2 * (tanh(aa * (Qsto - params.b * Dbar)) - tanh(cc * (Qsto - params.d * Dbar)) + 2);
    else
        params.kgut = params.kmax;
    end

    new_Qsto2 = x.Qsto2 + dt * (params.kmax * x.Qsto1 - params.kgut * x.Qsto2);

    % Intestine
    new_Qgut = x.Qgut + dt * (params.kgut * x.Qsto2 - params.kabs * x.Qgut);

end

function [new_Gp, new_Gt, new_Gsc] = glucose_subystem(x,params,dt)
    % Appearance rate of glucose in plasma
    Rat = params.f * params.kabs * x.Qgut / params.BW;

    % Endogenous glucose production
    EGPt = max(params.kp1 - params.kp2 * x.Gp - params.kp3 * x.Id, 0); % if insulin is secreted: -params.kp4*x.Ipo

    % Glucose utilization
    Uiit = params.Fcns;
    Vmt = params.Vm0 + params.Vmx * x.X;
    Kmt = params.Km0;
    Uidt = (Vmt * x.Gt) / (Kmt + x.Gt);

    % Glucose renal excretion
    Et = max(params.ke1 * (x.Gp - params.ke2), 0);

    % Glucose kinetics
    new_Gp = x.Gp + dt * (EGPt + Rat - Uiit - Et - params.k1 * x.Gp + params.k2 * x.Gt);
    if new_Gp < 0
        new_Gp = 0;
    end

    new_Gt = x.Gt + dt * (-Uidt + params.k1 * x.Gp - params.k2 * x.Gt);
    if new_Gt < 0
        new_Gt = 0;
    end

    % Subcutaneous glucose
    new_Gsc = x.Gsc + dt * (-params.ksc * x.Gsc + params.ksc * x.Gp);
    if new_Gsc < 0
        new_Gsc = 0;
    end

end

function [new_Il, new_Ip, new_Id, new_I1, new_X, new_Isc1, new_Isc2] = insulin_infusion_subsystem(x,u,params,dt)
    insulin = u.IIR * 6000 / params.BW;
    basal = params.u2ss * params.BW / 6000; % ask why

    % Liver Insulin kinetics
    new_Il = x.Il + dt * (-(params.m1 + params.m30) * x.Il + params.m2 * x.Ip);
    if new_Il < 0
        new_Il = 0;
    end

    % Subcutaneous insulin kinetics
    new_Isc1 = x.Isc1 + dt * (insulin - (params.kd + params.ka1) * x.Isc1);
    if new_Isc1 < 0
        new_Isc1 = 0;
    end

    new_Isc2 = x.Isc2 + dt * (params.kd * x.Isc1 - params.ka2 * x.Isc2);
    if new_Isc2 < 0
        new_Isc2 = 0;
    end
    
    % Appearance rate of insulin in plasma
    Rit = params.ka1 * x.Isc1 + params.ka2 * x.Isc2;
    
    % Plasma insulin kinetics (infusion)
    new_Ip = x.Ip + dt * (-(params.m2 + params.m4) * x.Ip + params.m1 * x.Il + Rit);
    It = x.Ip / params.VI;
    if new_Ip < 0
        new_Ip = 0;
    end

    % Insulin action on glucose utilization
    new_X = x.X + dt * (params.p2U * (-x.X + It - params.Ib));
    
    % Insulin action on glucose production
    new_I1 = x.I1 + dt * (-params.ki * (x.I1 - It));
    new_Id = x.Id + dt * (-params.ki * (x.Id - x.I1));

end

%% Secretion subsystem for normal/T2DM patient

function [new_Il, new_Ip, new_Id, new_I1, new_X, new_Ipo, new_Y] = insulin_secretion_subsystem(x, new_Gp,params,dt)
    Gt = x.Gp / params.VG;
    new_G = new_Gp / params.VG;
    
    % Beta-cells insulin
    if params.beta * (Gt - params.h) >= -params.Sb
        new_Y = x.Y + dt * (-params.alpha * (x.Y - params.beta * (Gt - params.h)));
    else
        new_Y = x.Y + dt * (-params.alpha * (x.Y - params.alpha * params.Sb));
    end
    
    if new_G > 0
        Spot = x.Y + params.K * new_G + params.Sb;
    else
        Spot = x.Y + params.Sb;
    end

    St = params.gamma * x.Ipo;
    new_Ipo = x.Ipo + dt * (-St + Spot);

    % Liver Insulin kinetics (secretion)
    HEt = -params.m5 * St + params.m6;
    m3t = HEt * params.m1 / (1 - HEt);
    new_Il = x.Il + dt * (-(params.m1 + m3t) * x.Il + params.m2 * x.Ip + St);
    if new_Il < 0
        new_Il = 0;
    end

    % Plasma insulin kinetics
    new_Ip = x.Ip + dt * (-(params.m2 + params.m4) * x.Ip + params.m1 * x.Il);
    It = x.Ip / params.VI;
    if new_Ip < 0
        new_Ip = 0;
    end
    
    % Insulin action on glucose utilization
    new_X = x.X + dt * (params.p2U * (-x.X + It - params.Ib));

    % Insulin action on glucose production
    new_I1 = x.I1 + dt * (-params.ki * (x.I1 - It));
    new_Id = x.Id + dt * (-params.ki * (x.Id - x.I1));

    
end


