% IT IS OUTDATED, USE SAME FORMULATION AS T1DM MODEL AND JUST REFER TO THIS
% MODEL FOR THE SECRETION SUBSYSTEM


%% States and inputs definition

state_fields = {'Qsto1','Qsto2','Qgut','Gp','Gt','Gsc','Il','Ip','Id','I1','X','Ipo','Y'};
c1 = cell(length(state_fields),1);
x = cell2struct(c1,state_fields);
x0 = cell2struct(c1,state_fields);

extra_state_fields = {'CHO_to_eat','D','lastQsto','is_eating'};
c2 = cell(length(extra_state_fields),1);
y = cell2struct(c2,extra_state_fields); 

input_fields = {'CHO', 'IIR'};
c3 = cell(length(input_fields),1);
u = cell2struct(c3,input_fields);

true_input_fields = {'CHO_consumed','IIR'};
c4 = cell(length(true_input_fields),1);
v = cell2struct(c4,true_input_fields);

dt = 0.1;

%% Non-linear model (for T2DM/Normal)

function [x_new, y_new] = model(x,y,u,params,dt)
    
    [CHO_consumed, IIR, new_CHO_to_eat, new_D, lastQsto, is_eating] = preprocess(x,y,u,params,dt);
    v = struct('CHO_consumed',CHO_consumed,'IIR',IIR);
    y_new = struct('CHO_to_eat',new_CHO_to_eat,'D',new_D,'lastQsto',lastQsto,'is_eating',is_eating);

    [new_Qsto1, new_Qsto2, new_Qgut] = gastro_intestinal_tract(x,y_new,v,params,dt);
    [new_Gp, new_Gt, new_Gsc] = glucose_subystem(x,params,dt);
    [new_Il, new_Ip, new_Id, new_I1, new_X, new_Ipo, new_Y] = insulin_secretion_subsystem(x,new_Gp,params,dt);

    x_new = struct('Qsto1',new_Qsto1,'Qsto2',new_Qsto2,'Qgut',new_Qgut,'Gp',new_Gp,'Gt',new_Gt,'Gsc',new_Gsc,'Il',new_Il,'Ip',new_Ip,'Id',new_Id,'I1',new_I1,'X',new_X,'Ipo',new_Ipo,'Y',new_Y);

end

function [CHO_consumed, IIR, new_CHO_to_eat, new_D, lastQsto, is_eating] = preprocess(x,y,u,params,dt)
    % What are the true inputs?
    if (y.CHO_to_eat >= params.eat_rate * dt) || (u.CHO >= params.eat_rate * dt && y.CHO_to_eat == 0)
        CHO_consumed = params.eat_rate * dt;
    elseif (u.CHO > 0) && (u.CHO < params.eat_rate * dt) && (y.CHO_to_eat == 0)
        CHO_consumed = u.CHO;
    else
        CHO_consumed = y.CHO_to_eat;
    end
    
    IIR = u.IIR;
    
    % Update extra states
    new_CHO_to_eat = u.CHO + y.CHO_to_eat - CHO_consumed;

    if CHO_consumed > 0
        new_D = y.D + CHO_consumed;
    else
        new_D = 0;    
    end

    if CHO_consumed > 0 && y.is_eating == false     % starts eating -> store last state of Qsto
        is_eating = true;
        lastQsto = x.Qsto1 + x.Qsto2;
    elseif CHO_consumed == 0 && y.is_eating == true     % stops eating -> restart updating lastQsto
        is_eating = false;
        lastQsto = x.Qsto1 + x.Qsto2;    
    else
        if y.is_eating
            lastQsto = y.lastQsto;
        else
            lastQsto = x.Qsto1 + x.Qsto2;
        end
        is_eating = y.is_eating;
    end
end

function [new_Qsto1, new_Qsto2, new_Qgut] = gastro_intestinal_tract(x,y,v,params,dt) 
    % Stomach
    new_Qsto1 = x.Qsto1 + dt * (-params.kgri * x.Qsto1) + v.CHO_consumed * 1000;
    
    Qsto = x.Qsto1 + x.Qsto2;
    Dbar = y.lastQsto + y.D;
    if Dbar > 0
        aa = 5 / (2 * (1 - params.b) * Dbar);
        cc = 5 / (2 * params.d * Dbar);
        kgut = params.kmin + (params.kmax - params.kmin) / 2 * (tanh(aa * (Qsto - params.b * Dbar)) - tanh(cc * (Qsto - params.d * Dbar)) + 2);
    else
        kgut = params.kmax;
    end

    new_Qsto2 = x.Qsto2 + dt * (params.kgri * x.Qsto1 - kgut * x.Qsto2);

    % Intestine
    new_Qgut = x.Qgut + dt * (kgut * x.Qsto2 - params.kabs * x.Qgut);

end

function [new_Gp, new_Gt, new_Gsc] = glucose_subystem(x,params,dt)
    % Appearance rate of glucose in plasma
    Rat = params.f * params.kabs * x.Qgut / params.BW;

    % Endogenous glucose production
    EGPt = max(params.kp1 - params.kp2 * x.Gp - params.kp3 * x.Id - params.kp4 * x.Ipo, 0);

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
    new_Gsc = x.Gsc + dt * (-1/params.Td * x.Gsc + 1/params.Td * x.Gp);
    if new_Gsc < 0
        new_Gsc = 0;
    end

end

function [new_Il, new_Ip, new_Id, new_I1, new_X, new_Ipo, new_Y] = insulin_secretion_subsystem(x,new_Gp,params,dt)
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
