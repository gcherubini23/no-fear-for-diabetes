%% Non-linear model (for T1DM)

classdef non_linear_model

    properties
        state_fields = {'Qsto1','Qsto2','Qgut','Gp','Gt','Gsc','Il','Ip','Id','I1','X','Isc1','Isc2'};
        extra_state_fields = {'CHO_to_eat','D','lastQsto','is_eating'}; 
        input_fields = {'CHO', 'IIR'};
        true_input_fields = {'CHO_consumed','IIR'};
    end

    methods
    
        function x_new = step(x,y,v,params,dt)
            
            [new_Qsto1, new_Qsto2, new_Qgut] = gastro_intestinal_tract(x,y,v,params,dt);
            [new_Gp, new_Gt, new_Gsc] = glucose_subystem(x,params,dt);
            [new_Il, new_Ip, new_Id, new_I1, new_X, new_Isc1, new_Isc2] = insulin_infusion_subsystem(x,v,params,dt);
        
            x_new = struct('Qsto1',new_Qsto1,'Qsto2',new_Qsto2,'Qgut',new_Qgut,'Gp',new_Gp,'Gt',new_Gt,'Gsc',new_Gsc,'Il',new_Il,'Ip',new_Ip,'Id',new_Id,'I1',new_I1,'X',new_X,'Isc1',new_Isc1,'Isc2',new_Isc2);
        
        end
        
        function [y, v] = preprocess(x,y_old,u,params,dt)
            % What are the true inputs?
            if (y_old.CHO_to_eat >= params.eat_rate * dt) || (u.CHO >= params.eat_rate * dt && y_old.CHO_to_eat == 0)
                CHO_consumed = params.eat_rate * dt;
            elseif (u.CHO > 0) && (u.CHO < params.eat_rate * dt) && (y_old.CHO_to_eat == 0)
                CHO_consumed = u.CHO;
            else
                CHO_consumed = y_old.CHO_to_eat;
            end
            
            IIR = u.IIR;
            
            % Update extra states
            new_CHO_to_eat = u.CHO + y_old.CHO_to_eat - CHO_consumed;
        
            if CHO_consumed > 0
                new_D = y_old.D + CHO_consumed;
            else
                new_D = 0;    
            end
        
            if CHO_consumed > 0 && y_old.is_eating == false     % starts eating -> store last state of Qsto
                is_eating = true;
                lastQsto = x.Qsto1 + x.Qsto2;
            elseif CHO_consumed == 0 && y_old.is_eating == true     % stops eating -> restart updating lastQsto
                is_eating = false;
                lastQsto = x.Qsto1 + x.Qsto2;    
            else
                if y_old.is_eating
                    lastQsto = y_old.lastQsto;
                else
                    lastQsto = x.Qsto1 + x.Qsto2;
                end
                is_eating = y_old.is_eating;
            end

            v = struct('CHO_consumed',CHO_consumed,'IIR',IIR);
            y = struct('CHO_to_eat',new_CHO_to_eat,'D',new_D,'lastQsto',lastQsto,'is_eating',is_eating);

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
            EGPt = max(params.kp1 - params.kp2 * x.Gp - params.kp3 * x.Id, 0);
        
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
        
        function [new_Il, new_Ip, new_Id, new_I1, new_X, new_Isc1, new_Isc2] = insulin_infusion_subsystem(x,v,params,dt)
            insulin = v.IIR * 6000 / params.BW;
            % basal = params.u2ss * params.BW / 6000;
        
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

    end
end
