%% Non-linear model (for T1DM)

classdef non_linear_model

    properties
        state_fields = {'Qsto1','Qsto2','Qgut','Gp','Gt','Gsc','Il','Ip','Id','I1','X','Isc1','Isc2'};
        extra_state_fields = {'CHO_to_eat','D','lastQsto','is_eating'}; 
        input_fields = {'CHO', 'IIR'};
        true_input_fields = {'CHO_consumed_rate','IIR'};
    end

    methods
    
        function dx_dt = step(x,y,v,params)
            
            [dQsto1_dt, dQsto2_dt, dQgut_dt] = gastro_intestinal_tract(x,y,v,params);
            [dGp_dt, dGt_dt, dGsc_dt] = glucose_subystem(x,params);
            [dIl_dt, dIp_dt, dI1_dt, dId_dt, dX_dt, dIsc1_dt, dIsc2_dt] = insulin_infusion_subsystem(x,v,params);
        
            dx_dt = struct('Qsto1',dQsto1_dt,'Qsto2',dQsto2_dt,'Qgut',dQgut_dt,'Gp',dGp_dt,'Gt',dGt_dt,'Gsc',dGsc_dt,'Il',dIl_dt,'Ip',dIp_dt,'Id',dId_dt,'I1',dI1_dt,'X',dX_dt,'Isc1',dIsc1_dt,'Isc2',dIsc2_dt);
        end
        
        function [y, v] = preprocess(x,y_old,u,params,dt)
            % What are the true inputs?
            if (y_old.CHO_to_eat / dt >= params.eat_rate) || (u.CHO / dt >= params.eat_rate && y_old.CHO_to_eat == 0)
                CHO_consumed_rate = params.eat_rate;
            elseif (u.CHO > 0) && (u.CHO / dt < params.eat_rate) && (y_old.CHO_to_eat == 0)
                CHO_consumed_rate = u.CHO / dt;
            else
                CHO_consumed_rate = y_old.CHO_to_eat / dt;
            end
            
            IIR = u.IIR;
            
            % Update extra states
            new_CHO_to_eat = u.CHO + y_old.CHO_to_eat - CHO_consumed_rate * dt;
        
            if CHO_consumed_rate > 0
                new_D = y_old.D + CHO_consumed_rate * dt;
            else
                new_D = 0;    
            end
        
            if CHO_consumed_rate > 0 && y_old.is_eating == false     % starts eating -> store last state of Qsto
                is_eating = true;
                lastQsto = x.Qsto1 + x.Qsto2;
            elseif CHO_consumed_rate == 0 && y_old.is_eating == true     % stops eating -> restart updating lastQsto
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

            v = struct('CHO_consumed_rate',CHO_consumed_rate,'IIR',IIR);
            y = struct('CHO_to_eat',new_CHO_to_eat,'D',new_D,'lastQsto',lastQsto,'is_eating',is_eating);

        end
        
        function [dQsto1_dt, dQsto2_dt, dQgut_dt] = gastro_intestinal_tract(x,y,v,params) 
            % Stomach
            dQsto1_dt = -params.kgri * x.Qsto1 + v.CHO_consumed_rate * 1000;
            
            Qsto = x.Qsto1 + x.Qsto2;
            Dbar = y.lastQsto + y.D;
            if Dbar > 0
                aa = 5 / (2 * (1 - params.b) * Dbar);
                cc = 5 / (2 * params.d * Dbar);
                kgut = params.kmin + (params.kmax - params.kmin) / 2 * (tanh(aa * (Qsto - params.b * Dbar)) - tanh(cc * (Qsto - params.d * Dbar)) + 2);
            else
                kgut = params.kmax;
            end
        
            dQsto2_dt = params.kgri * x.Qsto1 - kgut * x.Qsto2;
        
            % Intestine
            dQgut_dt = kgut * x.Qsto2 - params.kabs * x.Qgut;
        
        end
        
        function [dGp_dt, dGt_dt, dGsc_dt] = glucose_subystem(x,params)
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
            dGp_dt = max(0, EGPt + Rat - Uiit - Et - params.k1 * x.Gp + params.k2 * x.Gt);
            
            dGt_dt = max(0, -Uidt + params.k1 * x.Gp - params.k2 * x.Gt);
        
            % Subcutaneous glucose
            dGsc_dt = max(0, -1/params.Td * x.Gsc + 1/params.Td * x.Gp);
        
        end
        
        function [dIl_dt, dIp_dt, dId_dt, dI1_dt, dX_dt, dIsc1_dt, dIsc2_dt] = insulin_infusion_subsystem(x,v,params)
            insulin = v.IIR * 6000 / params.BW;
        
            % Liver Insulin kinetics
            dIl_dt = max(0, -(params.m1 + params.m30) * x.Il + params.m2 * x.Ip);
        
            % Subcutaneous insulin kinetics
            dIsc1_dt = max(0, insulin - (params.kd + params.ka1) * x.Isc1);
        
            dIsc2_dt = max(0, params.kd * x.Isc1 - params.ka2 * x.Isc2);
            
            % Appearance rate of insulin in plasma
            Rit = params.ka1 * x.Isc1 + params.ka2 * x.Isc2;
            
            % Plasma insulin kinetics (infusion)
            dIp_dt = max(0, -(params.m2 + params.m4) * x.Ip + params.m1 * x.Il + Rit);
            It = x.Ip / params.VI;
        
            % Insulin action on glucose utilization
            dX_dt = params.p2U * (-x.X + It - params.Ib);
            
            % Insulin action on glucose production
            dI1_dt = -params.ki * (x.I1 - It);
            dId_dt = -params.ki * (x.Id - x.I1);
        end

    end
end
