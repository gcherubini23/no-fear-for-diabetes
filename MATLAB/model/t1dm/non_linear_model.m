%% Non-linear model (for T1DM)

classdef non_linear_model

    properties
        tools;
        mA;
        mB;
        mD;
    end

    methods
        function obj = non_linear_model(tools)
            obj.tools = tools;
        end
    
        function dx_dt = step(obj,x,y,v,params)
                
                [dQsto1_dt, dQsto2_dt, dQgut_dt] = obj.gastro_intestinal_tract(x,y,v,params);
                [dGp_dt, dGt_dt, dGsc_dt] = obj.glucose_subystem(x,params);
                [dIl_dt, dIp_dt, dI1_dt, dId_dt, dX_dt, dIsc1_dt, dIsc2_dt] = obj.insulin_infusion_subsystem(x,v,params);
            
                dx_dt = [dQsto1_dt, dQsto2_dt, dQgut_dt, dGp_dt, dGt_dt, dGsc_dt, dIl_dt, dIp_dt, dI1_dt, dId_dt, dX_dt, dIsc1_dt, dIsc2_dt];
        end

        function obj = linearize(obj,x,y,params)
            % Partial derivative for Endogenous glucose production (EGPt)
            if (params.kp1 - params.kp2 * x(4) - params.kp3 * x(10) > 0)
                dEGPt_dGp = -params.kp2;
                dEGPt_dId = -params.kp3;
                bias_Gp_EGPt = params.kp1;
            else
                dEGPt_dGp = 0;
                dEGPt_dId = 0;
                bias_Gp_EGPt = 0;
            end
            
            % Partial derivative for renal excretion (Et)
            if (x(4) - params.ke2) > 0
                dEt_dGp = params.ke1;
                bias_Gp_Et = -params.ke1*params.ke2;
            else
                dEt_dGp = 0;
                bias_Gp_Et = 0;
            end
            
            % Non-linear term Gt
            dUidt_dGt = (params.Vm0 + params.Vmx * x(11)) * params.Km0 / ((params.Km0 + x(5))^2);
            dUidt_dX = params.Vmx * x(5) / (params.Km0 + x(5));
            
            % Non-linear term Kgut 
            Dbar = y(5) + y(4) * 1000;
            % if Dbar > 0
            %     aa = 5 / (2 * (1 - params.b) * Dbar);
            %     cc = 5 / (2 * params.d * Dbar);
            %     T1 = tanh(aa * (x(1) + x(2) - params.b * Dbar));
            %     T2 = tanh(cc * (x(1) + x(2) - params.d * Dbar));
            %     T3 = x(2) * (params.kmax - params.kmin) / 2 * (aa * (1 - T1^2) + cc * (1 - T2^2));
            % 
            %     dQsto2_dQsto1 = params.kgri - T3;
            %     dQsto2_dQsto2 = -params.kmin - (params.kmax - params.kmin) / 2 * (T1 - T2 + 2) - T3;
            % 
            %     dQgut_dQsto1 = T3;
            %     dQgut_dQsto2 = -dQsto2_dQsto2;
            %     dQgut_dQgut = -params.kabs;
            % else
            %     dQsto2_dQsto1 = params.kgri;
            %     dQsto2_dQsto2 = -params.kmax;
            % 
            %     dQgut_dQsto1 = 0;
            %     dQgut_dQsto2 = params.kmax;
            %     dQgut_dQgut = -params.kabs;
            % end

            dQsto2_dQsto1 = params.kgri;
            dQsto2_dQsto2 = -params.kmax;
    
            dQgut_dQsto1 = 0;
            dQgut_dQsto2 = params.kmax;
            dQgut_dQgut = -params.kabs;

            % Compute matrices
            obj.mA = [-params.kgri, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;  % Qsto1                                                             1
                     dQsto2_dQsto1, dQsto2_dQsto2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; % Qsto2                                                 2
                     dQgut_dQsto1, dQgut_dQsto2, dQgut_dQgut, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; % Qgut                                          3

                     0, 0, params.f*params.kabs/params.BW, -params.k1+dEGPt_dGp-dEt_dGp, params.k2, 0, 0, 0, 0, dEGPt_dId, 0, 0, 0; % Gp    4
                     0, 0, 0, params.k1, -params.k2-dUidt_dGt, 0, 0, 0, 0, 0, -dUidt_dX, 0, 0;  % Gt                                        5
                     0, 0, 0, 1/params.Td, 0, -1/params.Td, 0, 0, 0, 0, 0, 0, 0;    % Gpd                                                   6

                     0, 0, 0, 0, 0, 0, -params.m1-params.m30, params.m2, 0, 0, 0, 0, 0;    % Il                                             7
                     0, 0, 0, 0, 0, 0, params.m1, -params.m2-params.m4, 0, 0, 0, params.ka1, params.ka2;    % Ip                            8
                     0, 0, 0, 0, 0, 0, 0, params.ki/params.VI, -params.ki, 0, 0, 0, 0;  % I1                                                9
                     0, 0, 0, 0, 0, 0, 0, 0, params.ki, -params.ki, 0, 0, 0;    % Id                                                        10
                     0, 0, 0, 0, 0, 0, 0, params.p2U/params.VI, 0, 0, -params.p2U, 0, 0;    % X                                             11 

                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -params.kd-params.ka1, 0; % Isc1                                                      12
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, params.kd, -params.ka2];  % Isc2                                                      13

            obj.mD = [0;
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

            obj.mB = [1000, 0;
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

        end

    end

    methods(Static)        
        function [y_k, v_k] = preprocess(x_k,y_kminus1,u_k,params,dt)
            % What are the true inputs? 
            % Insulin
            if y_kminus1(1) <= 0
                IIR_dt = u_k(2);
            else
                IIR_dt = max(u_k(2), y_kminus1(2));    % assuming I do not inject new insulin if I still have some to inject
            end
            
            insulin_to_infuse = y_kminus1(1) + u_k(2);  % [U] as it is multiplied by 1 min
           
            % Check if we have enough insulin to infuse for this time step.
            if insulin_to_infuse < IIR_dt * dt && dt ~= 0
                % If there isn't enough insulin to infuse, just infuse whatever is left.
                IIR_dt = insulin_to_infuse / dt;
                insulin_to_infuse = 0;
            else
                insulin_to_infuse = insulin_to_infuse - IIR_dt * dt;
                if insulin_to_infuse < 1e-5
                    insulin_to_infuse = 0;
                end
            end
            
            % CHO
            if dt == 0
                CHO_consumed_rate = 0;
            elseif (y_kminus1(3) / dt >= params.eat_rate) || (u_k(1) / dt >= params.eat_rate && y_kminus1(3) == 0)
                CHO_consumed_rate = params.eat_rate;
            elseif (u_k(1) > 0) && (u_k(1) / dt < params.eat_rate) && (y_kminus1(3) == 0)
                CHO_consumed_rate = u_k(1) / dt;
            else
                CHO_consumed_rate = y_kminus1(3) / dt;
            end

            % Update extra states
            new_CHO_to_eat = u_k(1) + y_kminus1(3) - CHO_consumed_rate * dt;

            if y_kminus1(6)
                new_D = CHO_consumed_rate * dt + y_kminus1(4);
            elseif CHO_consumed_rate > 0
                new_D = CHO_consumed_rate * dt;
            else
                new_D = y_kminus1(4);
            end

            if CHO_consumed_rate > 0 && y_kminus1(6) == false     % starts eating -> store last state of Qsto
                is_eating = true;
                lastQsto = x_k(1) + x_k(2);
            elseif CHO_consumed_rate == 0 && y_kminus1(6) == true     % stops eating -> restart updating lastQsto
                is_eating = false;
                lastQsto = y_kminus1(5);
            else
                is_eating = y_kminus1(6);
                lastQsto = y_kminus1(5);    
            end

            v_k = [CHO_consumed_rate, IIR_dt];
            y_k = [insulin_to_infuse, IIR_dt, new_CHO_to_eat, new_D, lastQsto, is_eating];

        end
        
        function [dQsto1_dt, dQsto2_dt, dQgut_dt] = gastro_intestinal_tract(x,y,v,params) 
            % Stomach
            dQsto1_dt = (-params.kgri * x(1) + v(1) * 1000);

            Qsto = x(1) + x(2);
            Dbar = y(5) + y(4) * 1000;
            if Dbar > 0
                aa = 5 / (2 * (1 - params.b) * Dbar);
                cc = 5 / (2 * params.d * Dbar);
                kgut = params.kmin + (params.kmax - params.kmin) / 2 * (tanh(aa * (Qsto - params.b * Dbar)) - tanh(cc * (Qsto - params.d * Dbar)) + 2);
            else
                kgut = params.kmax;
            end

            % pause
            
            dQsto2_dt = (params.kgri * x(1) - kgut * x(2));

            % Intestine
            dQgut_dt = kgut * x(2) - params.kabs * x(3);

        end
        
        function [dGp_dt, dGt_dt, dGsc_dt] = glucose_subystem(x,params)
            % Appearance rate of glucose in plasma
            Rat = params.f * params.kabs * x(3) / params.BW;
        
            % Endogenous glucose production
            EGPt = max(params.kp1 - params.kp2 * x(4) - params.kp3 * x(10), 0);
        
            % Glucose utilization
            Uiit = params.Fcns;
            Vmt = params.Vm0 + params.Vmx * x(11);
            Kmt = params.Km0;
            Uidt = (Vmt * x(5)) / (Kmt + x(5));
         
            % Glucose renal excretion
            Et = max(0, params.ke1 * (x(4) - params.ke2));
        
            % Glucose kinetics
            dGp_dt = EGPt + Rat - Uiit - Et - params.k1 * x(4) + params.k2 * x(5); 
           
            dGt_dt = -Uidt + params.k1 * x(4) - params.k2 * x(5);
        
            % Subcutaneous glucose
            dGsc_dt = -1/params.Td * x(6) + 1/params.Td * x(4);

        end
        
        function [dIl_dt, dIp_dt, dI1_dt, dId_dt, dX_dt, dIsc1_dt, dIsc2_dt] = insulin_infusion_subsystem(x,v,params)
            insulin = v(2) * 6000 / params.BW;

            % Liver Insulin kinetics
            dIl_dt = (-(params.m1 + params.m30) * x(7) + params.m2 * x(8));
        
            % Subcutaneous insulin kinetics
            dIsc1_dt = (insulin - (params.kd + params.ka1) * x(12));

            dIsc2_dt = (params.kd * x(12) - params.ka2 * x(13));
            
            % Appearance rate of insulin in plasma
            Rit = params.ka1 * x(12) + params.ka2 * x(13);
            
            % Plasma insulin kinetics (infusion)
            dIp_dt = (-(params.m2 + params.m4) * x(8) + params.m1 * x(7) + Rit);
            if abs(dIp_dt) < 1e-5
                dIp_dt = 0;
            end

            It = x(8) / params.VI;
        
            % Insulin
            dX_dt = params.p2U * (-x(11) + It - params.Ib);  
            
            dI1_dt = params.ki * (It - x(9));

            dId_dt = params.ki * (x(9) - x(10));

        end

    end
end
