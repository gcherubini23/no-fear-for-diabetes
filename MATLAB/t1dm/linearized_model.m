%% Linearized model (for T1DM)

% x[k+1] = A[k] * x[k] + B[k] * v[k] + D[k]
% v[k] = g(x,y,u)
% z[k] = H * x[k]

classdef linearized_model
    properties
        A;
        B;
        D;
        state_fields = {'Qsto1','Qsto2','Qgut','Gp','Gt','Gsc','Il','Ip','I1','Id','X','Isc1','Isc2'};
        extra_state_fields = {'insulin_to_infuse','last_IIR','CHO_to_eat','D','lastQsto','is_eating'};
        input_fields = {'CHO', 'IIR'};
        true_input_fields = {'CHO_consumed_rate','IIR_dt'};
        tools;
    end

    methods
        function obj = linearized_model(tools)
            obj.tools = tools;
        end      
 
        function dx_dt = step(obj,x,y,v,params)
            % Linearize
            obj = obj.linearize(x,y,params);
            
            % Step using lin model
            vec_x = obj.tools.convert_to_vector(x);
            vec_v = obj.tools.convert_to_vector(v);
            vec_dx_dt = obj.A * vec_x + obj.B * vec_v + obj.D;  % to ensure states non negativity TBD
            dx_dt = obj.tools.convert_to_struct(vec_dx_dt);
        end
        
        function obj = linearize(obj,x,y,params)
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
            
            % Non-linear term Kgut 
            Dbar = y.lastQsto + y.D;
            if Dbar > 0
                aa = 5 / (2 * (1 - params.b) * Dbar);
                cc = 5 / (2 * params.d * Dbar);
                T1 = tanh(aa * (x.Qsto1 + x.Qsto2 - params.b * Dbar));
                T2 = tanh(cc * (x.Qsto1 + x.Qsto2 - params.d * Dbar));
                T3 = x.Qsto2 * (params.kmax - params.kmin) / 2 * (aa * (1 - T1^2) + cc * (1 - T2^2));
                
                dQsto2_dQsto1 = params.kgri - T3;
                dQsto2_dQsto2 = -params.kmin - (params.kmax - params.kmin) / 2 * (T1 - T2 + 2) - T3;
            
                dQgut_dQsto1 = T3;
                dQgut_dQsto2 = -dQsto2_dQsto2;
                dQgut_dQgut = -params.kabs;
            else
                dQsto2_dQsto1 = params.kgri;
                dQsto2_dQsto2 = -params.kmax;
        
                dQgut_dQsto1 = 0;
                dQgut_dQsto2 = params.kmax;
                dQgut_dQgut = -params.kabs;
            end

            % Compute matrices
            obj.A = [-params.kgri, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;  % Qsto1
                     dQsto2_dQsto1, dQsto2_dQsto2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; % Qsto2
                     dQgut_dQsto1, dQgut_dQsto2, dQgut_dQgut, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; % Qgut
                
                     0, 0, params.f*params.kabs/params.BW, -params.k1+dEGPt_dGp-dEt_dGp, params.k2, 0, 0, 0, dEGPt_dId, 0, 0, 0, 0; % Gp
                     0, 0, 0, params.k1, -params.k2-dUidt_dGt, 0, 0, 0, 0, 0, -dUidt_dX, 0, 0;  % Gt
                     0, 0, 0, 1/(params.Td * params.VG), 0, -1/params.Td, 0, 0, 0, 0, 0, 0, 0;    % Gsc
                
                     0, 0, 0, 0, 0, 0, -params.m1-params.m30, -params.m2, 0, 0, 0, 0, 0;    % Il
                     0, 0, 0, 0, 0, 0, params.m1, -params.m2-params.m4, 0, 0, 0, params.ka1, params.ka2;    % Ip
                     0, 0, 0, 0, 0, 0, 0, params.ki/params.VI, 0, -params.ki, 0, 0, 0;  % I1
                     0, 0, 0, 0, 0, 0, 0, 0, -params.ki, params.ki, 0, 0, 0;    % Id
                     0, 0, 0, 0, 0, 0, 0, params.p2U/params.VI, 0, 0, -params.p2U, 0, 0;    % X
                
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -params.kd-params.ka1, 0; % Isc1
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, params.kd, -params.ka2];  % Isc2
            
            obj.D = [0;
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
            
            obj.B = [1000, 0;
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
        function [y, v] = preprocess(x,y_old,u,params,dt)
            % What are the true inputs?
            
            % Insulin
            if y_old.insulin_to_infuse <= 0
                IIR_dt = u.IIR;
            else
                IIR_dt = y_old.last_IIR;
            end
            
            insulin_to_infuse = y_old.insulin_to_infuse + u.IIR;  
           
            % Check if we have enough insulin to infuse for this time step.
            if insulin_to_infuse < IIR_dt * dt
                % If there isn't enough insulin to infuse, just infuse whatever is left.
                IIR_dt = insulin_to_infuse / dt;
            end

            new_insulin_to_infuse = max(0, insulin_to_infuse - IIR_dt * dt);
            epsilon = 1e-5;
            if new_insulin_to_infuse <= epsilon
                new_insulin_to_infuse = 0;
            end
            
            % CHO
            if (y_old.CHO_to_eat / dt >= params.eat_rate) || (u.CHO / dt >= params.eat_rate && y_old.CHO_to_eat == 0)
                CHO_consumed_rate = params.eat_rate;
            elseif (u.CHO > 0) && (u.CHO / dt < params.eat_rate) && (y_old.CHO_to_eat == 0)
                CHO_consumed_rate = u.CHO / dt;
            else
                CHO_consumed_rate = y_old.CHO_to_eat / dt;
            end

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

            v = struct('CHO_consumed_rate',CHO_consumed_rate,'IIR_dt',IIR_dt);
            y = struct('insulin_to_infuse',new_insulin_to_infuse,'last_IIR',IIR_dt,'CHO_to_eat',new_CHO_to_eat,'D',new_D,'lastQsto',lastQsto,'is_eating',is_eating);

        end
    end
end

