classdef ekf
    properties
    H = [0,0,0,0,0,1,0,0,0,0,0,0,0];

    K = 0;
    Q;  % needs to be tuned depending the dt!!!
    R;
    
    dt;

    state_fields;
    true_input_fields;

    model;
    lin_model;
    tools;

    dx_dts;

    end

    methods
        function obj = ekf(non_linear_model, lin_model, tools, dt, Q, R)           
            obj.state_fields = tools.state_fields;
            obj.true_input_fields = tools.true_input_fields;
            obj.Q = Q;
            obj.R = R;
            obj.dt = dt;
            obj.model = non_linear_model;
            obj.lin_model = lin_model;
            obj.tools = tools;
        end      
        
        function [xp_k, Pp_k, y_kminus1, v_kminus1] = process_update(obj, x_kminus1, y_kminus2, u_kminus1, P_kminus1, params, euler, t)         
            % Euler
            if euler
                [xp_k, y_kminus1, v_kminus1] = obj.tools.euler_solve(obj.model,params,x_kminus1,y_kminus2,u_kminus1,obj.dt);
            else
                % RK4
                % [xp_k, y_kminus1, v_kminus1] = obj.tools.rk4_solve(obj.model,params,x_kminus1,y_kminus2,u_kminus1,obj.dt);

                [xp_k, y_kminus1, v_kminus1] = obj.tools.matlab_solve(obj.model,params,x_kminus1,y_kminus2,u_kminus1,t,obj.dt);
            end

            obj.lin_model = obj.lin_model.linearize(x_kminus1,y_kminus1,params);
            Pp_k = obj.lin_model.A * P_kminus1 * transpose(obj.lin_model.A) + obj.Q;
        end

        function [xm, Pm] = measurement_update(obj, xp, Pp, z)
            vec_xp = obj.tools.convert_to_vector(xp);
            
            obj.K = Pp * transpose(obj.H) * (obj.H * Pp * transpose(obj.H) + obj.R)^(-1);
            vec_xm = vec_xp + obj.K * (z - obj.H * vec_xp);
            Pm = (eye(length(vec_xp)) - obj.K * obj.H) * Pp * transpose(eye(length(vec_xp)) - obj.K * obj.H) + obj.K * obj.R * transpose(obj.K);

            xm = obj.tools.convert_to_struct(vec_xm);
        end

    end
end