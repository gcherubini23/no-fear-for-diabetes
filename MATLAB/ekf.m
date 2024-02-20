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
        
        function [xp, Pp, y, v] = process_update(obj, x_old, y_old, u, P_old, params, euler)         
            % Euler
            if euler
                [xp, y, v] = obj.tools.euler_solve(obj.model,params,x_old,y_old,u,obj.dt);
            else
                % RK4
                [xp, y, v] = obj.tools.rk4_solve(obj.model,params,x_old,y_old,u,obj.dt);
            end

            obj.lin_model = obj.lin_model.linearize(x_old,y_old,params);
            Pp = obj.lin_model.A * P_old * transpose(obj.lin_model.A) + obj.Q;
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