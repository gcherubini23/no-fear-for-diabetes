classdef ekf
    properties
    H;

    K = 0;
    Q;  % needs to be tuned depending the dt!!!
    R;
    
    dt;

    state_fields;
    true_input_fields;

    model;
    lin_model;
    tools;

    y_history = [];
    u_history = [];
    v_history = [];
    t_history = [];

    CGM_MARD = 15;

    end

    methods
        function obj = ekf(non_linear_model, lin_model, tools, params, dt, Q, R)           
            obj.state_fields = tools.state_fields;
            obj.true_input_fields = tools.true_input_fields;
            obj.Q = Q;
            obj.R = R;
            obj.dt = dt;
            obj.model = non_linear_model;
            obj.lin_model = lin_model;
            obj.tools = tools;
            obj.H = [0,0,0,0,0,1/params.VG,0,0,0,0,0,0,0];
        end      
        
        function [xp_k, Pp_k, y_kminus1, v_kminus1] = process_update(obj, x_kminus1, y_kminus2, u_kminus1, P_kminus1, dt, params)         
            % Euler
            [xp_k, y_kminus1, v_kminus1] = obj.tools.euler_solve(obj.model,params,x_kminus1,y_kminus2,u_kminus1,dt);
            obj.lin_model = obj.lin_model.linearize(x_kminus1,y_kminus1,params);
            Pp_k = obj.lin_model.A * P_kminus1 * transpose(obj.lin_model.A) + obj.Q;
        end

        function [xm, Pm, residual, innovation_covariance] = measurement_update(obj, xp, Pp, z)
            vec_xp = obj.tools.convert_to_vector(xp);
            innovation_covariance = obj.H * Pp * transpose(obj.H) + obj.R;
            obj.K = Pp * transpose(obj.H) * innovation_covariance^(-1);
            residual = z - obj.H * vec_xp;
            vec_xm = vec_xp + obj.K * residual;
            KRK = obj.K * obj.R * transpose(obj.K);
            Pm = (eye(length(vec_xp)) - obj.K * obj.H) * Pp * transpose(eye(length(vec_xp)) - obj.K * obj.H) + KRK;
            xm = obj.tools.convert_to_struct(vec_xm);
        end

        function [x_horizon, P_horizon, y_horizon_minus1, v_horizon_minus1] = predict(obj, x_k, y_kminus1, u_k, P_k, horizon, params)
            t = 0;
            x = x_k;
            y = y_kminus1;
            u = u_k;
            P = P_k;
            v.CHO_consumed_rate = 0;
            v.IIR_dt = y_kminus1.last_IIR;

            while t < horizon
                step_dt = min(obj.dt, horizon - t);             
                [x_new, P_new, y_new, v_new] = obj.process_update(x,y,u,P,step_dt,params);
                x = x_new;
                P = P_new;
                y = y_new;
                v = v_new;

                u.CHO = 0;
                u.IIR = 0;
                t = t + step_dt;

            end
   
            x_horizon = x;
            P_horizon = P;
            y_horizon_minus1 = y;
            v_horizon_minus1 = v;

        end

        function [x_horizon, P_horizon, y_horizon_minus1, v_horizon_minus1, obj] = predict_and_save(obj, x_k, y_kminus1, u_k, P_k, horizon, params, t0)
            t = 0;
            x = x_k;
            y = y_kminus1;
            u = u_k;
            P = P_k;
            v.CHO_consumed_rate = 0;
            v.IIR_dt = y_kminus1.last_IIR;
            
            time = t0;

            while t < horizon

                step_dt = min(obj.dt, horizon - t);             
                [x_new, P_new, y_new, v_new] = obj.process_update(x,y,u,P,step_dt,params);
               
                x = x_new;
                P = P_new;
                y = y_new;
                v = v_new;

                if isempty(obj.y_history)
                    obj.u_history = u;
                    obj.y_history = y;
                    obj.v_history = v;
                    obj.t_history = time;
                else
                    obj.u_history(end+1) = u;
                    obj.y_history(end+1) = y;
                    obj.v_history(end+1) = v;
                    obj.t_history(end+1) = time;
                end

                t = t + step_dt;
                time = time + minutes(step_dt);
                u.CHO = 0;
                u.IIR = 0;

            end
   
            x_horizon = x;
            P_horizon = P;
            y_horizon_minus1 = y;
            v_horizon_minus1 = v;

        end

        function obj = update_sensor_cov(obj, z)
            obj.R = (obj.CGM_MARD / 100 * z) * (obj.CGM_MARD / 100 * z);
        end

    end
end