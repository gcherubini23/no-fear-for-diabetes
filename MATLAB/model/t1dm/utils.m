classdef utils
    
    properties
        state_fields;
        extra_state_fields; 
        input_fields;
        true_input_fields;
        CGMs;
        IIRs;
        CHOs;
        BGs;
    end

    methods(Static)         
        function [x0, y_minus1] = init_conditions(params)
            x0.Qsto1 = 0;
            x0.Qsto2 = 0;
            x0.Qgut = 0;
            x0.Gp = params.Gpb;
            x0.Gt = params.Gtb;
            x0.Gpd = params.Gpb;
            x0.Il = params.Ilb;
            x0.Ip = params.Ipb;
            x0.I1 = params.Ib;
            x0.Id = params.Ib;
            x0.X = 0;
            x0.Isc1 = params.Isc1ss;
            x0.Isc2 = params.Isc2ss;

            y_minus1.insulin_to_infuse = 0;
            y_minus1.last_IIR = 0;
            y_minus1.CHO_to_eat = 0;
            y_minus1.D = 0;
            y_minus1.lastQsto = 0;
            y_minus1.is_eating = false;
        end

        function [x0, y_minus1] = rand_conditions(params)
            x0.Qsto1 = 0;
            x0.Qsto2 = 0;
            x0.Qgut = 0;
            x0.Gp = 103;
            x0.Gt = 49;
            x0.Gpd = 0;
            x0.Il = 40;
            x0.Ip = 120;
            x0.I1 = 0;
            x0.Id = 0;
            x0.X = 0;
            x0.Isc1 = params.Isc1ss;
            x0.Isc2 = params.Isc2ss;

            y_minus1.insulin_to_infuse = 0;
            y_minus1.last_IIR = 0;
            y_minus1.CHO_to_eat = 0;
            y_minus1.D = 0;
            y_minus1.lastQsto = 0;
            y_minus1.is_eating = false;
        end

    end

    methods
        function obj = utils(filename, state_fields, extra_state_fields, input_fields, true_input_fields)
            obj = obj.read_file(filename);
            obj.state_fields = state_fields;
            obj.extra_state_fields = extra_state_fields; 
            obj.input_fields = input_fields;
            obj.true_input_fields = true_input_fields;
        end

        function vec = convert_to_vector(obj, s)
            % Define field names based on the number of elements in the struct
            if numel(fieldnames(s)) == 2
                fields = obj.true_input_fields;
            else
                fields = obj.state_fields;
            end
            
            % Convert struct elements to a vector
            vec = zeros(length(fields), 1);
            for i = 1:length(fields)
                vec(i) = s.(fields{i});
            end
        end

        function s = convert_to_struct(obj,vec)
            % Define field names based on the length of the vector
            if length(vec) == 2
                fields = obj.true_input_fields;
            else
                fields = obj.state_fields;
            end
            
            % Convert vector elements back to a struct
            s = struct();
            for i = 1:length(fields)
                s.(fields{i}) = vec(i);
            end
        end

        function obj = read_file(obj, filename)
            dataTable = readtable(filename);
            obj.CGMs = dataTable.CGM; 
            obj.IIRs = dataTable.insulin;
            obj.CHOs = dataTable.CHO;
            obj.BGs = dataTable.BG;
        end

        function [x_next, y, v] = euler_solve(obj, model, params, x, y_old, u, dt)    
            [y, v] = model.preprocess(x,y_old,u,params,dt);
            dx_dt = model.step(x,y,v,params);
            vec_x = obj.convert_to_vector(x);
            vec_dx_dt = obj.convert_to_vector(dx_dt);
            vec_x_next = vec_x + dt * vec_dx_dt;
            x_next = obj.convert_to_struct(vec_x_next);
        end

        function [x_next, y, v] = rk4_solve(obj, model, params, x, y_old, u, dt)
            vec_x = obj.convert_to_vector(x);
            
            [y1, v1] = model.preprocess(x,y_old,u,params,dt);
            k1 = obj.convert_to_vector(model.step(x,y1,v1,params));
            
            x2 = obj.convert_to_struct(vec_x + dt / 2 * k1);
            [y2, v2] = model.preprocess(x2,y_old,u,params,dt / 2);
            k2 = obj.convert_to_vector(model.step(x2,y2,v2,params));

            x3 = obj.convert_to_struct(vec_x + dt / 2 * k2);
            [y3, v3] = model.preprocess(x3,y_old,u,params,dt / 2);
            k3 = obj.convert_to_vector(model.step(x3,y3,v3,params));

            x4 = obj.convert_to_struct(vec_x + dt * k3);
            [y4, v4] = model.preprocess(x4,y_old,u,params,dt);
            k4 = obj.convert_to_vector(model.step(x4,y4,v4,params));

            x_next = obj.convert_to_struct(vec_x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4));
            y = y1;
            v = v1;
        end

        function [x_next, y, v] = matlab_solve(obj, model, params, x0, y_old, u, t, dt)
            
            % This is the auxiliary function that ode45 will call.
            % It needs to accept time t and state x, and return the derivative dxdt.
            [y, v] = model.preprocess(x0, y_old, u, params, dt);
            odeFunction = @(t, x) odeFunctionWrapper(t, x, obj, model, params, y, v);
            
            % Set up the time span for the ODE solver
            tspan = [t, t+dt];
            
            % Call ode45
            vec_x0 = obj.convert_to_vector(x0);
            [t, x] = ode45(odeFunction, tspan, vec_x0);
            
            % The last state returned by ode45 will be at time dt, which is x(end,:)
            x_next = obj.convert_to_struct(x(end,:)');
            % pause
            
            % After the ODE solver, you may need to update y and v for the final state
            

        end

        function dxdt = odeFunctionWrapper(t, x_vec, obj, model, params, y, v)
            % Convert the vectorized state back into a struct
            x = obj.convert_to_struct(x_vec);
            
            % Preprocess to get updated y and v based on the current state
            % [y, v] = model.preprocess(x, y_old, u, params, t);
            
            % Compute the derivative using the model's step function
            dxdt = obj.convert_to_vector(model.step(x, y, v, params));

        end

        function x_new = struct_increment(obj, x, dx)
            vec_x = obj.convert_to_vector(x);
            vec_x_new = vec_x + dx;
            x_new = obj.convert_to_struct(vec_x_new);
        end

    end

end

