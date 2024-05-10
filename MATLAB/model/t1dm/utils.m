classdef utils
    
    properties
        CGMs;
        IIRs;
        CHOs;
        BGs;
        Time;
    end

    methods(Static)         
        function [x0, y_minus1] = init_conditions(params)
            x0 = [0,0,0,params.Gpb,params.Gtb,params.Gpb,params.Ilb,params.Ipb,params.Ib,params.Ib,0,params.Isc1ss,params.Isc2ss];
            y_minus1 = [0,0,0,0,0,false];
        end

        function [x0, y_minus1] = set_init_conditions(z, params)
            Gp = z * params.VG;
            Gt = (-params.Vm0 + params.k1 * Gp * params.Km0) / (params.k2 * params.Km0);
            x0 = [0,0,0,Gp,Gt,Gp,params.Ilb,params.Ipb,params.Ib,params.Ib,0,params.Isc1ss,params.Isc2ss];
            y_minus1 = [0,0,0,0,0,false];
        end

    end

    methods
        function obj = utils(filename)
            if ~strcmp(filename, 'none')
                obj = obj.read_file(filename);
            end
        end


        function obj = read_file(obj, filename)
            dataTable = readtable(filename);
            obj.CGMs = dataTable.CGM; 
            obj.IIRs = dataTable.insulin;
            obj.CHOs = dataTable.CHO;
            obj.BGs = dataTable.BG;
            obj.Time = dataTable.Time;
        end

        function [x_next, y, v] = euler_solve(obj, model, params, x, y_old, u, dt)    
            [y, v] = model.preprocess(x,y_old,u,params,dt);
            dx_dt = model.step(x,y,v,params);
            x_next = (x + dt * dx_dt);
            x_next = x_next .* (x_next >= 0);
        end

        % function [x_next, y, v] = rk4_solve(obj, model, params, x, y_old, u, dt)
        %     vec_x = obj.convert_to_vector(x);
        % 
        %     [y1, v1] = model.preprocess(x,y_old,u,params,dt);
        %     k1 = obj.convert_to_vector(model.step(x,y1,v1,params));
        % 
        %     x2 = obj.convert_to_struct(vec_x + dt / 2 * k1);
        %     [y2, v2] = model.preprocess(x2,y_old,u,params,dt / 2);
        %     k2 = obj.convert_to_vector(model.step(x2,y2,v2,params));
        % 
        %     x3 = obj.convert_to_struct(vec_x + dt / 2 * k2);
        %     [y3, v3] = model.preprocess(x3,y_old,u,params,dt / 2);
        %     k3 = obj.convert_to_vector(model.step(x3,y3,v3,params));
        % 
        %     x4 = obj.convert_to_struct(vec_x + dt * k3);
        %     [y4, v4] = model.preprocess(x4,y_old,u,params,dt);
        %     k4 = obj.convert_to_vector(model.step(x4,y4,v4,params));
        % 
        %     x_next = obj.convert_to_struct(vec_x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4));
        %     y = y1;
        %     v = v1;
        % end

        % function [x_next, y, v] = matlab_solve(obj, model, params, x0, y_old, u, t, dt)
        % 
        %     % This is the auxiliary function that ode45 will call.
        %     % It needs to accept time t and state x, and return the derivative dxdt.
        %     [y, v] = model.preprocess(x0, y_old, u, params, dt);
        %     odeFunction = @(t, x) odeFunctionWrapper(t, x, obj, model, params, y, v);
        % 
        %     % Set up the time span for the ODE solver
        %     tspan = [t, t+dt];
        % 
        %     % Call ode45
        %     vec_x0 = obj.convert_to_vector(x0);
        %     [t, x] = ode45(odeFunction, tspan, vec_x0);
        % 
        %     % The last state returned by ode45 will be at time dt, which is x(end,:)
        %     x_next = obj.convert_to_struct(x(end,:)');
        %     % pause
        % 
        %     % After the ODE solver, you may need to update y and v for the final state
        % 
        % 
        % end
        % 
        % function dxdt = odeFunctionWrapper(t, x_vec, obj, model, params, y, v)
        %     % Convert the vectorized state back into a struct
        %     x = obj.convert_to_struct(x_vec);
        % 
        %     % Compute the derivative using the model's step function
        %     dxdt = obj.convert_to_vector(model.step(x, y, v, params));
        % 
        % end
        % 
        % function x_new = struct_increment(obj, x, dx)
        %     vec_x = obj.convert_to_vector(x);
        %     vec_x_new = vec_x + dx;
        %     x_new = obj.convert_to_struct(vec_x_new);
        % end

    end

end

