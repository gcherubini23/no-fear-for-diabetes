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
        function [x0, y0] = init_conditions(params)
            x0.Qsto1 = 0;
            x0.Qsto2 = 0;
            x0.Qgut = 0;
            x0.Gp = params.Gpb;
            x0.Gt = params.Gtb;
            x0.Gsc = params.Gb;
            x0.Il = params.Ilb;
            x0.Ip = params.Ipb;
            x0.Id = params.Ib;
            x0.I1 = params.Ib;
            x0.X = 0;
            x0.Isc1 = params.Isc1ss;
            x0.Isc2 = params.Isc2ss;

            y0.insulin_to_infuse = 0;
            y0.last_IIR = params.basal;
            y0.CHO_to_eat = 0;
            y0.D = 0;
            y0.lastQsto = 0;
            y0.is_eating = false;
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

        function x_pred = euler_solve(obj, x, dx_dt, dt)
            vec_x = obj.convert_to_vector(x);
            vec_dx_dt = obj.convert_to_vector(dx_dt);

            % Euler method
            % vec_x_pred = max(0, vec_x + dt * vec_x_next);
            vec_x_pred = vec_x + dt * vec_dx_dt;


            x_pred = obj.convert_to_struct(vec_x_pred);
            
        end

        % function x_pred = rk4_solve(obj, model, params, x, y, v, dt)
        %     u.CHO = 0;
        %     u.IIR = 0;
        % 
        %     % k1 is the slope at the beginning of the interval, using the initial state.
        %     k1 = obj.convert_to_vector(model.step(x, y, v, params));
        % 
        %     % k2 is the slope at the midpoint, using x + dt/2 * k1.
        %     x_k1 = obj.struct_increment(x, 0.5 * dt * k1);
        %     [y_mid, v_mid] = model.preprocess(x_k1, y, u, params, dt / 2);
        %     k2 = obj.convert_to_vector(model.step(x_k1, y_mid, v_mid, params));
        % 
        %     % k3 is also the slope at the midpoint, but now using x + dt/2 * k2.
        %     x_k2 = obj.struct_increment(x, 0.5 * dt * k2);
        %     [y_mid, v_mid] = model.preprocess(x_k2, y, u, params, dt / 2);
        %     k3 = obj.convert_to_vector(model.step(x_k2, y_mid, v_mid, params));
        % 
        %     % k4 is the slope at the end of the interval, using x + dt * k3.
        %     x_k3 = obj.struct_increment(x, dt * k3);
        %     [y_end, v_end] = model.preprocess(x_k3, y, u, params, dt);
        %     k4 = obj.convert_to_vector(model.step(x_k3, y_end, v_end, params));
        % 
        %     % Combine the slopes to estimate the state at the next timestep.
        %     x_pred = obj.struct_increment(x, (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4));
        % end


        function x_new = struct_increment(obj, x, dx)
            vec_x = obj.convert_to_vector(x);
            vec_x_new = vec_x + dx;
            x_new = obj.convert_to_struct(vec_x_new);
        end

    end

end

