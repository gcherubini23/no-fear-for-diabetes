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
            x0.I1 = params.Ib;
            x0.Id = params.Ib;
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

        function [x, y, v] = euler_solve(obj, model, params, x_old, y_old, u, dt)    
            [y, v] = model.preprocess(x_old,y_old,u,params,dt);
            dx_dt = model.step(x_old,y,v,params);
            vec_x_old = obj.convert_to_vector(x_old);
            vec_dx_dt = obj.convert_to_vector(dx_dt);
            vec_x_pred = vec_x_old + dt * vec_dx_dt;
            x = obj.convert_to_struct(vec_x_pred);
        end

        function [x, y, v] = rk4_solve(obj, model, params, x_old, y_old, u, dt)
            vec_x_old = obj.convert_to_vector(x_old);
            
            [y1, v1] = model.preprocess(x_old,y_old,u,params,dt);
            k1 = obj.convert_to_vector(model.step(x_old,y1,v1,params));
            
            x2 = obj.convert_to_struct(vec_x_old + dt / 2 * k1);    
            [y2, v2] = model.preprocess(x2,y_old,u,params,dt / 2);
            k2 = obj.convert_to_vector(model.step(x2,y2,v2,params));

            % pause

            x3 = obj.convert_to_struct(vec_x_old + dt / 2 * k2);
            [y3, v3] = model.preprocess(x3,y_old,u,params,dt / 2);
            k3 = obj.convert_to_vector(model.step(x3,y3,v3,params));

            x4 = obj.convert_to_struct(vec_x_old + dt * k3);
            [y4, v4] = model.preprocess(x4,y_old,u,params,dt);
            k4 = obj.convert_to_vector(model.step(x4,y4,v4,params));

            x = obj.convert_to_struct(vec_x_old + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4));
            y = y1;
            v = v1;

            % euler_step = dt * k1
            % rk4_step = dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            % pause

        end



        function x_new = struct_increment(obj, x, dx)
            vec_x = obj.convert_to_vector(x);
            vec_x_new = vec_x + dx;
            x_new = obj.convert_to_struct(vec_x_new);
        end

    end

end

