classdef utils
    
    properties
        state_fields = {'Qsto1','Qsto2','Qgut','Gp','Gt','Gsc','Il','Ip','Id','I1','X','Isc1','Isc2'};
        true_input_fields = {'CHO_consumed_rate','IIR'};
    end

    methods(Static)         
        function [x0, y0] = init_conditions(params)
            x0.Qsto1 = 0;
            x0.Qsto2 = 0;
            x0.Qgut = 0;
            x0.Gp = params.Gpb;
            x0.Gt = params.Gtb;
            x0.Gsc = params.Gpb;
            x0.Il = params.Ilb;
            x0.Ip = params.Ipb;
            x0.Id = params.Ib;
            x0.I1 = params.Ib;
            x0.X = 0;
            x0.Isc1 = params.Isc1ss;
            x0.Isc2 = params.Isc2ss;

            y0.CHO_to_eat = 0;
            y0.D = 0;
            y0.lastQsto = 0;
            y0.is_eating = false;
        end
        
        function z = sense(filename)
        
        end

        function CHO = announce_meal(filename)

        end

        function IIR = insulin_injection(filename)

        end
    end

    methods
        function obj = utils()
        end
        
        function x_pred = euler_solve(obj, x, x_next, dt)
            vec_x = obj.convert_to_vector(x);
            vec_x_next = obj.convert_to_vector(x_next);

            % Euler method
            vec_x_pred = vec_x + dt * vec_x_next;

            x_pred = obj.convert_to_struct(vec_x_pred);
        end

        function x_pred = rk2_solve()

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


    end

end

