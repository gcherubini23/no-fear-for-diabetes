classdef ekf
    properties
    A;
    D;
    B;
    H = [0,0,0,0,0,1,0,0,0,0,0,0,0];

    K = 0;
    Q;  % needs to be tuned depending the dt!!!
    R;

    state_fields;
    true_input_fields = {'CHO_consumed_rate','IIR'};

    end

    methods
        function obj = ekf(params, state_fields, R)           
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
                     0, 0];

            obj.state_fields = state_fields;
            obj.R = R;
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

        function obj = update_matrices(obj, A, D)
           obj.A = A;
           obj.D = D;
        end

        function obj = set_process_noise(obj, Q)
            obj.Q = Q;
        end

        function [xp_new, Pp] = process_update(obj, x, v, P)
            vec_x = obj.convert_to_vector(x);
            vec_v = obj.convert_to_vector(v);

            vec_x_new = obj.A * vec_x + obj.B * vec_v + obj.D;
            Pp = obj.A * P * transpose(obj.A) + obj.Q;
            
            xp_new = obj.convert_to_struct(vec_x_new);
        end

        function [xm, Pm] = measurement_update(obj, xp, Pp, z)
            vec_xp = obj.convert_to_vector(xp);
            
            obj.K = Pp * transpose(obj.H) * (obj.H * Pp * transpose(obj.H) + obj.R)^(-1);
            vec_xm = vec_xp + obj.K * (z - obj.H * vec_xp);
            Pm = (eye(length(vec_xp)) - obj.K * obj.H) * Pp * transpose(eye(length(vec_xp)) - obj.K * obj.H) + obj.K * obj.R * transpose(obj.K);

            xm = obj.convert_to_struct(vec_xm);
        end

    end
end