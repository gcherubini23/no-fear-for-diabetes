classdef patient_1002
    properties
        % Patient features
        BW = 65; % body weight [kg]
        eat_rate = 5; % [g/min]
        
        % To be set by physician
        basal;
        u2ss; % steady state insulin rate (IIRb)    
        
        % To be measured at basal state
        % Gb = 237.6;
        Gb = 100;
        
        % Sensor features
        Td = 10; % glucose sensor delay
        
        % Glucose Kinetics
        VG = 1.88;
        k1 = 0.065;
        k2 = 0.079;
        Gpb; % basal glucose in plasma
        % Insulin Kinetics
        VI = 0.05;
        HEb = 0.6;
        CL;
        m1 = 0.190;
        m2 = 0.484;
        m30;
        m4 = 0.194;
        m5 = 0.0304;
        Ipb;  % basal insulin in plasma
        Ilb;
        Ib = 78;
        % Rate of Appearance
        kmax = 0.0558;
        kmin = 0.0080;
        kabs = 0.0570;
        kgri;
        f = 0.9;
        b = 0.82;
        d = 0.010;
        % Endogenous Glucose Production
        kp1 = 2.70;
        kp2 = 0.0021;
        kp3 = 0.009;
        ki = 0.01;
        EGPb;
        % Utilization
        Fcns = 1;
        Gtb;  % basal glucose in slowly equilibrating tissues
        Km0 = 225.59;
        Vm0 = 2.5;
        Vmx = 0.047;
        p2U = 0.0331;
        % Insulin Infusion
        kd = 0.0164;
        ka1 = 0.0018;
        ka2 = 0.0182;
        Isc1ss;
        Isc2ss;
        % Renal Excretion
        ke1 = 0.0005;
        ke2 = 339;
    end

    methods
        function obj = patient_1002()
            obj.m30 = obj.m1 * obj.HEb / (1 - obj.HEb);
            obj.CL = 0.0242 * obj.BW;
            obj.m2 = 3/5 * obj.CL / (obj.HEb * obj.VI * obj.BW);
            obj.m4 = 2/5 * obj.CL / (obj.VI * obj.BW);
            % obj.Ipb = obj.u2ss / (obj.m2 + obj.m4 - obj.m1 * obj.m2 / (obj.m1 + obj.m30));  % basal insulin in plasma
            obj.kgri = obj.kmax;

            % obj.Ib = obj.Ipb / obj.VI;
            obj.Ipb = obj.Ib * obj.VI;
            obj.Ilb = obj.m2 / (obj.m1 + obj.m30) * obj.Ipb;
            obj.u2ss = obj.Ipb * (obj.m2 + obj.m4 - obj.m1 * obj.m2 / (obj.m1 + obj.m30));
            
            obj.Gpb = obj.Gb * obj.VG; % basal glucose in plasma

            obj.EGPb = obj.kp1 - obj.kp2 * obj.Gpb - obj.kp3 * obj.Ib;
            obj.Gtb = 1 / obj.k2 * (obj.Fcns - obj.EGPb + obj.k1 * obj.Gpb);
            obj.Vm0 = (obj.EGPb - obj.Fcns) * (obj.Km0 + obj.Gtb) / obj.Gtb;
            obj.Isc1ss = obj.u2ss / (obj.kd + obj.ka1);
            obj.Isc2ss = obj.Isc1ss * obj.kd / obj.ka2;

        end

        function obj = recompile(obj)

            obj.m30 = obj.m1 * obj.HEb / (1 - obj.HEb);
            obj.m2 = 3/5 * obj.CL / (obj.HEb * obj.VI * obj.BW);
            obj.m4 = 2/5 * obj.CL / (obj.VI * obj.BW);
            % obj.Ipb = obj.u2ss / (obj.m2 + obj.m4 - obj.m1 * obj.m2 / (obj.m1 + obj.m30));  % basal insulin in plasma
            obj.kgri = obj.kmax;

            % obj.Ib = obj.Ipb / obj.VI;
            obj.Ipb = obj.Ib * obj.VI;
            obj.Ilb = obj.m2 / (obj.m1 + obj.m30) * obj.Ipb;
            obj.u2ss = obj.Ipb * (obj.m2 + obj.m4 - obj.m1 * obj.m2 / (obj.m1 + obj.m30));
            
            obj.Gpb = obj.Gb * obj.VG; % basal glucose in plasma

            obj.EGPb = obj.kp1 - obj.kp2 * obj.Gpb - obj.kp3 * obj.Ib;
            obj.Gtb = 1 / obj.k2 * (obj.Fcns - obj.EGPb + obj.k1 * obj.Gpb);
            obj.Vm0 = (obj.EGPb - obj.Fcns) * (obj.Km0 + obj.Gtb) / obj.Gtb;
            obj.Isc1ss = obj.u2ss / (obj.kd + obj.ka1);
            obj.Isc2ss = obj.Isc1ss * obj.kd / obj.ka2;
            
        end

        function obj = set_params(obj, params_to_estimate, p)
            fields = params_to_estimate; 
            for i = 1:numel(fields)
                fieldName = fields{i}; 
                if isprop(obj, fieldName) 
                    obj.(fieldName) = p(i); 
                else
                    warning('Property %s does not exist in patient_11.', fieldName);
                end
            end

            obj = obj.recompile();

        end

    end

end
