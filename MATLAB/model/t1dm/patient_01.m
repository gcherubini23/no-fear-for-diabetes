%% patient_01 true parameters

classdef patient_01
    properties
        % Patient features
        BW = 102.32; % body weight [kg]
        eat_rate = 5; % [g/min]
        
        % To be set by physician
        basal;
        u2ss = 1.2386244136; % steady state (basal) insulin rate (IIRb)
        % u2ss;      
        
        % To be measured at basal state
        Gb = 138.56;
        
        % Sensor features
        Td = 1 / 0.0766; % glucose sensor delay
        
        % Glucose Kinetics
        VG = 1.9152;
        k1 = 0.058138;
        k2 = 0.087114;
        Gpb; % basal glucose in plasma
        % Insulin Kinetics
        VI = 0.054906;
        HEb = 0.6;
        m1 = 0.15446;
        m2 = 0.225;
        m30;
        m4 = 0.09;
        m5 = 0.027345;
        Ipb;  % basal insulin in plasma
        Ilb;
        Ib;
        % Rate of Appearance
        kmax = 0.046122;
        kmin = 0.0037927;
        kabs = 0.08906;
        kgri;
        f = 0.9;
        b = 0.70391;
        d = 0.21057;
        % Endogenous Glucose Production
        kp1 = 4.7314;
        kp2 = 0.00469;
        kp3 = 0.01208;
        ki = 0.0046374;
        EGPb;
        % Utilization
        Fcns = 1;
        Gtb;  % basal glucose in slowly equilibrating tissues
        Km0 = 253.52;
        Vm0;
        Vmx = 0.031319;
        p2U = 0.0278;
        % Insulin Infusion
        kd = 0.0152;
        ka1 = 0.0019;
        ka2 = 0.0078;
        Isc1ss;
        Isc2ss;
        % Renal Excretion
        ke1 = 0.0005;
        ke2 = 339;
    end

    methods
        function obj = patient_01(basal)
            obj.basal = basal;
            obj.u2ss = basal * 6000 / obj.BW;

            % obj.basal = obj.u2ss * obj.BW / 6000;

            obj.m30 = obj.m1 * obj.HEb / (1 - obj.HEb);
            obj.Ipb = obj.u2ss / (obj.m2 + obj.m4 - obj.m1 * obj.m2 / (obj.m1 + obj.m30));  % basal insulin in plasma
            obj.Ilb = obj.m2 / (obj.m1 + obj.m30) * obj.Ipb;
            obj.Ib = obj.Ipb / obj.VI;
            obj.kgri = obj.kmax;
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
                    warning('Property %s does not exist in patient_00.', fieldName);
                end
            end
        end


    end

end