classdef patient_11b
    properties
        % Patient features
        BW = 86.1; % body weight [kg]
        eat_rate = 5; % [g/min]
        
        % To be set by physician
        basal;
        u2ss; % steady state insulin rate (IIRb)
        % u2ss;      
        
        % To be measured at basal state
        Gb = 121.56;
        
        % Sensor features
        Td = 10; % glucose sensor delay
        
        % Glucose Kinetics
        VG = 1.5;
        k1 = 0.0852818455967453;
        k2 = 0.495649287405968;
        Gpb; % basal glucose in plasma
        % Insulin Kinetics
        VI = 0.05;
        HEb = 0.6;
        CL;
        m1 = 0.193036980628340;
        m2;
        m30;
        m4;
        m5 = 0.0304;
        Ipb;  % basal insulin in plasma
        Ilb;
        Ib;
        % Rate of Appearance
        kmax = 0.0586392094261030;
        kmin = 0.00747315367504439;
        kabs = 0.0548721699886073;
        kgri;
        % f = 0.9;
        f = 1;
        b = 0.82;
        d = 0.010;
        % Endogenous Glucose Production
        kp1 = 2.01025779577412;
        kp2 = 0.00158976840701995;
        kp3 = 0.000199488177267483;
        ki = 0.00481304394026455;
        EGPb;
        % Utilization
        Fcns = 1;
        Gtb;  % basal glucose in slowly equilibrating tissues
        Km0 = 225.59;
        Vm0 = 2.5;
        Vmx = 0.0912778113975294;
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
        function obj = patient_11b(dailyBasal)
            obj.basal = dailyBasal / 1440;
            obj.u2ss = obj.basal * 6000 / obj.BW;

            obj.m30 = obj.m1 * obj.HEb / (1 - obj.HEb);
            obj.CL = 0.0242 * obj.BW;
            obj.m2 = 3/5 * obj.CL / (obj.HEb * obj.VI * obj.BW);
            obj.m4 = 2/5 * obj.CL / (obj.VI * obj.BW);
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

        function obj = recompile(obj)

            obj.m30 = obj.m1 * obj.HEb / (1 - obj.HEb);
            obj.m2 = 3/5 * obj.CL / (obj.HEb * obj.VI * obj.BW);
            obj.m4 = 2/5 * obj.CL / (obj.VI * obj.BW);
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
                    warning('Property %s does not exist in patient_11b.', fieldName);
                end
            end

            obj = obj.recompile();

        end

    end

end
