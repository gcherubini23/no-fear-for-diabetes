classdef utils
    methods          
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
        
        function [z] = sense(filename)
        
        end

        function [CHO] = announce_meal(filename)

        end

        function [IIR] = insulin_injection(filename)

        end

    end
end

