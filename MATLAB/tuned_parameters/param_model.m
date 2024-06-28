%% Comparing to GT
if ~use_true_patient
    p = [0.00862715370302021	0.0869838631707173	0.468280280606665	5.60279336220474	0.00417337033678547	0.00274618001524637	0.0312571412812858	0.00600973432471946	0.241557076584449	0.0384071080739675	0.0614822407901443];
    params = params.set_params({'kp2','k1','k2','kp1','ki','ke1','kmax','kmin','kabs','kp3','Vmx'}, p);
end

%% True patient 11

if use_true_patient

    if use_anderson

        if patient_ID == 11
        p = [0.000369257569682681	0.0436103631924742	0.473351358459643	1.43174530155256	0.0473727284641901	9.76543232258312e-05	0.0361346005353632	0.00565322966894751	0.00950053481341468	0.00153454251792645	0.262451413035015];
        params = params.set_params({'kp2','k1','k2','kp1','ki','ke1','kmax','kmin','kabs','kp3','Vmx'}, p);

        end

        if patient_ID == 17
            % p = [1.99999675759867	0.397817947807300	1.50000000000000	0.0999416778642028	0.0100004785726379	300	0.199994301068677	0.00251429780584819	0.0133863719548739	0.00175663037962863	0.300000000000000];
            % p = [1.99997604566632	0.386300458054730	1.50000000000000	0.0591939233400531	0.0100044504129369	300	0.199969316332053	0.000100000000000000	0.0154578198347359	0.00174612055042329	0.300000000000000];
            % params = params.set_params({'VG','m1','CL','Vmx','k1','Km0','k2','kp2','kmax','kmin','kabs'}, p);
        


            % p = [0.0463151386139017	0.00850370689850918	0.00217204346863193	6	0.0124838977320622	1.00000000000000e-05	0.131405836603298	0.201301739455772	0.00364350133515147	1.00000000000000e-05	0.161006468573114];
            % params = params.set_params({'kp2','k1','k2','kp1','ke1','kmax','kmin', 'kgri','kabs','kp3','Vmx'}, p);



            p = [0.00499027344446903	4.77351115092225	0.00678620955008476	1.00000000000000e-05	0.0329961574486569	1.00000000000000e-05];
            params = params.set_params({'kp2','kp1','ki','ke1','kp3','Vmx'}, p);
        
        end



    end

    if use_tmoore
        % p = [2.65565440901474	0.649430243127675	3	0.00918951553502747	0.190116930370657	100	0.0250000000000000	0.00708837537237668	0.115296002949592	0.0200000000000000	0.0391769477649835	0.00200000000000000];
        % params = params.set_params({'VG','m1','CL','Vmx','k1','Km0','k2','kp2','kmax','kmin','kabs','ki'}, p);
     
        p = [1.50000000000000	0.399808524339688	1.50000000000000	0.0999913629537717	0.0146629675849774	200.003945291855	0.199668784611321	0.0100000000000000	0.00173116215659319	0.0100000000000000	0.300000000000000	0.00400000000000000	0.100000000000000];
        params = params.set_params({'VG','m1','CL','Vmx','k1','Km0','k2','kp2','kmax','kmin','kabs','ki','u2ss'}, p);
    
    end

    if use_shanghai
        if patient_ID == 1007
            p = [0.000165177452227904	0.0272022659129592	0.014427888711385	2.71420308225136	0.0104530747227785	1	0.0567647966896068	0.130379185678025	0.0261874493041442	0.0262491594641587	1.00000000000000e-05	1.93586988140253];
            params = params.set_params({'kp2','k1','k2','kp1','ki','ke1','kmax','kmin','kabs','kp3','Vmx','VG'}, p);

            % p = [0.0286216133026989	0.00821610437289417	0.000178806562253576	6	0.000287616376199055	0.281305786740562	0.00225241566656545	1	0.0393831522310934	0.0864613072938534	0.227259141040084];
            % params = params.set_params({'kp2','k1','k2','kp1','ke1','kmax','kmin', 'kgri','kabs','kp3','Vmx'}, p);

            % p = [0.00454426005719953	0.00154281342913677	1	1.72621557678894	0.00504531149754731	0.000500000000000000	0.0457286914335321	0.000373505839527540	0.00300167203992036	0.0106301750260069	0.00127114565303802];
            % params = params.set_params({'kp2','k1','k2','kp1','ke1','kmax','kmin', 'kgri','kabs','kp3','Vmx'}, p);
            % 
            % p = [1.00000000000000e-05	2.28306200392377	0.00843164022062208	0.000500000000000000	0.163995429750905	0.210223703553112	0.00290254495338659	0.00672318515678398	1.00000000000000e-05];
            % params = params.set_params({'kp2','kp1','ke1','kmax','kmin', 'kgri','kabs','kp3','Vmx'}, p);


            % p = [0.0160401217083955	2.75283000510846	1.00000000000000e-05	0.000500000000000000	1.00000000000000e-05	1.00000000000000e-05];
            % params = params.set_params({'kp2','kp1','ki','ke1','kp3','Vmx'}, p);
        end

        if patient_ID == 1002
            % p = [0.0176553540787431	0.148760818269538	1	6	0.00721792557853543	1	0.125241211983798	0.0196931318796094	0.0252920489257595	0.0187552276526636	0.557010182947336	2.84711585647146];
            % params = params.set_params({'kp2','k1','k2','kp1','ki','ke1','kmax','kmin','kabs','kp3','Vmx','VG'}, p);
        
            p = [1.00000000000000e-05	6	0.00683702629198257	0.190076978580308	0.0190981901870210	1];
            params = params.set_params({'kp2','kp1','ki','ke1','kp3','Vmx'}, p);
        
        end

        if patient_ID == 1010
            p = [0.0469407068907540	0.235472471650335	1	2.69917883678794	1.00000000000000e-05	1	0.498307026956838	0.911108818575275	1.00000000000000e-05	1.00000000000000e-05	0.162261765274725];
            params = params.set_params({'kp2','k1','k2','kp1','ki','ke1','kmax','kmin','kabs','kp3','Vmx'}, p);
        end

    end
end



