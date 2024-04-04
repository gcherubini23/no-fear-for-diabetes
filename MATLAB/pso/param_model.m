%% Comparing to GT
if ~use_true_patient
    % % Best using MSE
    % p = [0.00542767432046846,	0.406509811384502,	0.500000000000000,	4.48796062636850,	0.00480323418050043,	0.000104733132703769,	0.0934117134735254,	0.00934850839462120,	0.0977611506438287,	0.0265190704919319,	0.0352422437487089]; % BEST
    % p = [0.00692969041535363,	0.221175855482405,	0.289354129999467,
    % 4.88915947935168,	0.00443847724112982,	0.000100000000000000,
    % 0.0783288398117539,	0.00758782832475949,	0.0990034968223399,
    % 0.0293280078067010,	0.0313402517781819]; % BEST using longer experiment
    
    % % Best using RMSE
    % p = [0.00867171835016086,	0.394605478835008,	0.499442578412629,	5.95546142920697,	0.00428850132999850,	0.000913071382032137,	0.0997867581817115,	0.00703582396230473,	0.0961175442588771,	0.0386044391603998,	0.0224839553944235];
    
    %% Comparing to CGM
    % Best using MSE
    % p = [0.000790653343271446,	0.405640265724335,	0.490795840246288,	5.42022112811144,	0.00436200018973209,	0.000828901354548169,	0.100000000000000,	0.00796998392872804,	0.0950721054238413,	0.0270949155469932,	0.0317580189853282]; % BEST
    
    % Best using logMSE
    % p = [0.00146357466634368,	0.218798034701857,	0.332790084881620,	4.46551510843811,	0.00539047403268261,	0.000107668440747057,	0.0203292704426594,	0.00580387922110643,	0.100000000000000,	0.0290688859079165,	0.0319965142022406];
    
    % Best using log-likelihood
    % p = [0.00701668087717302,	0.369162325140569,	0.484425825956220,	4.88463681499870,	0.00440027384718579,	0.000188360321539952,	0.0928464063965499,	0.00767504547844108,	0.0876475567984728,	0.0293042078281975,	0.0313281318108638]; % BEST
    
    % p = [0.00716998351386864,	0.390062188468490,	0.482703397483656,	4.83804721507048,	0.00452135095508243,	0.000899234025720681,	0.0955592571364963,	0.00807954101149883,	0.100000000000000,	0.0300281789239006,	0.0314811981490732];
    % 
    % params = params.set_params({'kp2','k1','k2','kp1','ki','ke1','kmax','kmin','kabs','kp3','Vmx'}, p);
    
    
    %% If also Gb is estimated

    p = [0.00627984690700475,	0.273317762642556,	0.370499484862642,	4.64478309221182,	0.00451523454939320,	0.000328020658339412,	0.0886392094261030,	0.00747315367504439,	0.0948721699886073,	0.0276480508914665,	0.0323549128391175,	136.359246268003];
    
    params = params.set_params({'kp2','k1','k2','kp1','ki','ke1','kmax','kmin','kabs','kp3','Vmx','Gb'}, p);

end

%% True patient

if use_true_patient

    % p = [0.000171310660656453,	0.0445522261577045,	0.235559622154112,	3.55348628855907,	0.00994643442661664,	0.000174891811102616,	0.00496210677784614,	0.0100000000000000,	0.0840336530309484,	0.000100000000000000,	0.00159107213337244,	121.486862235716];
    
    % p = [0.0221123740335257,	0.0544318151208622,	0.409294347812267,	5.57537740205567,	0.000104418016028150,	0.000484289413385518,	0.0950786941768456,	0.00634647467555318,	0.0353369891388144,	0.00836496141759953,	0.00384529576640333,	102.936057045932];
    
    % p = [0.0221123740335257,	0.0544318151208622,	0.409294347812267,	5.57537740205567,	0.4804418016028150,	0.000484289413385518,	0.0558,	0.0080,	0.00353369891388144,	0.00836496141759953,	0.00384529576640333,	102.936057045932];
    
    p = [0.0140623315440670,	0.0247406166060108,	0.160203457794628,	2.75052935423773,	0.00428433560924936,	0.000275243603208159,	0.0100000000000000,	0.00459064722470484,	0.138063373785979,	0.00289568157586046, 0.00192591185929972, 64.8920869007200];
    
    params = params.set_params({'kp2','k1','k2','kp1','ki','ke1','kmax','kmin','kabs','kp3','Vmx','Gb'}, p);
    
end