state_fields = {'Qsto1','Qsto2','Qgut','Gp','Gt','Gpd','Il','Ip','I1','Id','X','Isc1','Isc2'};
extra_state_fields = {'insulin_to_infuse','last_IIR','CHO_to_eat','D','lastQsto','is_eating'};
input_fields = {'CHO', 'IIR'};
true_input_fields = {'CHO_consumed_rate','IIR_dt'};

% QSTO1 = 1;
% QSTO2 = 2;
% QGUT = 3;
% GP = 4;
% GT = 5;
% GSC = 6;
% IL = 7;
% IP = 8;
% I1 = 9;
% ID = 10;
% X = 11;
% ISC1 = 12;
% ISC2 = 13;
% 
% INS_TO_INF = 1;
% LAST_IIR = 2;
% CHO_TO_EAT = 3;
% D = 4;
% LASTQSTO = 5;
% IS_EATING = 6;
% 
% CHO = 1;
% IIR = 2;
% 
% CHO_RATE = 1;
% IIR_DT = 2;
