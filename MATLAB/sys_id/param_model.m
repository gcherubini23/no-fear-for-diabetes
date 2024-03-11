% % Best using mse
% p = [0.00613673067477461,	0.400287164971832,	0.494617050929902,	4.92072844013362,	0.00382498896706389,	0.000826933211308203,	0.0998386175768517,	0.00854364392451055,	0.0997034961579507,	0.0309196769142360,	0.0389878346812570];

% % Best using rmse
% p = [0.00867171835016086,	0.394605478835008,	0.499442578412629,	5.95546142920697,	0.00428850132999850,	0.000913071382032137,	0.0997867581817115,	0.00703582396230473,	0.0961175442588771,	0.0386044391603998,	0.0224839553944235];

% Best using iterative mse
% p = [0.00832754941533622,	0.293957585844301,	0.363002205725996,	5.67368376491201,	0.00491783516666805,	0.000998846078555672,	0.0969893309513250,	0.00832285939257755,	0.0950049642840006,	0.0339534702044607,	0.0299321732306261];

% p = [0.00917483832670998,	0.0384422693118361,	0.0490699549522641,	5.99348314075621,	0.00383190671798864,	0.000335050490729362,	0.0272088100868198,	0.00894939293412607,	0.100000000000000,	0.0372923137615464,	0.0294514500733146];
% p = [0.00542767432046846,	0.406509811384502,	0.500000000000000,	4.48796062636850,	0.00480323418050043,	0.000104733132703769,	0.0934117134735254,	0.00934850839462120,	0.0977611506438287,	0.0265190704919319,	0.0352422437487089];


% Best
p = [0.00478484624671444,	0.379141730715929,	0.484131364169999,	5.68420836558207,	0.00471459706318235,	0.000240077378575619,	0.0950974847790787,	0.00603144000958310,	0.0675899945978740,	0.0342488487071233,	0.0249178888461712];

params = params.set_params({'kp2','k1','k2','kp1','ki','ke1','kmax','kmin','kabs','kp3','Vmx'}, p);
