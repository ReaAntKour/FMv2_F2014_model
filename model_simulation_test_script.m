%% interactive script to run the model for different cases

temp = 20.5;
rise = 0;
set = 12;
co2 = 42;
light = 145;

%Calling for meteorological data
%_______________________________  

%hour,T,sunrise,sunset,CO2,PAR
load('weather.mat')

hour=weather(:,1);
T=temp*weather(:,2);
sunrise=rise*weather(:,3);
sunset=set*weather(:,4); 
CO2=co2*weather(:,5); %CO2 partial pressure (Pa)
PAR=light*weather(:,6); %total absorbed PAR per unit leaf area (micromol m-2 s-1)
Photoperiod = sunset(end) - sunrise(end);



%Specifying the genotype for the clock and starch models
experiment_genotype_index = 1;
switch experiment_genotype_index
    case 1
        w = 0.91; %water content
        mf_use = 0.7; %malate+fumarate turnover
    case 2
        w = 0.89; %water content
        mf_use = 0.7; %malate+fumarate turnover
    case 3
        w = 0.89; %water content
        mf_use = 0.7; %malate+fumarate turnover
    case 4
        w = 0.89; %water content
        mf_use = 0.25; %malate+fumarate turnover
end

% w = 0.9;
d = 1-w; %dry matter


load('parameter.mat')

p=parameter;

clock_genotype = {''};
% clock_genotype = {'prr9','prr7'};
starch_genotype = {''};
% starch_genotype = {'lsf1'};

clock_parameters = P2011_parameter_call(clock_genotype);
nCP = length(clock_parameters);

starch_parameters = starch_parameter_call(starch_genotype);
nSP = length(starch_parameters);
% starch_parameters.turnover_override = true;
% starch_parameters.turnover = 0.90;

run_phenology_model = 1; % Code edited to run phenology model (changed from 0 to 1)

% [output,sim_data] = simulate_FM(hour,T,sunrise,sunset,CO2,PAR,Photoperiod,clock_parameters,starch_parameters,p,d,mf_use,run_phenology_model,clock_dynamics_model_i); % Code edited here to add "clock_dynamics_model_i"
simulate_FM_simp(hour,T,sunrise,sunset,clock_parameters,clock_dynamics,clock_dynamics_model_i,p,run_phenology_model);
