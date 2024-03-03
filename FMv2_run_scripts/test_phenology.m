temp = 20.5;
rise = 0;
set = 12;
co2 = 42;
light = 145;

addpath('PIF_CO_FT_model')

geno = 2;

global p d mf_use

load('parameter.mat')

p=parameter;

load('weather.mat')

hour=weather(:,1);
T=temp*weather(:,2);
sunrise=rise*weather(:,3);
sunset=set*weather(:,4); 
CO2=co2*weather(:,5); %CO2 partial pressure (Pa)
PAR=light*weather(:,6); %total absorbed PAR per unit leaf area (micromol m-2 s-1)
Photoperiod = sunset(end) - sunrise(end);


clock_genotype = {''};
% clock_genotype = {'prr9','prr7'};

%initialise clock and FT module to starting conditions
FT_module_state = ones(1,18);
clock_state_0=[1.0151 0.956 0.0755 0.0041 0.506 0.0977 0.0238 0.0731 0.0697 0.0196 0.0435 0.2505 0.0709 0.1017 0.0658 0.4016 0.1167 0.1012 0.207 0.0788 0.3102 0.0553 0.2991 0.1503 0.0286 0.65 0.2566 0.1012 0.576 0.3269]; %12:12 wt;
for i = 1:5
    clock_output=circadian_module(sunrise(1),sunset(1),clock_state_0,clock_genotype);
    %work out the clock state at ZT24 i.e. at the end of the previous day
    clock_state_0 = interp1q(clock_output.T,clock_output.Y,24);
    [~,FT_module_state] = simulate_PIF_CO_FT_model(sunrise(1),sunset(1),clock_output,clock_genotype,FT_module_state);
end

starch_module_state = [1,1,1,1];

day_idx = 1;
N_max_days = 120; %100 days max simulation
has_flowered = false;
has_emerged = false;
CumPhenThrm=0;

while day_idx <= N_max_days && ~(has_flowered)
    %initial timepoint
    t = (day_idx-1)*24+1;
    
    %run clock model for this day
    clock_output=circadian_module(sunrise(t),sunset(t),clock_state_0,clock_genotype);
    %work out the clock state at ZT24 i.e. at the end of the previous day
    clock_state_0 = interp1q(clock_output.T,clock_output.Y,24);
    
    %run phenology model
    [DayPhenThrm,FT_module_state] = phen(hour,T,sunrise(t),sunset(t),geno,clock_genotype,clock_output,FT_module_state);
    CumPhenThrm = DayPhenThrm+CumPhenThrm;
    has_flowered = flowering_threshold_test(CumPhenThrm,geno);
    
    if ~(has_flowered)
        day_idx = day_idx + 1;
    end
end

[Bolting_point] = phenology (hour,T,sunrise,sunset,geno)
day_idx