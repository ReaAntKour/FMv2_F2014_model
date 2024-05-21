function [dailyFTarea,FT_module_state,PIF_output] = simulate_PIF_CO_FT_model(sunrise, sunset, clock_output, clock_dynamics_model_i, FT_module_state, temperature, options) % Code edited here to add new inputs and outputs
% Code edited from original code to simulate the Arabidopsis Framework
% model v2, described in Chew et al, 2017 [https://doi.org/10.1101/105437]. 
% Copyright 2018 Yin Hoon Chew, Daniel Seaton, Andrew Millar, and The University of Edinburgh
% by Rea L Antoniou-Kourounioti (University of Glasgow) as indicated
%%
assert(length(sunrise)==1)
assert(length(sunset)==1)

% Code edited here to remove temperature and flowering genotype information. 
% These are now inputs to the function. Temperature may only take values of 22 (default parameter settings) or 27

% Load light conditions into 'c' for common light function
c.period = 24;
c.phase = 0;
c.dawn = sunrise;
c.photoperiod = sunset-sunrise;

% % Include model folders in path
% addpath('P2011_model')
% addpath('PIF_CO_FT_model')

%% Code edited here to allow different circadian model versions to be used as input
% extract parameters from options
flowering_genotype = options.genotype;
hp = options.hypocotyl_parameters;
paramSet = options.paramSet;
YHB = options.YHB; % 0.5 (33%)Y, 1 (50%), 3 (75%)Y or 9 (90%)

% Convert clock dynamics into input struct for PIF_CO_FT model
if clock_dynamics_model_i<3
	clock_parameters = load_P2011_parameters(flowering_genotype,YHB);
	u = wrap_P2011_model_dynamics(clock_output.T,clock_output.Y,clock_parameters);
else
	clock_parameters = load_F2014_parameters(flowering_genotype,paramSet,YHB);
	u = wrap_F2014_model_dynamics(clock_output.T,clock_output.Y,clock_parameters);
end
%% Code edit ends here

% Load PIF_CO_FT parameters
parameters.PIF_CO_FT = load_PIF_CO_FT_parameters(flowering_genotype,temperature);


% Run CO-PIF-FT model
% Simulate for one day
[T,Y] = ode15s(@(t,y) PIF_CO_FT_dynamics(t,y,parameters.PIF_CO_FT,u,c),[0 c.period],FT_module_state);
FT_module_state = Y(end,:)';


dailyFTarea = trapz(T,Y(:,15));


%% Code edited here to add hypocotyl length calculation and function output
% calculate area under ATHB2m as proxy for hypocotyl length
Hyp_length=hp.a1*trapz(T,min(Y(:,8),hp.a3)-hp.a2);

PIF_output = struct();
PIF_output.T = T;
PIF_output.Y = Y;
PIF_output.FTarea = dailyFTarea;
PIF_output.Hyp_length = Hyp_length;