function [dailyFTarea,FT_module_state,PIF_output] = simulate_PIF_CO_FT_model(sunrise, sunset, clock_output, clock_dynamics_model_i, FT_module_state, temperature, flowering_genotype, hp) % Code edited here and many places
assert(length(sunrise)==1)
assert(length(sunset)==1)

% Load light conditions into 'c' for common light function
c.period = 24;
c.phase = 0;
c.dawn = sunrise;
c.photoperiod = sunset-sunrise;

% % Include model folders in path
% addpath('P2011_model')
% addpath('PIF_CO_FT_model')


% Convert clock dynamics into input struct for PIF_CO_FT model
if clock_dynamics_model_i==1
	clock_parameters = load_P2011_parameters(flowering_genotype);
	u = wrap_P2011_model_dynamics(clock_output.T,clock_output.Y,clock_parameters); %%%%%%%%%%%%%%%%%%%%%%% Only call to wrap_P2011 function
elseif clock_dynamics_model_i==2
	clock_parameters = load_F2014_parameters(flowering_genotype);
	u = wrap_F2014_model_dynamics(clock_output.T,clock_output.Y,clock_parameters); %%%%%%%%%%%%%%%%%%%%%%% Only call to wrap_F2014 function
end

% Load PIF_CO_FT parameters
parameters.PIF_CO_FT = load_PIF_CO_FT_parameters(flowering_genotype,temperature);% temperature may only take values of 22 (default parameter settings) or 27


% Run CO-PIF-FT model
% Simulate for one day
[T,Y] = ode15s(@(t,y) PIF_CO_FT_dynamics(t,y,parameters.PIF_CO_FT,u,c),[0 c.period],FT_module_state);
FT_module_state = Y(end,:)';

dailyFTarea = trapz(T,Y(:,15));


% calculate area under ATHB2m as proxy for hypocotyl length
Hyp_length=hp.a1*trapz(T,min(Y(:,8),hp.a3)-hp.a2);

PIF_output = struct();
PIF_output.T = T;
PIF_output.Y = Y;
PIF_output.Hyp_length = Hyp_length;