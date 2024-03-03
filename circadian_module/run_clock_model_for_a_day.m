function [clock_output]=run_clock_model_for_a_day(parameters,clock_dynamics,sunrise,sunset,y0)
	% Run model for 27 hours:
	Tout = 0:0.05:27;
	[T,Y] = ode15s(@(t,y) clock_dynamics(t,y,parameters,sunrise,sunset),Tout,y0);
	
	clock_output = struct();
	clock_output.T = T;
	clock_output.Y = Y;
end