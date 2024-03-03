function [clock_output]=entrain_clock_model(parameters,clock_dynamics,sunrise,sunset,y0)
	% Initialise clock model (12 days)
	[~,Ycinit] = ode15s(@(t,y) clock_dynamics(t,y,parameters,sunrise,sunset),0:0.1:(12*24),y0);
	y0 = Ycinit(end,:)';
	
	% Run model for 27 hours:
	Tout = 0:0.05:27;
	[T,Y] = ode15s(@(t,y) clock_dynamics(t,y,parameters,sunrise,sunset),Tout,y0);
	
	clock_output = struct();
	clock_output.T = T;
	clock_output.Y = Y;
end