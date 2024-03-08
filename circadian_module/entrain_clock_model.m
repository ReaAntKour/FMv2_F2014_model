function [clock_output]=entrain_clock_model(parameters,clock_dynamics,sunrise,sunset,y0)
% 
%   Copyright 2024 Rea L Antoniou-Kourounioti and The University of Glasgow
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
%
%       http://www.apache.org/licenses/LICENSE-2.0
%
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.
% 
% 
% File is based on the original FMv2 file "circadian_module.m", Copyright
% 2018 Yin Hoon Chew, Daniel Seaton, Andrew Millar, and The University of
% Edinburgh, licensed under the Apache License, Version 2.0


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