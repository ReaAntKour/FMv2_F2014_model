function [clock_output] = F2014_circadian_module(sunrise,sunset,clock_state_0,P) % Code editied here to change name
%% ODE model of the circadian clock, from Fogelmark & Troein (2014)
%
% Inputs:
%   sunrise
%   sunset
%   clock_state_0 - vector of molecular clock state at start of the day (t=0, relative to when sunrise and sunset occur)
%   P - vector of clock parameters (reaction rates, binding constants, etc. - not to be confused with parameters used for overall FM)
%
% Output:
%   clock_output - vector of molecular clock state at the end of the day
%
%   Copyright 2018 Yin Hoon Chew, Daniel Seaton, Andrew Millar, and The University of Edinburgh
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

 
param=load_F2014_parameters();% Code edited here to load parameters from new file

Tout = [0:0.05:27];
% Run model for 27 hours:
[T,Y]=ode15s(@(t,y) F2014_dynamics(t,y,param,sunrise,sunset),Tout,clock_state_0);% code edited here to update to new model

clock_output = struct();
clock_output.T = T;
clock_output.Y = Y;

% Code edited here to remove old model