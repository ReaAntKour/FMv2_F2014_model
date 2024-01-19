function [new_starch_module_state,starch_consumption] = starch_module(time,starch_module_state,sta_c,leaf_c,is_light,clock_output,starch_parameters)
% Input:
%   time - current timepoint
%   starch_module_state - initial state vector of the starch module in this timestep
%   sta_c - starch carbon
%   leaf_c - leaf carbon (used to normalise starch content to plant size)
%   is_light - boolean, presence or absence of light in this timestep
%   clock_output - dynamic vector of the circadian clock model, which provides input to the starch model
%   starch_parameters - model parameters for the starch model dynamics (kinetic parameters, binding constants, etc.)
%
% Output:
%   new_starch_module_state - state vector of the starch module at the end of this timestep
%   starch_consumption - the quantity of starch consumed in this timestep
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

% Note: Dynamics (dydt) are directly manipulated to convert units from per
% hour (this model) to per day (framework model).
% Note: Starch concentration is systematically rescaled to put it into the
% scale for which the model was originally parameterised

starch_factor = 60;

y0 = starch_module_state;
sta_c_conc = starch_factor*sta_c/leaf_c;
y0(3) = sta_c_conc;

tspan = [time-1,time];

% Run model for one hour:
[T,Y]=ode15s(@(t,y)P2011_starch_combined_dynamics(t,y,starch_parameters,is_light,clock_output),tspan,y0);

new_starch_module_state = Y(end,:)';
starch_consumption = -leaf_c*(Y(end,3) - y0(3))/starch_factor;

function dydt = P2011_starch_combined_dynamics(t,y,starch_parameters,is_light,clock_output)

P_starch = starch_parameters;

y_P2011 = interp1q(clock_output.T,clock_output.Y,t);
y_starch = y;

% Define light conditions:
L=is_light;

% Signalling from clock to starch by modulation of starch model parameters:
g = 2;
P_starch.ksT = P_starch.ksT2*(y_P2011(5)+y_P2011(20)).^g/(P_starch.KsT2.^g + (y_P2011(5)+y_P2011(20)).^g) + P_starch.ksT1;

% Submodel dynamics:
dy_starch = starch_dynamics(t,y_starch,P_starch,L);

% Combined derivative:
dydt = dy_starch;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model variant 2 dynamics
function Func = starch_dynamics(t,y,par,L)
    % This function only considers starch degradation
    X = y(1);
    Y = y(2);
    S = y(3);
    T = y(4);

    dXdt = (par.ksX*S - par.kdX*X*T);
    dYdt = 0;
    dSdt = - (1-L).*par.kdS*(X/(par.KdSX + X)) .* S ./ (S+par.KdS);
    dTdt = par.ksT*L  - par.kdT1*L*T - par.kdT2*T/(par.KdT2 + T);

    Func = [dXdt; dYdt; dSdt;dTdt];
end

end

end
