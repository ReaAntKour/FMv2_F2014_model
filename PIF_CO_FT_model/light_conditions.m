function L = light_conditions(t,c)
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

% Code edited here to move twilightPer definition below

LightOffset = 0; %Shifts light function up or down.
twilightPer=0.00005; %The duration of time between value of force in dark and value of force in light.% Code edited here to make twilight preriod shorter.
LightAmp = 1; %The amplitude of the light wave.

t = mod(t-15*twilightPer, c.period); % Code edited here to start twilight at ZT=0

if (c.photoperiod == 0)
    L = 0*t;
elseif (c.photoperiod == c.period)
    L = 1+0*t;
else
    L = LightOffset + 0.5*LightAmp*(1 + tanh((c.period/twilightPer)*((t+c.dawn)/c.period - floor(floor(t+c.dawn)/c.period)))) - 0.5*LightAmp*(1 + tanh((c.period/twilightPer)*((t+c.dawn)/c.period - floor(floor(t+c.dawn)/c.period)) - c.photoperiod/twilightPer)) + 0.5*LightAmp*(1 + tanh((c.period/twilightPer)*((t+c.dawn)/c.period - floor(floor(t+c.dawn)/c.period)) - c.period/twilightPer));
end
