function [ data ] = issf( time, day, night, dayLength, twilight)
% Calculates 24 hour periodic data uses to simulate environmental inputs
%   Specification taken from http://jbr.sagepub.com/content/27/4/328.full
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

  
    % Move the center of the transition down by 0.5 hours to allow going from high to low in one hour
    % Also remove 1 because the simulation starts at time 1 but modular arithmetic has base 0 
    tPrime = mod(time-1+0.5,24);
    theta0 = night;
    theta1 = day - night;
    Tp = dayLength;
    Tc = 24;
    T = max(0.0001,twilight);
    
    data = theta0 + 0.5 * theta1 * ( ( 1 + tanh(tPrime/T)) - (1+ tanh((tPrime-Tp)/T)) + (1 + tanh((tPrime-Tc)/T)));

end
