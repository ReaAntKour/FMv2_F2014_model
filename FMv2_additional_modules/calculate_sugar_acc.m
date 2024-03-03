function [efficiency,growth_capacity]=calculate_sugar_acc(t,p,eme,efficiency,growth_capacity,is_light)
%% calculate sugar accumulation
%
% Input:
%   t - current time
%   p - vector of model parameters
%   eme - time of emergence
%   efficiency - timeseries vector of photosynthetic efficiency
%   growth_capacity - timeseries vector of growth capacity
%   is_light - timeseries vector specifying when there is light
%
% Output:
%   efficiency - timeseries vector of photosynthetic efficiency
%   growth_capacity - timeseries vector of growth capacity
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

if t < eme + 24 + 1
    efficiency(t) = p(84);
    growth_capacity(t) = p(63);
else

    growth_capacity(t) = p(63);


    if is_light(t) == 1 && is_light(t-1) == 0

        efficiency(t) = p(84);
    else
        efficiency(t) = efficiency(t-1);
    end
end    
