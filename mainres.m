function [leaf_res,root_res] = mainres(Tleaf,leaf_c,root_c,suc_c_perplant,rosette_area,timestep,p)
%% calculates rate of maintenance respiration
%
% Inputs:
%    Tleaf - Leaf temperature (degrees Centigrade)
%    leaf_c - Leaf carbon (gC/plant)
%    root_c - Root carbon (gC/plant)
%    suc_c - Sucrose carbon
%    rosette_area - Rosette area (m2/plant)
%    timestep - 1 h by model-wide, hard-wired default
%    p - Vector of parameters
%
% Outputs:
%    leaf res - leaf respiration per plant per time step (gC per plant per
%    time step)
%    root_res - root respiration per plant per time step (gC per plant per
%    time step
%
%  Functions called:
%    -
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


%Calculation for maintenance respiration
%_______________________________________
%_______________________________________

act_ener_res20 = p(56); %activation energy for rl 20 (kJ/mol)

suc_conc_s = suc_c_perplant/rosette_area;  

rl20_leafres = (p(57)*suc_conc_s + p(58))*24; %Leaf respiration at 20 (g CO2 C/m2/day)

numres = act_ener_res20*(Tleaf - 20);
denomres = 293*p(39)*(Tleaf + 273);

if  suc_conc_s <= 0
    
    rl_leafres = 0;
else
    rl_leafres = rl20_leafres*exp(numres/denomres); %leaf respiration at leaf temperature (g CO2 C/m2/day)
end

leaf_res = rl_leafres*rosette_area*timestep; %leaf respiration per plant per time step (gC per plant per time step)

root_res = leaf_res*root_c/leaf_c; %gC/plant/time step
