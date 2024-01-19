function [suc_sta_base,suc_MF_base,sta_use,MF_use,suc_equi,al_suc,suc_c_disp,suc_c_interm,al_pt_plant_assim]...
    = assimilation_MF_clock(daylength,is_light,net_rate,timestep,leaf_res,...
    root_res,suc_c_perplant,rosette_area,sta_c_endday,MF_c_endday,leaf_c,starch_consumption,p,mf_use)
%% rates of assimilation and carbon partitioning between starch and sucrose
%
% Inputs:
%   daylength - Duration of daylight
%   is_light - Whether lights are on or off in the current timestep
%   net_rate - net_rate - Net rate of leaf photosynthesis (micromol CO2 m-2 leaf s-1)
%   timestep - 1 h by model-wide, hard-wired default
%   leaf res - leaf respiration per plant per time step (gC per plant per
%   time step)
%   root_res - root respiration per plant per time step (gC per plant per
%   time step
%   suc_c_perplant - Sucrose carbon (gC/plant)
%   rosette_area - Rosette area (m2/plant)
%   sta_c_endday - Starch carbon at the end of the previous day
%   MF_c_endday - MF carbon at the end of the previous day
%   leaf_c - Leaf carbon (gC/plant)
%   starch_consumption - Rate of starch consumption (gC/plant/timestep)
%   p - Vector of parameters
%   mf_use - Fraction of malate and fumarate used across the night
%
% Outputs:
%   suc_sta_base - baseline starch conversion
%   suc_MF_base - baseline malate+fumarate conversion
%   sta_use - conversion of starch to sugar
%   MF_use - conversion of malate+fumarate to sugar
%   suc_equi - equilibrium sucrose plus hexose concentration in leaves (g C/m2 leaf)
%   al_suc - allocation remaining for sucrose after accounting for
%   conversion to starch and m+f
%   suc_c_disp - Amount of sugar available for growth
%   suc_c_interm - Intermediate sucrose carbon calculation
%   al_pt_plant_assim - Assimilatory flux per plant (gC/plant/timestep)
%
% Called functions:
%   -
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

convert_to_gC = timestep*p(59)*10^(-6)*12; %conversion factor for umol/m2 leaf/sec to gC/m2 leaf/timestep

%Baseline conversion coefficient based on AGPase activity data
sta_base = p(60)*(-0.0296*daylength + 0.7157);

mf_syn_frac = p(85);

%Malate and fumarate conversion coefficient
MF_base = mf_syn_frac*sta_base;
  
MF_convert_night = mf_use;

if      is_light == 1
    
        sta_use = 0; %conversion of starch to sugar
        MF_use = 0; %conversion of fumarate and malate to sugar
        al_pt_plant_assim = net_rate*rosette_area*convert_to_gC; %Assimilatory flux per plant (gC/plant/timestep)
        suc_sta_base = al_pt_plant_assim*sta_base; %baseline starch conversion
        suc_MF_base = al_pt_plant_assim*MF_base; %baseline fm conversion
        al_suc =al_pt_plant_assim - suc_sta_base - suc_MF_base;
else
        al_pt_plant_assim = 0;
        sta_use = starch_consumption; %Dark conversion of starch to sugar.
        MF_use = ((MF_c_endday*MF_convert_night)/(24-daylength))*timestep*24;
        suc_sta_base = 0;
        suc_MF_base = 0;
        al_suc = 0;
end

suc_equi = p(62); %equilibrium sucrose plus hexose concentration in leaves (g C/m2 leaf)

rosette_DW = leaf_c/p(30); %rosette DW calculated from leaf C and carbon content
suc_c_interm = suc_c_perplant + sta_use + MF_use + al_suc - root_res - leaf_res;

current_value = suc_c_interm - (suc_equi*rosette_area);


if      current_value <= 0
        suc_c_disp = 0;
else
        suc_c_disp = current_value; %Amount of sugar available for growth
end
