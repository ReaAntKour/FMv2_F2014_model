function [output,sim_data] = simulate_FM_simp(hour,T,sunrise,sunset,clock_parameters,clock_dynamics,clock_dynamics_model_i,p,run_phenology_model,options)
% Code significantly edited by Rea L Antoniou-Kourounioti (University of
% Glasgow) to simplify and remove modules not used in this version, to add
% new outputs, and to include both circadian models. 
% Otherwise based on below: 
%% simulate the Framework Model v2. This function is the core function,
%
% Input:
%   hour - vector of timesteps
%   T - vector of temperatures over timesteps given by hour
%   sunrise - daily time of sunrise (integer, in hours)
%   sunset - daily time of sunset (integer, in hours)
%   clock_parameters - vector of P2011 circadian clock model parameters
%   p - vector of Framework Model parameters (beyond those used for the clock and starch models)
%   run_phenology_model - boolean variable determining whether to run the phenology model
%   fileID - file to write output data to
%
% Output:
%   output - subset of output data - fresh weight, assimilation, respiration, and starch content on day 38
%   sim_data - structure containing nested files with complete simulation data
%
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

%specify Col accession for phenology model threshold
flowering_thresh_geno = 2;

% Adding PIF model path
addpath('PIF_CO_FT_model')

% Initialise clock to starting conditions
clock_output=entrain_clock_model(clock_parameters,clock_dynamics,sunrise(1),sunrise(1)+options.entrain,options.y0);
%work out the clock state at ZT24 i.e. at the start of the next day
clock_state_0 = interp1q(clock_output.T,clock_output.Y,24);

FT_module_state =ones(1,18);

% set whether the phenology model controls the end of the simulation
N_max_days = 40; %default number of days max simulation
if run_phenology_model
    %upper ceiling on the number of days max simulation, if the phenology model is being run
    N_max_days = 90;
end


day_idx = 1;
has_flowered = false;
CumPhenThrm=0;
CumPhenThrmAll=zeros(90,1);

while day_idx <= N_max_days && ~(has_flowered)
    %initial timepoint
    t = (day_idx-1)*24+1;

    %run clock model for this day
	clock_output=run_clock_model_for_a_day(clock_parameters,clock_dynamics,sunrise(t),sunset(t),clock_state_0);
    %work out the clock state at ZT24 i.e. at the start of the next day
    clock_state_0 = interp1q(clock_output.T,clock_output.Y,24);
    
    %run phenology model
    if run_phenology_model
	    [DayPhenThrm,FT_module_state,PIF_output] = phen(T,t,sunrise(t),sunset(t),flowering_thresh_geno,clock_dynamics_model_i,clock_output,FT_module_state,p,options.genotype,options.hypocotyl_parameters);
    else
        DayPhenThrm = 0;
    end
    CumPhenThrm = DayPhenThrm+CumPhenThrm;
    has_flowered = flowering_threshold_test(CumPhenThrm,flowering_thresh_geno,p);
    CumPhenThrmAll(day_idx)=CumPhenThrm; % Keep record of CumPhenThrm over time

    %Update day index
    day_idx = day_idx+1;
end


%-------------------------------End-of-model-------------------------------

sim_data = struct();
output = [];

try
    %Try to instantiate a complete output as far as possible
    sim_data.metadata.hour = hour;
    sim_data.metadata.T = T;
    sim_data.metadata.sunrise = sunrise;
    sim_data.metadata.sunset = sunset;
    sim_data.metadata.Photoperiod = sunset - sunrise;
    sim_data.metadata.clock_parameters = clock_parameters;
    sim_data.metadata.p = p;
    sim_data.metadata.run_phenology_model = run_phenology_model;
    
	sim_data.clock_output = clock_output;
	sim_data.PIF_output = PIF_output;
	sim_data.Fl = day_idx - 1;
	sim_data.Hyp = PIF_output.Hyp_length;
    sim_data.CumPhenThrm = CumPhenThrmAll;
    output = [sim_data.clock,sim_data.PIF,sim_data.Hyp,sim_data.Fl];
catch
    %Do nothing
end
