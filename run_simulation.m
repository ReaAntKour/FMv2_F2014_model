function [] = run_simulation(resultFileName,tempDay,tempNight,tempDayLength,tempTwilight,co2Day,co2Night,co2DayLength,co2Twilight,lightDay,lightNight,lightDayLength,lightTwilight,genotype)
%% wrapper function to simulate the model
%
% Input:
%   resultFileName - file path to write results to
%   [temp,co2,light]Day - daytime [temp,co2,light]
%   [temp,co2,light]Night - nighttime [temp,co2,light]
%   [temp,co2,light]DayLength - duration of day [temp,co2,light]
%   [temp,co2,light]Twilight - duration of transition between day and night [temp,co2,light]
%   genotype - genotype to be simulated ('WT' for wildtype, else prr7prr9 is assumed)
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


  fileID = fopen(resultFileName,'w');

  %hour,T,sunrise,sunset,CO2,PAR
  load('weather.mat')

  hour=weather(:,1);

  rise = 1;
  set = rise+lightDayLength;

  T=issf(hour, tempDay, tempNight, tempDayLength, tempTwilight);
  sunrise=rise*weather(:,3);
  sunset=set*weather(:,4); 
  CO2=issf(hour, co2Day, co2Night, co2DayLength, co2Twilight);
  PAR=issf(hour, lightDay, lightNight, lightDayLength, lightTwilight);
  Photoperiod = sunset(end) - sunrise(end);

  %Specifying the genotype for the clock and starch models

  if strcmp(genotype,'WT')
    w = 0.91; %water content
    mf_use = 0.7; %malate+fumarate turnover
    clock_genotype = {''};
  else  % assume prr7prr9 double mutant
    w = 0.89; %water content
    mf_use = 0.25; %malate+fumarate turnover
    clock_genotype = {'prr9','prr7'};
  end

  d = 1-w; %dry matter

  load('parameter.mat')

  % Could specify mutants here, e.g.: clock_genotype = {'prr9','prr7'};
  clock_parameters = P2011_parameter_call(clock_genotype);

  starch_genotype = {''};
  starch_parameters = starch_parameter_call(starch_genotype);


  simulate_FM(hour,T,sunrise,sunset,CO2,PAR,Photoperiod,clock_parameters,starch_parameters,parameter,d,mf_use,true,fileID);

  % Close the result file
  fclose(fileID);

end
