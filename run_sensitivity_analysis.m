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

temp = 20.5;
rise = 0;
set = 12;
co2 = 42;
light = 145;

%Calling for meteorological data
%_______________________________  

%hour,T,sunrise,sunset,CO2,PAR
load('weather.mat')

hour=weather(:,1);
T=temp*weather(:,2);
sunrise=rise*weather(:,3);
sunset=set*weather(:,4); 
CO2=co2*weather(:,5); %CO2 partial pressure (Pa)
PAR=light*weather(:,6); %total absorbed PAR per unit leaf area (micromol m-2 s-1)
Photoperiod = sunset(end) - sunrise(end);

output_filename = input(['Enter name for datafile output \n(e.g. sens_analysis_results_WT):'], 's');


%WT water content and MF
w = 0.91; %water content
mf_use = 0.7; %malate+fumarate turnover
d = 1-w; %dry matter

load('parameter.mat')

p=parameter;
p0 = p;

deltaP = 0.01;

nP = length(p);
% nP = 2;

clock_genotype = {''};
% clock_genotype = {'prr9','prr7'};
starch_genotype = {''};
clock_parameters0 = P2011_parameter_call(clock_genotype);
clock_parameters = clock_parameters0;
nCP = length(clock_parameters);

starch_parameters0 = starch_parameter_call(starch_genotype);
starch_parameters = starch_parameters0;
starch_parameter_names = fieldnames(starch_parameters);
nSP = length(starch_parameter_names);

run_phenology_model = 0;

[output_basal,~] = simulate_FM(hour,T,sunrise,sunset,CO2,PAR,Photoperiod,clock_parameters,starch_parameters,p0,d,mf_use,run_phenology_model);

nD = length(output_basal); %number of dimensions of output

%Add two additional parameters for water content and mf_use

sens = zeros(2+nP+nCP+nSP,nD);
errors = [];

p = p0;
%water content
w_perturbed = w*(1+deltaP);
d_perturbed = 1-w_perturbed;
[output,~] = simulate_FM(hour,T,sunrise,sunset,CO2,PAR,Photoperiod,clock_parameters,starch_parameters,p,d_perturbed,mf_use,run_phenology_model);
sens(1,:) = ((output-output_basal)/deltaP)./output_basal;

mf_use_perturbed = mf_use*(1+deltaP);
[output,~] = simulate_FM(hour,T,sunrise,sunset,CO2,PAR,Photoperiod,clock_parameters,starch_parameters,p,d,mf_use_perturbed,run_phenology_model);
sens(2,:) = ((output-output_basal)/deltaP)./output_basal;


for i = 3:nP+2
    i
    p = p0;
    p(i-2) = p(i-2)*(1+deltaP);
    try
        [output,~] = simulate_FM(hour,T,sunrise,sunset,CO2,PAR,Photoperiod,clock_parameters,starch_parameters,p,d,mf_use,run_phenology_model);
        sens(i,:) = ((output-output_basal)/deltaP)./output_basal;
    catch
        errors = [errors;i];
    end
end

p = p0;
for i = 2+nP+1:2+nP+nCP
    i
    clock_parameters = clock_parameters0;
    clock_parameters(i-nP-2) = clock_parameters(i-nP-2)*(1+deltaP);
    try
        [output,~] = simulate_FM(hour,T,sunrise,sunset,CO2,PAR,Photoperiod,clock_parameters,starch_parameters,p,d,mf_use,run_phenology_model);
        sens(i,:) = ((output-output_basal)/deltaP)./output_basal;
    catch
        errors = [errors;i];
    end
end

p = p0;
clock_parameters = clock_parameters0;
for i = 2+nP+nCP+1:2+nP+nCP+nSP
    i
    starch_parameters = starch_parameters0;
    starch_parameters.(starch_parameter_names{i-nP-nCP-2}) = starch_parameters.(starch_parameter_names{i-nP-nCP-2})*(1+deltaP);
    try
        [output,~] = simulate_FM(hour,T,sunrise,sunset,CO2,PAR,Photoperiod,clock_parameters,starch_parameters,p,d,mf_use,run_phenology_model);
        sens(i,:) = ((output-output_basal)/deltaP)./output_basal;
    catch
        errors = [errors;i];
    end
end

save(output_filename,'output_basal','temp','rise','set','co2','light','p0','clock_parameters0','clock_genotype','deltaP','errors','sens')
