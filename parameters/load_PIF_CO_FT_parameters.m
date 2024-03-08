function P = load_PIF_CO_FT_parameters(genotype,temperature)
% Code edited by Rea L Antoniou-Kourounioti (University of Glasgow) as indicated
%% returns a vector of parameter values for a given genotype (e.g. {"pif4", "pif5"} is the pif4pif5 double knockout)
% 
% Input:
%   genotype - structure of character strings specifying knockout/overexpression genes
%   temperature - either 22 or 27 (different EC parameters at the two temperatures)
%
% Output:
%   P - vector of Parameters
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

if nargin >2
    error('Too many input parameters to load_PIF_CO_FT_parameters')
elseif nargin > 1
    % Temperature has been input by the user
elseif nargin == 1
    % Default temperature value
    temperature = 22;
elseif nargin == 0
    % Default temperature and genotype
    temperature = 22;
    genotype = {''};
end

load('PIF_CO_FT_parameters','P')

% Modify parameters according to genotype
if ismember('co',genotype)
    P.n5 = 0;
    P.n4 = 0;
    P.Bco = 0;
end
if ismember('COox',genotype)
    P.n5 = 0;
    P.n4 = 0;
    P.Bco = 2;
end
if ismember('cdf1',genotype)
    P.n2 = 0;
    P.n1 = 0;
end
if ismember('CDF1ox',genotype)
    P.n2 = 0;
    P.n1 = 2;
    P.g2 = 100000;
end
if ismember('fkf1',genotype)
    P.q1 = 0;
    P.n3 = 0;
end
if ismember('pif4',genotype)
    P.n7 = 0;
    P.n6 = 0;
end
if ismember('pif5',genotype)
    P.n8 = 0;
    P.n9 = 0;
end
if ismember('PIF4ox',genotype)
    P.g7 = P.g7*10000;
    P.n6 = P.n6*2;
end
if ismember('PIF5ox',genotype)
    P.g8 = P.g8*10000;
    P.n8 = P.n8*2;
end
if ismember('delta1',genotype)
    % No CDF1 destabilisation by FKF1
    P.p2 = 0;
end
if ismember('delta2',genotype)
    % No CO stabilisation by FKF1
    %make michaelis constant very high to reduce FKF1 effect to negligible level
    P.k2 = 100000;
end
if ismember('cop1',genotype)
    % No role for cop1 in this pathway
    P.m7 = 0;
    P.n5 = 0;
end

%% Code edited here to add YHB mutant
if ismember('YHB',genotype)
    % constitutively active PhyB, not light-dependent
    P.YHB = 1; % does not depend on yhb parameter of clock model. Only phyB here, not all red light inputs.
else
	% wild-type has light-dependent PhyB activity
    P.YHB = 0;
end
%% Code edit ends here

% Modify parameters according to temperature
if temperature == 22;
    %do nothing
elseif temperature == 27;
    % Relieve EC inhibition of PIF4 transcription
    P.g7 = P.g7*4;
    % FT activation at the higher temperature
    f = 3.24;
    P.n14 = P.n14*f;
    P.n15 = P.n15*f;
else
    %raise exception
    error('Invalid temperature selected: Temperature must be 22oC or 27oC')
end
