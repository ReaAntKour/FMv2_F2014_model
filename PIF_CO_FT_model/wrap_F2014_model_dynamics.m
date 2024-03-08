function u = wrap_F2014_model_dynamics(T,Y,P)
% Code edited by Rea L Antoniou-Kourounioti (University of Glasgow) as indicated
%%
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

u.T = T;
nT = length(T);
varnames = {'LHY','EC','PRR9','PRR7','PRR5','TOC1','cP','COP1n_n','GIn'};
nV = length(varnames);
for i = 1:nV
    u.(varnames{i}) = zeros(nT,1);
end


%% Code edited to reflect new model beyond this point
% Code:
%   Copyright 2024 Rea L Antoniou-Kourounioti and The University of Glasgow
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

% Define Parameters
param = num2cell(P);
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,...
	~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,...
	~,~,~,~,~,~,~,~,~,f1,f2,f3,f4,~,f6,~,~,~,~,...
	~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,...
	~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,...
	~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,...
	~,~,~,~,~,~,~,~,~,~,~,~] = deal(param{:});
    

% Define the variables
LUXp = Y(:,23);
NOXp = Y(:,33);
ELF34 = Y(:,21);
ELF3p = Y(:,20);
u.EC = ((LUXp + f6 * NOXp) .* (ELF34 + f1 * ELF3p)) ./ (1 + f3 * (LUXp + f2 * NOXp) + f4 * (ELF34 + f1 * ELF3p));

u.COP1n_n = Y(:,25);
u.LHY = (Y(:,2) + Y(:,4))/1.561; % LHY + CCA1
u.PRR9 = Y(:,7);
u.PRR7 = Y(:,9)/2.6754;
u.cP = Y(:,5);
u.GIn = Y(:,31)*40.9;
u.PRR5 = Y(:,12)*0.841;% nuclear
u.TOC1 = Y(:,14)*1.21;% nuclear


