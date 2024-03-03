function has_flowered = flowering_threshold_test(CumPhenThrm,geno,p)
%% determine whether the flowering threshold has been passed
%
% Input:
%   CumPhenThrm - cumulative thermal units, adjusted by the phenology model depending on the photoperiod
%   geno - genotype specified. Different values result in different thresholds for flowering (for Col or Ler)
%   p - vector of FM parameters
%
% Output:
%   has_flowered - boolean, indicating whether the flowering threshold has been passed
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


if  geno ==1
    
    Threshold = p(13);%2907;%3259
    
else
    Threshold = p(17);%3212;%4218

end

if CumPhenThrm > Threshold
    has_flowered = true;
else
    has_flowered = false;
end
