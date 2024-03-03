function [rlc_pt1,suc_sta,root_gro1,rrc_pt1,leaf_gro1] = allocation(rrc_pt,totdem,rlc_pt,root_growth,leaf_growth,suc_c_disp,is_light)
%% allocation of carbon to growth and respiration
% Input:
%   rrc_pt - Root growth respiration In g C/plant/timestep
%   totdem - sum of all growth demand, including respiration
%   rlc_pt - Leaf growth respiration In g C/plant/timestep (required)
%   root_growth - Root growth in gC/plant/timestep
%   leaf_growth - Leaf growth in gC/plant/timestep
%   suc_c_disp - Amount of sugar available for growth
%   is_light - binary, whether it is light or dark
%   
%
% Output:
%   suc_growth - sucrose for growth (g C/plant)
%   rrc_ptl - Root growth respiration In g C/plant/timestep (actually used)
%   root_gro1 - root growth actually achieved
%   rlc_ptl - Leaf growth respiration In g C/plant/timestep (actually used)
%   leaf_gro1 - leaf growth actually achieved
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
%%

if      totdem < suc_c_disp
    
        suc_growth = totdem;
        rrc_pt1 = rrc_pt;
        root_gro1 = root_growth;
        rlc_pt1 = rlc_pt;
        leaf_gro1 = leaf_growth;
        
        if  is_light == 1
            suc_sta = suc_c_disp - suc_growth;
        else
            suc_sta = 0;
        end    

else
        suc_growth = suc_c_disp; %sucrose for growth (g C/plant)
        rrc_pt1 = rrc_pt*(suc_growth/totdem); % Root growth respiration In g C/plant/timestep (actually used)
        root_gro1 = root_growth*(suc_growth/totdem); %root growth actually achieved
        rlc_pt1 = rlc_pt*(suc_growth/totdem); % Leaf growth respiration In g C/plant/timestep (actually used)
        leaf_gro1 = leaf_growth*(suc_growth/totdem); %leaf growth actually achieved (g C/plant/timestep)
        suc_sta = 0;
end


