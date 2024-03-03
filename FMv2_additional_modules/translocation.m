function [leaf_trans,root_trans] = translocation(rosette_area,suc_c_interm,suc_equi,leaf_c,root_c)
%% transolcation of sugar between compartments if sucrose falls below equilibrium values
%
% Input:
%   rosette_area - Total leaf area of the rosette
%   suc_c_interm - Intermediate sucrose value
%   suc_equi - equilibrium sucrose plus hexose concentration in leaves (g C/m2 leaf)
%   leaf_c - Leaf carbon (gC/plant)
%   root_c - Root carbon (gC/plant)
%
% Output:
%   leaf_trans - Translocation from leaves (gC/plant)
%   root_trans - Translocation from roots (gC/plant)
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


suc_equi_plant = suc_equi*rosette_area; %Equilibrium sucrose plus hexose mass for whole plant (g C/plant)
root_and_leaf_c = root_c + leaf_c; %Total root and leaf C (g C/plant)

if      suc_c_interm <= suc_equi_plant %translocation needed
        
        scalingl = leaf_c/root_and_leaf_c;
        leaf_trans = (suc_equi_plant - suc_c_interm)*scalingl; %Translocation from leaves (gC/plant)
        scalingr = root_c/root_and_leaf_c;
        root_trans = (suc_equi_plant - suc_c_interm)*scalingr; %Translocation from root (g C/plant)
        
else
        leaf_trans = 0;
        root_trans = 0;
end        
    
