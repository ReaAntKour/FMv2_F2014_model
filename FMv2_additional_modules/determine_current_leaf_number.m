function [Appear,Leaf_mass,Si,Leaf_no,GC,GCstart] = determine_current_leaf_number(t,probability,CumThrm,phyllochron,Appear,Leaf_mass,Si,Leaf_no,GC,GCstart)
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

if  (CumThrm(t) - CumThrm(GCstart(end))) >= phyllochron % A growth cycle has been reached

    GC = GC+1;
    GCstart(GC+1) = t; % the timepoint when a new growth cycle starts

    Leaf_no(t) = Leaf_no(t-1) + binornd(1,probability); % A new leaf may/may not be initiated

    if Leaf_no(t) ~= Leaf_no(t-1) % A new leaf is initiated
       Appear(Leaf_no(t)) = t;
       Leaf_mass(t-1,Leaf_no(t))=0;
       Si(t-1,Leaf_no(t))=0;
    end              

else
    Leaf_no(t) = Leaf_no(t-1);        
end          
