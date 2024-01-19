function u = wrap_P2011_model_dynamics(T,Y,P) % only called from "simulate_PIF_CO_FT_model"
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

u=struct;
u.T = T;
nT = length(T);
varnames = {'LHY','EC','PRR9','PRR7','PRR5','TOC1','cP','COP1n_n','GIn'};
nV = length(varnames);
for i = 1:nV
    u.(varnames{i}) = zeros(nT,1);
end

P = num2cell(P);
[q1,q2,q3,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,...
    p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,m1,m2,m3,m4,m5,m6,...
    m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,...
    m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39,n1,n2,n3,...
    n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,...
    g11,g12,g13,g14,g15,g16,a,b,c,d,e,f] = deal(P{:});
    
    u.LHY = Y(:,5);
    u.COP1n_n = Y(:,26);
    u.EC = Y(:,13);
    u.PRR9 = Y(:,7);
    u.PRR7 = Y(:,17);
    u.PRR5 = Y(:,9);
    u.TOC1 = Y(:,11);
    u.cP = Y(:,2);
    u.GIn = p28.*Y(:,10)./(p29+m19+p17.*Y(:,25));
