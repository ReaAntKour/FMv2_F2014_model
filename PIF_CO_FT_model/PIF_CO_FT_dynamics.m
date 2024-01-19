function dydt = PIF_CO_FT_dynamics(t,y,P,u,c)
% y(1) PIF4m
% y(2) PIF5m
% y(3) PIF
% y(4) phyB
% y(5) PR
% y(6) INT
% y(7) IAA29m
% y(8) ATHB2m
% y(9) CDF1m
% y(10) FKF1m
% y(11) FKF1
% y(12) CDF1
% y(13) COm
% y(14) CO
% y(15) FTm
% y(16) InducedC1
% y(17) InducedC2
% y(18) RepressedC1
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

%interpret light conditions at time t
L = light_conditions(t,c);

LHY = interp1(u.T,u.LHY,mod(t,c.period));
COP1n_n = interp1(u.T,u.COP1n_n,mod(t,c.period));
EC = interp1(u.T,u.EC,mod(t,c.period));
EC_adv = interp1(u.T,u.EC,mod(t+P.advance,c.period));
PRR9 = interp1(u.T,u.PRR9,mod(t,c.period));
PRR7 = interp1(u.T,u.PRR7,mod(t,c.period));
PRR5 = interp1(u.T,u.PRR5,mod(t,c.period));
TOC1 = interp1(u.T,u.TOC1,mod(t,c.period));
cP = interp1(u.T,u.cP,mod(t,c.period));
GIn = interp1(u.T,u.GIn,mod(t,c.period));


% Unpack and name model variables
y = num2cell(y);
[PIF4m,PIF5m,PIF,phyB,PR,INT,IAA29m,ATHB2m,...
    CDF1m,FKF1m,FKF1,CDF1,COm,CO,FTm,InducedC1,InducedC2,RepressedC1] = deal(y{:});

PIFtot = PIF;
PIFact_1 = PIF^2/(PIF^2+(P.g14*INT)^2); %IAA29
PIFact_2 = PIF^2/(PIF^2+(P.g15*INT)^2); %ATHB2
PIFact_3 = PIF^2/(PIF^2+(P.g16*INT)^2); %InducedC1
PIFact_4 = PIF^2/(PIF^2+(P.g17*INT)^2); %InducedC2
PIFact_5 = PIF^2/(PIF^2+(P.g18*INT)^2); %RepressedC1

% PIF model
dPIF4mdt = P.n7+P.n6*P.g7^P.e/(P.g7^P.e+EC_adv^P.e)- P.m8*PIF4m;
dPIF5mdt = P.n9+P.n8*P.g8^P.f/(P.g8^P.f+EC_adv^P.f)- P.m9*PIF5m;
dPIFdt = (P.p13*PIF5m+P.p12*PIF4m) - (P.m13 + P.m14*phyB)*PIF;
dphyBdt = P.p8*L*(1-phyB) - P.m10*phyB;
dPRdt = P.p9*L*(1-PR) - P.m11*PR;
dINTdt = P.p10+ P.p11*PR - P.m12*INT;
dIAA29mdt = P.n10+P.n11*PIFact_1^2/(P.g9^2+PIFact_1^2) - P.m15*IAA29m;
dATHB2mdt = P.n12+P.n13*PIFact_2^2/(P.g10^2+PIFact_2^2) - P.m16*ATHB2m;

% CO model
dCDF1mdt= (P.n1 + P.n2*(LHY^P.a)/((P.g1^P.a)+(LHY^P.a)))*(P.g2^P.b)/((P.g2^P.b) + (PRR9+PRR7+PRR5+TOC1)^P.b) - P.m1*CDF1m;
dFKF1mdt= P.q1*L*cP + P.n3*((P.g3^P.c)/((P.g3^P.c)+(LHY^P.c)))*(P.g4/(P.g4+EC)) - P.m3*FKF1m;
dFKF1dt= P.p4*FKF1m - P.p5*(P.m4 - L*(GIn /(P.k1 + GIn)))*FKF1;
dCDF1dt= P.p1*CDF1m - P.m2*(P.p2*FKF1*GIn + P.p3*GIn + 1)*CDF1;
dCOmdt= P.Bco + (P.g5^P.d)/((P.g5^P.d) + (CDF1^P.d))*(P.n4 + P.n5*(1-L)*(COP1n_n/(P.g6 + COP1n_n))) - P.m5*COm;
dCOdt= P.p6*COm - P.p7*(P.m6 + P.m7*(1-L)*COP1n_n - L*FKF1/(FKF1 + P.k2))*CO;

% FT model
dFTmdt = (P.n14 + P.n15*PIFtot/(P.g11 + PIFtot))*...
    (P.n16 + P.n17*P.g12/(P.g12 + CDF1))*(CO^P.h)/((CO^P.h)...
    + (P.g13^P.h)) - P.m17*FTm;

% Additional PIF targets
dInducedC1dt = P.n18 + P.n19*PIFact_3^2/(P.g19^2+PIFact_3^2) - P.m18*InducedC1;
dInducedC2dt = P.n20 + P.n21*PIFact_4^2/(P.g20^2+PIFact_4^2) - P.m19*InducedC2;
dRepressedC1dt = P.n22 + P.n23*P.g21^2/(P.g21^2+PIFact_5^2) - P.m20*RepressedC1;

% output
dydt_all = [dPIF4mdt
    dPIF5mdt
    dPIFdt
    dphyBdt
    dPRdt
    dINTdt
    dIAA29mdt
    dATHB2mdt
    dCDF1mdt
    dFKF1mdt
    dFKF1dt
    dCDF1dt
    dCOmdt
    dCOdt
    dFTmdt
    dInducedC1dt
    dInducedC2dt
    dRepressedC1dt];

% only simulate intended species
dydt = zeros(length(dydt_all),1);
dydt(P.sim_Sidx) = dydt_all(P.sim_Sidx);

end
