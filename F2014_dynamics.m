function out=F2014_dynamics(t,y,param,sunrise,sunset)
% F2014 clock model dynamics

% Model
% Fogelmark K, Troein C (2014) Rethinking Transcriptional Activation in the Arabidopsis Circadian Clock. PLoS Comput Biol 10(7): e1003705. doi:10.1371/journal.pcbi.1003705

% Code:
%   Copyright 2023 Rea L Antoniou-Kourounioti, and The University of Glasgow
%%%%%%%%%%% Decide on Licence and update copyright and licence accordingly
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

out=zeros(35,1);

% Light conditions
cond.period=24;
cond.dawn=sunrise;
cond.photoperiod = sunset-sunrise;
L = light_conditions(t,cond);

% Define the variables
LHYm = y(1);
LHYp = y(2);
CCA1m = y(3);
CCA1p = y(4);
P = y(5);
PRR9m = y(6);
PRR9p = y(7);
PRR7m = y(8);
PRR7p = y(9);
PRR5m = y(10);
PRR5c = y(11);
PRR5n = y(12);
TOC1m = y(13);
TOC1n = y(14);
TOC1c = y(15);
ELF4m = y(16);
ELF4p = y(17);
ELF4d = y(18);
ELF3m = y(19);
ELF3p = y(20);
ELF34 = y(21);
LUXm = y(22);
LUXp = y(23);
COP1c = y(24);
COP1n = y(25);
COP1d = y(26);
ZTL = y(27);
ZG = y(28);
GIm = y(29);
GIc = y(30);
GIn = y(31);
NOXm = y(32);
NOXp = y(33);
RVE8m = y(34);
RVE8p = y(35);


% Define Parameters
param = num2cell(param);
[a3,a4,a5,a6,a7,a8,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,...
	r16,r17,r18,r19,r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r30,r31,...
	r32,r33,r34,r35,r36,r37,r38,r40,n1,f1,f2,f3,f4,f5,f6,t5,t6,t7,t8,...
	t9,m1,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m22,...
	m23,m24,m25,m26,m28,m29,m30,m32,m34,m35,m36,m37,m38,m39,m42,m43,m44,...
	m45,m46,m47,p11,p16,p23,p25,p28,p29,q1,q3,m19,m20,m21,m27,m31,m33,...
	n5,n6,n14,p6,p7,p10,p12,p13,p14,p15,yhb] = deal(param{:});
BL=L;
L=max(L+yhb,1);
D=1-L;

% Equations
LC = (LHYp + f5 * CCA1p);
LCcommon = (q1 * L * P + n1) / (1 + (r1 * PRR9p)^2 + (r2 * PRR7p)^2 + (r3 * PRR5n)^2 + (r4 * TOC1n)^2);
EC = ((LUXp + f6 * NOXp) * (ELF34 + f1 * ELF3p)) / (1 + f3 * (LUXp + f2 * NOXp) + f4 * (ELF34 + f1 * ELF3p));
P5trans = t5 * PRR5c - t6 * PRR5n;
Ttrans = t7 * TOC1c - (t8) / (1 + m37 * PRR5n) * TOC1n;
E34prod = p25 * ELF3p * ELF4d;
E3deg = m30 * COP1d + m29 * COP1n + m9 + m10 * GIn;
ZGprod = p12 * ZTL * GIc - (p13 * (1-BL) + p10 * (BL)) * ZG;
ELF3tot = ELF3p + ELF34;
Gtrans = p28 * GIc - (p29) / (1 + t9 * ELF3tot) * GIn;


% LHYm
out(1) = (LCcommon) / (1 + (r11 * LC)^2) - m1 * LHYm;

% LHYp
out(2) = (L + m4 * D) * LHYm - m3 * LHYp;

% CCA1m
out(3) = LCcommon - m1 * CCA1m;

% CCA1p
out(4) = (L + m4 * D) * CCA1m - m3 * CCA1p;

% P
out(5) = p7 * D * (1 - P) - m11 * P * L;

% PRR9m
out(6) = q3 * P * L - m12 * PRR9m + (1 + a3 * r33 * RVE8p) / ((1 + r33 * RVE8p) * (1 + (r5 * LC)^2) * (1 + (r6 * EC)^2) * (1 + (r7 * TOC1n)^2) * (1 + (r40 * PRR5n)^2));

% PRR9p
out(7) = PRR9m - m13 * PRR9p;

% PRR7m
out(8) = 1/((1 + (r8 * LC)^2) * (1 + (r9 * EC)^2) * (1 + (r10 * TOC1n)^2) * (1 + (r40 * PRR5n)^2)) - m14 * PRR7m;

% PRR7p
out(9) = PRR7m - (m15 + m23 * D) * PRR7p;

% PRR5m
out(10) = (1 + a4 * r34 * RVE8p) / ((1 + r34 * RVE8p) * (1 + (r12 * LC)^2) * (1 + (r13 * EC)^2) * (1 + (r14 * TOC1n)^2)) - m16 * PRR5m;

% PRR5c
out(11) = PRR5m - (m17 + m24 * ZTL) * PRR5c - P5trans;

% PRR5n
out(12) = P5trans - m42 * PRR5n;

% TOC1m
out(13) = (1 + a5 * r35 * RVE8p) / ((1 + r35 * RVE8p) * (1 + (r15 * LC)^2) * (1 + (r16 * EC)^2) * (1 + (r17 * TOC1n)^2)) - m5 * TOC1m;

% TOC1n
out(14) = Ttrans - (m43) / (1 + m38 * PRR5n) * TOC1n;

% TOC1c
out(15) = TOC1m - (m8 + m6 * ZTL) * TOC1c - Ttrans;

% ELF4m
out(16) = (1 + a6 * r36 * RVE8p) / ((1 + r36 * RVE8p) * (1 + (r18 * EC)^2) * (1 + (r19 * LC)^2) * (1 + (r20 * TOC1n)^2)) - m7 * ELF4m;

% ELF4p
out(17) = p23 * ELF4m - m35 * ELF4p - ELF4p^2;

% ELF4d
out(18) = ELF4p^2 - m36 * ELF4d - E34prod;

% ELF3m
out(19) = 1/(1 + (r21 * LC)^2) - m26 * ELF3m;

% ELF3p
out(20) = p16 * ELF3m - E34prod - E3deg * ELF3p;

% ELF34
out(21) = E34prod - m22 * ELF34 * E3deg;

% LUXm
out(22) = (1 + a7 * r37 * RVE8p) / ((1 + r37 * RVE8p) * (1 + (r22 * EC)^2) * (1 + (r23 * LC)^2) * (1 + (r24 * TOC1n)^2)) - m34 * LUXm;

% LUXp
out(23) = LUXm - m39 * LUXp;

% COP1c
out(24) = n5 - p6 * COP1c - m27 * COP1c * (1 + p15 * L);

% COP1n
out(25) = p6 * COP1c - (n14 + n6 * L * P) * COP1n - m27 * COP1n * (1 + p15 * L);

% COP1d
out(26) = (n14 + n6 * L * P) * COP1n - m31 * (1 + m33 * D) * COP1d;

% ZTL
out(27) = p14 - ZGprod - m20 * ZTL;

% ZG
out(28) = ZGprod - m21 * ZG;

% GIm
out(29) = (1 + a8 * r38 * RVE8p) / ((1 + r38 * RVE8p) * (1 + (r25 * EC)^2) * (1 + (r26 * LC)^2) * (1 + (r27 * TOC1n)^2)) - m18 * GIm;

% GIc
out(30) = p11 * GIm - ZGprod - Gtrans - m19 * GIc;

% GIn
out(31) = Gtrans - m19 * GIn - m25 * ELF3tot * (1 + m28 * COP1d + m32 * COP1n) * GIn;

% NOXm
out(32) = 1/((1 + (r28 * LC)^2) * (1 + (r29 * PRR7p)^2)) - m44 * NOXm;

% NOXp
out(33) = NOXm - m45 * NOXp;

% RVE8m
out(34) = 1/(1 + (r30 * PRR9p)^2 + (r31 * PRR7p)^2 + (r32 * PRR5n)^2) - m46 * RVE8m;

% RVE8p
out(35) = RVE8m - m47 * RVE8p;
