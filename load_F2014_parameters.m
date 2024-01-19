function param=load_F2014_parameters(genotype)
% F2014 clock model parameters

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

% Mutant variations:
% possible_mutants = {'cca1 lhy','cca1','CCA1OX','elf3','wt'}; % 'YHB' still to be added!

%% Define Parameters

% varying parameters in fitting (set 1 of 8)
a3 = 5.314;
a4 = 7.175;
a5 = 1.96;
a6 = 8.874;
a7 = 1.528;
a8 = 3.014;

r1 = 5.673;
r2 = 1.395;
r3 = 8.556;
r4 = 3.422;
r5 = 0.1532;
r6 = 25.77;
r7 = 12.27;
r8 = 1.753;
r9 = 30.75;
r10 = 9.981;
r11 = 1.941;
r12 = 4.987;
r13 = 44.99;
r14 = 3.904;
r15 = 5.344;
r16 = 9.74;
r17 = 2.25;
r18 = 41.94;
r19 = 15.93;
r20 = 0.109;
r21 = 3.336;
r22 = 56.5;
r23 = 3.28;
r24 = 4.315;
r25 = 41.19;
r26 = 4.774;
r27 = 1.34;
r28 = 5.91;
r29 = 0.2745;
r30 = 1.411;
r31 = 0.03975;
r32 = 6.967;
r33 = 0.7146;
r34 = 2.874;
r35 = 0.06041;
r36 = 0.08923;
r37 = 9.958;
r38 = 3.887;
r40 = 1.422;
% r41 = 1.856;

n1 = 1;

f1 = 0.06021;
f2 = 0.02732;
f3 = 0.2654;
f4 = 0.2687;
f5 = 0.2899;
f6 = 0.08977;

t5 = 0.2518;
t6 = 0.112;
t7 = 3.772;
t8 = 28.34;
t9 = 2.543;

m1 = 0.6127;
m3 = 0.6154;
m4 = 0.4322;
m5 = 1.869;
m6 = 0.451;
m7 = 1.277;
m8 = 0.1;
m9 = 0.1886;
m10 = 0.01001;
m11 = 0.7611;
m12 = 2.57;
m13 = 0.6654;
m14 = 0.5015;
m15 = 0.1791;
m16 = 0.5383;
m17 = 0.07499;
m18 = 4.505;
m22 = 0.3007;
m23 = 0.08544;
m24 = 1.5;
m25 = 0.6526;
m26 = 1.032;
m28 = 0.0382;
m29 = 0.01001;
m30 = 4.949;
m32 = 9.998;
m34 = 0.1561;
m35 = 0.9413;
m36 = 0.506;
m37 = 0.01001;
m38 = 1.774;
m39 = 0.2001;
m42 = 0.9229;
m43 = 1.139;
m44 = 0.6769;
m45 = 0.8039;
m46 = 5.273;
m47 = 0.2466;

p11 = 1.912;
p16 = 0.1211;
p23 = 1.011;
p25 = 1.003;
p28 = 1.061;
p29 = 10.18;

q1 = 0.2607;
q3 = 0.4659;

% non-varying parameters
m19 = 0.2;
m20 = 1.8;
m21 = 0.1;
m27 = 0.1;
m31 = 0.3;
m33 = 13;
n5 = 0.23;
n6 = 20;
n14 = 0.1;
p6 = 0.6;
p7 = 0.3;
p10 = 0.2;
p12 = 8;
p13 = 0.7;
p14 = 0.3;
p15 = 3;

%% Mutants 

% modify according to genotype:
for i = 1:length(genotype)
    if strcmp(genotype{i},'cca1 lhy')
        q1 = 0;
        n1 = 0;
	elseif strcmp(genotype{i},'cca1')
        f5 = 0;
    elseif strcmp(genotype{i},'CCA1OX')
        % CCA1-OX genotype{i}:
        q1 = 0;
        n1 = 5;
    elseif strcmp(genotype{i},'elf3')
        p16 = 0; % no ELF3 protein production. (mRNA is not directly affected)
    elseif strcmp(genotype{i},'YHB')
        error('still to be added!')
	elseif strcmp(genotype{i},'wt')
		% Do nothing
	elseif strcmp(genotype{i},'wildtype')
		% Do nothing
	elseif strcmp(genotype{i},'wild-type')
		% Do nothing
	else
		error('No such mutant')
    end
end

%% Define param
param = [a3,a4,a5,a6,a7,a8,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r30,r31,r32,r33,r34,r35,r36,r37,r38,r40,n1,f1,f2,f3,f4,f5,f6,t5,t6,t7,t8,t9,m1,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m22,m23,m24,m25,m26,m28,m29,m30,m32,m34,m35,m36,m37,m38,m39,m42,m43,m44,m45,m46,m47,p11,p16,p23,p25,p28,p29,q1,q3,m19,m20,m21,m27,m31,m33,n5,n6,n14,p6,p7,p10,p12,p13,p14,p15];
