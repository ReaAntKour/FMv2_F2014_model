function Func = P2011_dynamics(t,y,P,sunrise,sunset)

    P = num2cell(P);
    [q1,q2,q3,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,...
        p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,m1,m2,m3,m4,m5,m6,...
        m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,...
        m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39,n1,n2,n3,...
        n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,...
        g11,g12,g13,g14,g15,g16,a,b,c,d,e,f,yhb] = deal(P{:});

    Func = zeros(30, 1);

    cond.period=24;
    cond.dawn=sunrise;
    cond.photoperiod = sunset-sunrise;
    L = light_conditions(t,cond);
    YHBL=(yhb+L)/(yhb+1);

    %ODEs
	LHYm = y(1); % LHY mRNA                   
	P= y(2); % P
	GIZ = y(3); % GI-ZTL 
	GIE = y(4); % GI-ELF3 cytoplasm
	LHYp = y(5); % LHY prot
	TOCm = y(6); % TOC1 mRNA
	PRR9p = y(7); % PRR9 prot
	PRR5m = y(8); % PRR5 (NI) mRNA
	PRR5p = y(9); % PRR5 (NI) prot
	GIp = y(10); % GI prot cytoplasm
	TOCp = y(11); % TOC1 prot
	ZTL = y(12); % ZTL
	EC = y(13); % EC
	GIm = y(14); % GI mRNA
	PRR9m = y(15); % PRR9 mRNA
	PRR7m = y(16); % PRR7 mRNA
	PRR7p = y(17); % PRR7 prot
	ELF4m = y(18); % ELF4 mRNA
	ELF4p = y(19); % ELF4 prot
	LHYpm = y(20); % LHY prot modif.
	HY5 = y(21); % HY5 prot used only for fitting of COP1 parameters
	HFR1 = y(22); % HFR1 prot used only for fitting of COP1 parameters
	ELF3m = y(23); % ELF3 mRNA
	ELF3pc = y(24); % ELF3 cytoplasm
	ELF3pn = y(25); % ELF3 nuclear
	COP1pnn = y(26); % COP1 nuclear night
	COP1pnd = y(27); % COP1 nuclear day
	LUXm = y(28); % LUX mRNA
	LUXp = y(29); % LUX prot
	COP1pc = y(30); % COP1 cytoplasm

    Gn=p28*GIp/(p29+m19+p17*ELF3pn);
    EGn=(p18*GIE+p17*ELF3pn*Gn)/(m9*COP1pnn+m10*COP1pnd+p31);
    e34=p25*ELF4p*ELF3pn/(p26*LUXp+p21+m36*COP1pnn+m37*COP1pnd);

    Func(1) = 1* (q1*L*P+n1*g1^a/(g1^a+(PRR9p+PRR7p+PRR5p+TOCp)^a))-LHYm*(m1*L+m2*(1-L)); % LHY mRNA 
    Func(2) = p7*(1-L)*(1-P)-m11*P*L; % P
    Func(3) = p12*L*ZTL*GIp-p13*GIZ*(1-L)-m21*GIZ; % GI-ZTL 
    Func(4) = p17*ELF3pc*GIp-m9*GIE*COP1pc-p18*GIE+p31*EGn; % GI-ELF3 cytoplasm
    Func(5) = (p2+p1*L)*LHYm-m3*LHYp-p3*LHYp^c/(LHYp^c+g3^c); % LHY prot
    Func(6) = 1*n2*g4/(g4+EC)*(g5^e/(g5^e+LHYp^e))-TOCm*m5; % TOC1 mRNA
    Func(7) = p8*PRR9m-(m13+m22*(1-L))*PRR9p; % PRR9 prot
    Func(8) = 1*(n10*LHYpm^e/(g12^e+LHYpm^e)+n11*PRR7p^b/(g13^b+PRR7p^b))-m16*PRR5m; % PRR5 (NI) mRNA
    Func(9) = p10*PRR5m-(m17+m24*(1-L))*PRR5p; % PRR5 (NI) prot
    Func(10)= p11*GIm-m19*GIp-p12*L*ZTL*GIp+p13*GIZ*(1-L)-p17*ELF3pc*GIp-p28*GIp+p29*Gn; % GI prot cytoplasm
    Func(11)= p4*TOCm-m8*TOCp-(m6+m7*(1-L))*TOCp*(p5*ZTL+GIZ); % TOC1 prot
    Func(12)= 1*p14-m20*ZTL-p12*L*ZTL*GIp+p13*GIZ*(1-L); % ZTL
    Func(13)= p26*LUXp*e34-m36*EC*COP1pnn-m37*EC*COP1pnd-m32*EC*(1+p24*L*(EGn+Gn)^d/(g7^d+(EGn+Gn)^d)); % EC
    Func(14)= 1*(q2*L*P+g15^e/(g15^e+LHYp^e)*g14/(g14+EC)*n12)-GIm*m18; % GI mRNA
    Func(15)= 1*q3*L*P+g8/(g8+EC)*(n4+n7*LHYp^e/(LHYp^e+g9^e))-m12*PRR9m; % PRR9 mRNA
    Func(16)= 1*(n8*(LHYp+LHYpm)^e/(g10^e+(LHYp+LHYpm)^e)+n9*PRR9p^f/(g11^f+PRR9p^f))-m14*PRR7m; % PRR7 mRNA
    Func(17)= p9*PRR7m-PRR7p*(m15+m23*(1-L)); % PRR7 prot
    Func(18)= n13*(g6^e/(g6^e+LHYp^e))*g2/(g2+EC)-ELF4m*m34; % ELF4 mRNA
    Func(19)= p23*ELF4m-m35*ELF4p-p25*ELF3pn*ELF4p+p21*e34; % ELF4 prot
    Func(20)= p3*LHYp^c/(LHYp^c+g3^c)-m4*LHYpm; % LHY prot modif.
    Func(21)= p22-m38*HY5*COP1pnd-m25*HY5*COP1pnn; % HY5 prot used only for fitting of COP1 parameters
    Func(22)= p30-m28*HFR1*COP1pnn; % HFR1 prot used only for fitting of COP1 parameters
    Func(23)= 1*n3*g16^e/(g16^e+LHYp^e)-m26*ELF3m; % ELF3 mRNA
    Func(24)= p16*ELF3m-m9*ELF3pc*COP1pc-p17*ELF3pc*GIp-p19*ELF3pc+p20*ELF3pn; % ELF3 cytoplasm
    Func(25)= p19*ELF3pc-p20*ELF3pn-m29*ELF3pn*COP1pnn-m30*ELF3pn*COP1pnd-p25*ELF3pn*ELF4p+p21*e34-p17*ELF3pn*Gn; % ELF3 nuclear
    Func(26)= p6*COP1pc-n6*YHBL*P*COP1pnn-n14*COP1pnn-m27*COP1pnn*(1+p15*YHBL); % COP1 nuclear night
    Func(27)= 1*(n14*COP1pnn+n6*YHBL*P*COP1pnn)-m31*(1+m33*(1-YHBL))*COP1pnd; % COP1 nuclear day
    Func(28)= n13*(g6^e/(g6^e+LHYp^e))*g2/(g2+EC)-LUXm*m34; % LUX mRNA
    Func(29)= p27*LUXm-m39*LUXp-p26*LUXp*e34; % LUX prot
    Func(30)= 1*n5-p6*COP1pc-m27*COP1pc*(1+p15*YHBL); % COP1 cytoplasm
end