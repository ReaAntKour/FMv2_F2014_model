function [output,sim_data] = simulate_FM(hour,T,sunrise,sunset,CO2,PAR,Photoperiod,clock_parameters,starch_parameters,p,d,mf_use,run_phenology_model,fileID)
%% simulate the Framework Model v2. This function is the core function,
%
% Input:
%   hour - vector of timesteps
%   T - vector of temperatures over timesteps given by hour
%   sunrise - daily time of sunrise (integer, in hours)
%   sunset - daily time of sunset (integer, in hours)
%   CO2 - CO2 concentration
%   PAR - incident photosynthetically active radiation
%   Photoperiod - Obsolete parameter, kept in input for backwards compatibility
%   clock_parameters - vector of P2011 circadian clock model parameters
%   starch_parameters - structure of starch model parameters
%   p - vector of Framework Model parameters (beyond those used for the clock and starch models)
%   d - dry weight proportion
%   mf_use - fraction of malate and fumarate turned over during the night
%   run_phenology_model - boolean variable determining whether to run the phenology model
%   fileID - file to write output data to
%
% Output:
%   output - subset of output data - fresh weight, assimilation, respiration, and starch content on day 38
%   sim_data - structure containing nested files with complete simulation data
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

%specify Col accession for phenology model threshold
flowering_thresh_geno = 2;

%Calling for parameters
%______________________

addpath('PIF_CO_FT_model')

probability = 1; %deterministic


%Parameters for growth model
%__________________________________________

Tb = p(1); %Base temperature for thermal-time calculation
TT0 = p(19); %degree days required from seed germination until plant emergence (stage 1.0 in Boyes et al. 2001)


TS = p(26); %degree days of leaf lifespan (from the point of leaf initiation)
TLcot = p(27); %cotyledon expansion period
TLtrue = p(28); %expansion period of true leaves
Root_LC = p(29); %root expansion period 
leaf_factor = p(30); %carbon-dry mass conversion
root_factor = p(31); 
seed_input = p(32);
cot_area = p(33); %when they are fully opened
sucrose_initial = p(34); %initial sucrose content per unit area
starch_initial = p(35); %initial starch content per initial sucrose content
petiole = p(36); %maximum petiole elongation factor 
min_alpha = p(65); %minimum zenithal angle
inc_alpha = p(66); %maximum increase in zenithal angle from the minimum value
juv_TT = p(67); %thermal time since plant emergence for the juvenile-to-adult transition
juv_rate = p(68); %rate of leaf appearance at juvenile stage
juv_interc = p(81);%y-intercept for the juvenile rate of leaf appearance
ad_rate = p(69); %rate of leaf appearance at adult stage
ad_interc = p(70); %y-intercept for the adult rate of leaf appearance
SLA_cot = p(71); %SLA for cotyledons (in m2/g dry mass) Fig. 3 Christophe et al (2008)
SLA_exp = p(72); %curve shape for SLA against thermal time
juv_phyllochron = p(82); %30.3
ad_phyllochron = p(83); %11.9
photosyn_efficiency = p(84); %0.88
ss_efficiency = p(60); %0.6

%Change SLA depending on photoperiod, such that SLA_cot takes a value of
%0.88 in 18h photoperiods, or has the default value for photoperiods of 12h
%or less.
SLA_cot = SLA_cot*(1-(max(Photoperiod,12)-12)/6*0.12);


%Expansion period for the root system
%____________________________________

PR=p(23); %root expansion
aR=p(24);
bR=p(25);

TR = 1000;%Root_LC*CumThrm(end);
MaxfRpoint = TR/(1+(bR-1)/(aR-1))-0.5; %determining the maximum point by solving the function derivative
MaxfR = ((MaxfRpoint+0.5)/TR)^(aR-1)*(1-(MaxfRpoint+0.5)/TR)^(bR-1);  

%Maximum values of leaf sink variation
%_____________________________________

PL=p(20); %leaf expansion
aL=p(21);
bL=p(22);

Maxfcotpoint = TLcot/(1+(bL-1)/(aL-1))-0.5; %determining the maximum point by solving the function derivative
Maxfcot = ((Maxfcotpoint+0.5)/TLcot)^(aL-1)*(1-(Maxfcotpoint+0.5)/TLcot)^(bL-1);

MaxfLpoint = TLtrue/(1+(bL-1)/(aL-1))-0.5; %determining the maximum point by solving the function derivative
MaxfL = ((MaxfLpoint+0.5)/TLtrue)^(aL-1)*(1-(MaxfLpoint+0.5)/TLtrue)^(bL-1);  
          
%Initialise matrices that will fill up
fL = [];
fR = [];
rsratio = [];
efficiency = [];
growth_capacity = [];



%%%This should be a separate function
%initialise clock to starting conditions
clock_state_0=[1.0151 0.956 0.0755 0.0041 0.506 0.0977 0.0238 0.0731 0.0697 0.0196 0.0435 0.2505 0.0709 0.1017 0.0658 0.4016 0.1167 0.1012 0.207 0.0788 0.3102 0.0553 0.2991 0.1503 0.0286 0.65 0.2566 0.1012 0.576 0.3269]; %12:12 wt;
for i = 1:5
    clock_output=circadian_module(sunrise(1),sunset(1),clock_state_0,clock_parameters);
    %work out the clock state at ZT24 i.e. at the end of the previous day
    clock_state_0 = interp1q(clock_output.T,clock_output.Y,24);
end

starch_module_state = [1,1,1,1];
FT_module_state =ones(1,18);

% set whether the phenology model controls the end of the simulation
N_max_days = 40; %default number of days max simulation
if run_phenology_model
    %upper ceiling on the number of days max simulation, if the phenology model is being run
    N_max_days = 90;
end
input_size = N_max_days*24;


day_idx = 1;
has_flowered = false;
has_emerged = false;
CumPhenThrm=0;

while day_idx <= N_max_days && ~(has_flowered)
    %initial timepoint
    t = (day_idx-1)*24+1;
    
    %run clock model for this day
    clock_output=circadian_module(sunrise(t),sunset(t),clock_state_0,clock_parameters);
    %work out the clock state at ZT24 i.e. at the end of the previous day
    clock_state_0 = interp1q(clock_output.T,clock_output.Y,24);
    
    %run phenology model
    if run_phenology_model
    [DayPhenThrm,FT_module_state] = phen(T,t,sunrise(t),sunset(t),flowering_thresh_geno,clock_parameters,clock_output,FT_module_state,p);
    else
        DayPhenThrm = 0;
    end
    CumPhenThrm = DayPhenThrm+CumPhenThrm;
    has_flowered = flowering_threshold_test(CumPhenThrm,flowering_thresh_geno,p);
    
    for hour_idx = 1:24
        %timepoint
        t = (day_idx-1)*24+hour_idx;
        
        Thrm(t)=max(0,(T(t)-Tb))/24; %in degree days   
        CumThrm(t) = sum(Thrm(1:t)); %cumulative thermal time
        
        % These may be overwritten - if not then they are zero
        Total_leaf_mass(t) = 0;
        Root_mass(t)       = 0;
        Starch_carbon(t)   = 0;
        MF_carbon(t)       = 0;
        Leaf_no(t)         = 0; 

        if ~(has_emerged)
            %Plant has not yet emerged
            %______________________________________________

            if CumThrm(t) >= TT0
                eme = t;
                has_emerged = 1;
                %For the first growth cycle (from germination to TT0):
                %____________________________________________________
                %%%FUNC: initialise seedling state at emergence

                Seed_mass = seed_input; %Seed dry mass (g), emptied after the first cycle to reach stage 1.0
                Q(eme) = Seed_mass;
                fR = ((CumThrm(eme)+0.5)/TR)^(aR-1)*(1-(CumThrm(eme)+0.5)/TR)^(bR-1)/MaxfR;  
                D(eme) = PR*fR + 2*PL;
                Root_mass(eme) = PR*fR*Q(eme)/D(eme); %root dry mass (g)

                Si(eme,1:2) = cot_area; %area for EACH cotyledon (m2), cycle 1 for leaf ranked 1 and 2
                Leaf_mass(eme,1:2) = Si(eme,1)/SLA_cot; %dry mass for EACH cotyledon (g)
                Hypocotyl_mass = Seed_mass - Root_mass(eme) - sum(Leaf_mass(eme,:));
                Total_shoot(eme) = Seed_mass - Root_mass(eme);
                Destructive_area(eme) = sum(Si(eme,:));
                Leaf_no(eme) = 2;
                Appear(1:2) = eme; %The timepoint when the cotyledons start expanding using photosynthate


                %Initialisation
                %______________

                S(eme) = sum(Si(eme,1:2));
                S_intercept(eme) = S(eme)*cos(10/180*pi); 
                Func_area(eme)=S(eme);
                %%%This could give the impression that another timestep can
                %%%be chosen
                steps_perhour = 1;
                steps_perday = 24*steps_perhour;
                timestep = 1/steps_perday; % Time step length used to solve the model

                Leaf_carbon(eme) = leaf_factor*sum(Leaf_mass(eme,:));
                Root_carbon(eme) = root_factor*Root_mass(eme);    
                Sucrose_carbon(eme) = sucrose_initial*S(eme); %Initial amount of sucrose carbon (g/plant)
                Starch_carbon(eme) = starch_initial*Sucrose_carbon(eme); %Initial amount of starch carbon (g/plant)
                MF_carbon(eme) = 0.4*Starch_carbon(eme); %Initial amount of malate+fumarate carbon (g/plant)
                rsratio(eme) = Root_carbon(eme)/Leaf_carbon(eme);
                growth_capacity(eme) = p(63);

                if      hour(eme) < sunrise(eme)
                        sta_c_endday = Starch_carbon(eme)/(p(61)*((sunrise(eme)-hour(eme))/(sunrise(eme)+24-sunset(eme-1)))+1-p(61));
                elseif  hour(eme) > sunset(eme)
                        sta_c_endday = Starch_carbon(eme)/(p(61)*((24-hour(eme)+sunrise(eme+1))/(24-sunset(eme)+sunrise(eme+1)))+1-p(61));
                else    
                        sta_c_endday = 0;
                end        

                MF_c_endday = 0.4*sta_c_endday;

                [rlc_pt1(eme),rrc_pt1(eme),leaf_res(eme),root_res(eme),rgtot(eme),rmtot(eme),totalC(eme),assim(eme)] = ...
                                                      ini_carbon_balance_MF(T(eme),CO2(eme),PAR(eme),...
                                                      sunrise(eme),sunset(eme),rsratio(eme),...
                                                      S_intercept(eme),Leaf_carbon(eme),...
                                                      Root_carbon(eme),Sucrose_carbon(eme),...
                                                      Starch_carbon(eme),MF_carbon(eme),...
                                                      timestep,growth_capacity(eme),p,mf_use);
                %%%END FUNC: initialise seedling state at emergence

                GC = 1; %number of growth cycle elapsed                                  
                GCstart(GC+1) = eme; %The timepoint when the second growth cycle starts
            end
        else
            %Plant is growing
            %________________

            %To determine the rosette area at the previous time point for light interception
            %_________________________________________________________________________________


            %Determining for each leaf rank i (functioning/senesced) at the previous timepoint:

            for i = 1:Leaf_no(t-1)
                if CumThrm(t-1) < (CumThrm(Appear(i))+TS)
                   fS(i) = 1;
                else
                   fS(i) = 0;
                end
            end


            %To determine rosette structure (at the previous time point):

            [Max_area,i_max] = max(Si(t-1,:));
            i_max = max(i_max,2);

            %%%DS note - change loop or remove mult by indicator variable fS(i)
            for i=1:Leaf_no(t-1) %for each leaf rank i

                if      i<=i_max

                        alpha(i)=min_alpha; %zenithal angle
                else 
                        alpha(i)=min_alpha+inc_alpha*(i-i_max)/(Leaf_no(t-1)-i_max); 
                end
               
                S(t-1,i)=Si(t-1,i)*cos(alpha(i)/180*pi)*fS(i); %Projected functioning leaf area 
                Sf(t-1,i)=Si(t-1,i)*fS(i); %Functioning leaf area 

            end


            %Rosette area (at the previous time point):      
            [S_intercept,Func_area] = calculate_area(t,S,petiole,input_size,Appear,Leaf_no,Sf,S_intercept,Func_area);


            %To determine whole-plant carbon balance (at the current time point):
            %___________________________________________________________________


            if CumThrm(t) < (juv_TT+TT0)
               phyllochron = juv_phyllochron*probability;
            else
               phyllochron = ad_phyllochron*probability;
            end

            %To determine current leaf number:
            [Appear,Leaf_mass,Si,Leaf_no,GC,GCstart] = determine_current_leaf_number(...
                t,probability,CumThrm,phyllochron,Appear,Leaf_mass,Si,Leaf_no,GC,GCstart);

            %To calculate sink variation for leaf:
            [fL,fR,rsratio] = calculate_root_shoot_sink_variation(t,Leaf_no,TLcot,Maxfcot,...
                TLtrue,MaxfL,aL,bL,CumThrm,Appear,TR,MaxfR,aR,bR,PR,PL,rsratio,...
                leaf_factor,root_factor,fL,fR);
            %Calculating root-to-shoot ratio (in g Carbon):
            leafconversion = leaf_factor;     % Gorsuch et al 2010a and 2010b

            rootconversion = root_factor; %g Carbon per g dry mass
                             %from Kumar et al (Agroforest Syst 2010) around 30-35%
                             %from UK BioChar Research Centre (UoE) around 35 % in
                             %wheat, cited as Prendergast-Miller M and Sohi SP 2010. 
                             %Investigating biochar impacts on plant roots and root carbon.
                             %Poster presented in the Organic Matter Stabilization
                             %and Ecosystem Functions session at Soil Organic
                             %Matter Conference, Cote d'Azur, France (Sept 2010)

            % To determine light status:                
            %islight = determine_light_status(t,hour,sunrise,sunset)
            if  hour(t) > sunrise(t) && hour(t) <= sunset(t)            
                is_light(t) = 1;            
            else
                is_light(t) = 0;
            end   


            %To determine starch status:
            %sta_c_endday = determine_starch_status(t,sta_c_endday,Starch_carbon,is_light)
            if      is_light(t) == 1     
                    sta_c_endday = 0;

            elseif  is_light(t) == 0 && is_light(t-1) == 1  
                    sta_c_endday = Starch_carbon(t-1);
            else   
                    sta_c_endday = sta_c_endday;
            end


            %To determine malate and fumarate (MF) status:
            %MF_c_endday = determine_MF_status(t,MF_c_endday,MF_carbon,is_light)

            if      is_light(t) == 1     
                    MF_c_endday = 0;

            elseif  is_light(t) == 0 && is_light(t-1) == 1  
                    MF_c_endday = MF_carbon(t-1);
            else   
                    MF_c_endday = MF_c_endday;
            end


            %To determine the amount of sugar accumulated over a day
            [efficiency,growth_capacity]=calculate_sugar_acc(t,p,eme,...
                efficiency,growth_capacity,is_light);

            
            [rlc_pt1(t),rrc_pt1(t),leaf_res(t),root_res(t),...
             Leaf_carbon(t),Root_carbon(t),Sucrose_carbon(t),Starch_carbon(t),MF_carbon(t),rgtot(t),rmtot(t),totalC(t),assim(t),...
             suc_sta(t),net_rate(t),new_starch_module_state] = plant_carbon_balance_MF_clock(hour_idx,T(t),CO2(t),PAR(t),sunrise(t),sunset(t),is_light(t),...
                   rsratio(t),S_intercept(t-1),Leaf_carbon(t-1),Root_carbon(t-1),...
                   Sucrose_carbon(t-1),Starch_carbon(t-1),MF_carbon(t-1),rgtot(t-1),rmtot(t-1),timestep,sta_c_endday,MF_c_endday,efficiency(t),...
                   growth_capacity(t),clock_output,starch_module_state,starch_parameters,p,mf_use);

            starch_module_state = new_starch_module_state;

            %Root mass:
            Root_mass(t) = Root_carbon(t)/rootconversion; %root dry mass (g)
            %Actual root C growth
            CarbonR = Root_carbon(t) - Root_carbon(t-1);    

            %Total dry biomass accumulation for leaves (g): 
            CarbonL = Leaf_carbon(t) - Leaf_carbon(t-1);
            BiomassL = CarbonL/leafconversion; 


            %Leaf demand:
            for     i = 1:Leaf_no(t)
                    Ld(t,i) = PL*fL(t,i);
            end           

            Totalleafdemand = sum(Ld(t,:));


            %To calculate SLA:        
            SLA(t) = SLA_cot*exp(-SLA_exp*(CumThrm(t)-TT0)); %SLA (in m2/g dry mass) Fig. 3 Christophe et al (2008)


            %Leaf mass and area:
            for     i = 1:Leaf_no(t)
                    Leaf_mass(t,i) = Leaf_mass(t-1,i) + Ld(t,i)/Totalleafdemand*BiomassL;
                    Si(t,i) = max(Leaf_mass(t,i)*SLA(t),Si(t-1,i));
            end        

            Destructive_area(t) = sum(Si(t,1:Leaf_no(t)));        
            Total_leaf_mass(t) = sum(Leaf_mass(t,:));
            Total_shoot(t) = Total_leaf_mass(t) + Hypocotyl_mass; 


            Cstatus(t) = (CarbonL+CarbonR)/(Totalleafdemand+fR(t)*PR);
            GroRes(t) = rgtot(t)-rgtot(t-1);
            MainRes(t) = rmtot(t)-rmtot(t-1);
            Total_plant(t) = Total_shoot(t) + Root_mass(t);  

            Normalisedgrowth(t) = (Total_plant(t)-Total_plant(t-1))/Total_plant(t-1);
            Normalisedcost(t) = net_rate(t)/Total_plant(t-1);

        end

        % Output data as used by the web app
        if exist('fileID','var')
            fprintf(fileID,'t:%d,rosette:%f,rosette_fw:%f,root:%f,starch:%f,mf:%f,leaf_num:%f,temp:%f,co2:%f,light:%f\n', t, Total_leaf_mass(t), Total_leaf_mass(t)/d, Root_mass(t), Starch_carbon(t), MF_carbon(t),Leaf_no(t),T(t),CO2(t),PAR(t) );
        end
    end
    
    
    %Update day index
    day_idx = day_idx+1;
end

% Output that the simulation has finished
if exist('fileID','var')
    fprintf(fileID, 'Finished\n');
end

%-------------------------------End-of-model-------------------------------

sim_data = struct();
output = [];

try
    %Try to instantiate a complete output as far as possible
    sim_data.metadata.hour = hour;
    sim_data.metadata.sunrise = sunrise;
    sim_data.metadata.sunset = sunset;
    sim_data.metadata.CO2 = CO2;
    sim_data.metadata.PAR = PAR;
    sim_data.metadata.Photoperiod = sunset - sunrise;
    sim_data.metadata.clock_parameters = clock_parameters;
    sim_data.metadata.starch_parameters = starch_parameters;
    sim_data.metadata.p = p;
    sim_data.metadata.d = d;
    sim_data.metadata.mf_use = mf_use;
    sim_data.metadata.run_phenology_model = run_phenology_model;
    
    sim_data.input_size = input_size;

    sim_data.S_intercept = S_intercept;
    sim_data.totgrowth = sum(Normalisedgrowth(eme+25:end));
    sim_data.totcost = sum(Normalisedcost(eme+25:end));
    sim_data.minratio = sim_data.totcost/sim_data.totgrowth;
    sim_data.final_leaf_no = Leaf_no(end);
    sim_data.final_S_intercept = S_intercept(end);
    sim_data.final_root_mass = Root_mass(end);

    sim_data.Total_shoot = Total_shoot;
    sim_data.Starch_carbon = Starch_carbon;
    sim_data.MF_carbon = MF_carbon;

    sim_data.Fresh_weight = Total_shoot(end)/d; %g FW rosette
    sim_data.starch = Starch_carbon; %g C
    sim_data.sucrose = Sucrose_carbon; %g C
    sim_data.NPP = assim - rlc_pt1 - rrc_pt1 - leaf_res - root_res; %g C/h
    sim_data.Gas_exchange_perFW = sim_data.NPP./(Total_shoot/d); %g C/g FW/h
    sim_data.starch_perFW = sim_data.starch./(Total_shoot/d); %g C/g FW
    sim_data.sucrose_perFW = sim_data.sucrose./(Total_shoot/d); %g C/g FW


    sim_data.FW_ED = Total_shoot(input_size+sunset-24)/d; %g
    sim_data.FW_EN = Total_shoot(input_size+sunrise-24)/d; %g
    sim_data.AD = sim_data.Gas_exchange_perFW(input_size+sunrise+2-24); %g C/g FW/h
    sim_data.AR1 = sim_data.Gas_exchange_perFW(input_size-24+sunset+5); %g C/g FW/h
    sim_data.AR2 = sim_data.Gas_exchange_perFW(input_size-24+sunrise-1); %g C/g FW/h


    sim_data.FW = Total_shoot/d;
    sim_data.starchC6 = Starch_carbon/12/6*10^6./sim_data.FW;
    sim_data.malatefumarateC4 = MF_carbon/12/4*10^6./sim_data.FW;
    sim_data.sugar = Sucrose_carbon/12/12*10^6./sim_data.FW;
    sim_data.gFW_22 = sim_data.FW(22*24);
    sim_data.gFW_28 = sim_data.FW(28*24);
    sim_data.gFW_29 = sim_data.FW(29*24);
    sim_data.gFW_38 = sim_data.FW(38*24);
    sim_data.starchC6_ED_28 = sim_data.starchC6(28.5*24);
    sim_data.starchC6_EN_28 = sim_data.starchC6(29*24);
    sim_data.AperFW = sim_data.Gas_exchange_perFW(38*24+1)/12*10^6; %micromol CO2/g FW/h
    sim_data.RperFW = sim_data.Gas_exchange_perFW(38*24-1)/12*10^6; %micromol CO2/g FW/h
    sim_data.Aperarea = sim_data.NPP(38*24+1)/S_intercept(38*24)/12*10^6/10000; %micromol CO2/cm2/h
    sim_data.Rperarea = sim_data.NPP(38*24-1)/S_intercept(38*24-2)/12*10^6/10000; %micromol CO2/cm2/h
    output = [sim_data.gFW_38,sim_data.AperFW,sim_data.RperFW,sim_data.starchC6_ED_28,sim_data.starchC6_EN_28];
catch
    %Do nothing
end
