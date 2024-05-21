function [DayPhenThrm,FT_module_state,PIF_output] = phen(T,t,sunrise,sunset,geno,clock_dynamics_model_i,clock_output,FT_module_state,p,options) % Code edited here to add option of multiple models, flowring_genotype, hypocotyl_parameters and output gene levels
% Code edited by Rea L Antoniou-Kourounioti (University of Glasgow) as indicated
%% phenology model to determine flowering time
% Inputs:
%   T - Temperature (celcius)
%   t - time (hours)
%   sunrise - sunrise (hours)
%   sunset - sunset (hours)
%   geno - genotype index
%   clock_parameters - vector of clock parameters
%   clock_output - output timeseries for the clock model (concentrations of model components over time)
%   FT_module_state - state vector of the FT model (concentrations of model components at initial timepoint)
%   p - parameter vector
%
% Outputs:
%   DayPhenThrm - Accumulated modified thermal time units
%   FT_module_state - state vector of the FT model at the end of the timestep
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

% Fixed parameters
%_________________

Tb = p(1);
Tvmin = p(2);
Tvmax = p(3);
Vsat = p(4);
v = p(5);
sigma = p(6);
m = p(7);
Dld = p(8);
CSDL = p(9);
CLDL = p(10);

if  geno ==1
    
    Fb = p(11);
    Dsd = p(12);
    Threshold = p(13);%2907;%3259
    Night = p(14);
    Phot_a = p(73);
    Phot_b = p(74);
    Phot_c = p(75);
    Phot_n = p(76);
    
else

    Fb = p(15);
    Dsd = p(16);
    Threshold = p(17);%3212;%4218
    Night = p(18);
    Phot_a = p(77);
    Phot_b = p(78);
    Phot_c = p(79);
    Phot_n = p(80);

end

%Overwrite photoperiod response parameters
Phot_a = 1;
Phot_b = -0.4016;
Phot_c = 2.7479;
Phot_n = 3;

%Calculation begins
%__________________
Thrm = zeros(24,1);

[dailyFTarea,FT_module_state,PIF_output] = simulate_PIF_CO_FT_model(sunrise,sunset,clock_output,clock_dynamics_model_i,FT_module_state,T,options); % Code edited here to add option of multiple models, temperature, flowering_genotype, hypocotyl_parameters and output gene levels

hour = 1:24;

for i = 1:24

    %Calculating thermal component
    %_____________________________

    if      sunrise>=hour(i) || sunset<=hour(i)-1
            fraction_light(i)=0;
    elseif  sunrise<=hour(i)-1 && sunset>hour(i)
            fraction_light(i)=1;
    elseif  sunrise>=hour(i)-1
            fraction_light(i)=hour(i)-sunrise;
    else
            fraction_light(i)=sunset-hour(i)+1;
    end



    if      fraction_light(i)==0
            Thrm(i)=Night*max(T(t+i-1)-Tb,0);
    elseif  fraction_light(i)==1    
            Thrm(i)=max(T(t+i-1)-Tb,0);
    else 
            Thrm(i)=max(0,(T(t+i-1)-Tb)*fraction_light(i)) + Night*max(0,(T(t+i-1)-Tb)*(1-fraction_light(i)));
    end




    %Calculating photoperiod component
    %_________________________________

%     dl(i) = sunset(i) - sunrise(i);

%     if  i == 1
% 
%         [FTarea24,yo] = link(dl(1),sunrise(1)); %initialise
%         [FTarea1,yo] = sublink(hour(i),dl(i),sunrise(i),yo);
%         FTarea24 = [FTarea24(2:end) FTarea1];
%         dailyFTarea(i) = sum(FTarea24);                                       
%     else
%         [FTarea1,yo] = sublink(hour(i),dl(i),sunrise(i),yo);
%         FTarea24 = [FTarea24(2:end) FTarea1];
%         dailyFTarea(i) = sum(FTarea24);
%     end

%     dailyFTarea(i) = 0;
    Phot(i) = Phot_a + Phot_b*Phot_c^Phot_n/(Phot_c^Phot_n+dailyFTarea^Phot_n);
%     Phot(i)



    %Calculating vernalization effectiveness
    %_______________________________________


    if	T(t+i-1)>=Tvmin && T(t+i-1)<=Tvmax

        Ve(i)=exp(v)*(T(t+i-1)-Tvmin)^m*(Tvmax-T(t+i-1))^sigma*1;
    else	
        Ve(i)=0;
    end



    %Calculating cumulative vernalization hours
    %__________________________________________


    Vh(i)=sum(Ve(1:i));




    %Calculating vernalization fraction
    %__________________________________

    if	Vh(i)<=Vsat

        Vern(i) = Fb + Vh(i)*(1 - Fb)/Vsat;
    else
        Vern(i)= 1;	
    end



    %Calculating modified photothermal unit
    %______________________________________

    mptu(i) = Vern(i)*Phot(i)*Thrm(i);

end


%Output:
%-------
DayPhenThrm=sum(mptu);
