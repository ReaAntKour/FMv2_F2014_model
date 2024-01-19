% load('simulation_data_file','sim_data')
%sim_data contains all output simulation data
temp = 20.5;
rise = 0;
set = 12;
co2 = 42;
light = 145;

%Calling for meteorological data
%_______________________________  

%hour,T,sunrise,sunset,CO2,PAR
load('weather.mat')

hour=weather(:,1);
T=temp*weather(:,2);
sunrise=rise*weather(:,3);
sunset=set*weather(:,4); 
CO2=co2*weather(:,5); %CO2 partial pressure (Pa)
PAR=light*weather(:,6); %total absorbed PAR per unit leaf area (micromol m-2 s-1)
Photoperiod = sunset(end) - sunrise(end);


Measured_day = [22 28 29 38];
Measured_hour = [0 12 24];

run_phenology_model = 0;

for experiment_genotype_index = [1,2,3,4]

    %Specifying the genotype for the clock and starch models
    switch experiment_genotype_index
        case 1
            w = 0.91; %water content
            mf_use = 0.7; %malate+fumarate turnover
            clock_genotype = {''};
            starch_genotype = {''};
        case 2
            w = 0.89; %water content
            mf_use = 0.7; %malate+fumarate turnover
            clock_genotype = {''};
            starch_genotype = {'lsf1'};
        case 3
            w = 0.89; %water content
            mf_use = 0.7; %malate+fumarate turnover
            clock_genotype = {'prr9','prr7'};
            starch_genotype = {''};
        case 4
            w = 0.89; %water content
            mf_use = 0.25; %malate+fumarate turnover
            clock_genotype = {'prr9','prr7'};
            starch_genotype = {''};
    end
    
    % w = 0.9;
    d = 1-w; %dry matter


    load('parameter.mat')

    p=parameter;

    switch experiment_genotype_index
        case 1
            Measured_biomass = [0.034 0.1262 0.2192 0.95014]; %Col
            Measured_sd_biomass = [0.004092676	0.017887426	0.032882518	0.14407886]; %Col
            Measured_starch = [2.41 42.70 4.55]; %Col
            Measured_OA = [4.85 17.87 5.38]; %Col
            sd_starch = [0.35 2.37 0.31]; %Col
            sd_OA = [0.51 1.24 0.71]; %Col

        case 2
            Measured_biomass = [0.01866	0.08546	0.0968	0.48772]; %lsf
            Measured_sd_biomass = [0.003036939	0.007267255	0.012679511	0.046257616]; %lsf
            Measured_starch = [45.28 90.19 44.25]; %lsf
            Measured_OA = [6.62 24.18 7.85]; %lsf    
            sd_starch = [1.05 4.61 2.08]; %lsf
            sd_OA = [1.18 1.07 0.27]; %lsf

        case 3
            Measured_biomass = [0.02102	0.09238	0.10086	0.5743]; %prr7/prr9
            Measured_sd_biomass = [0.001413153	0.013399515	0.008260327	0.061144337]; %prr7/prr9
            Measured_starch = [9.93 67.31 23.77]; %prr7/prr9
            Measured_OA = [18.65 41.80 25.31]; %prr7/prr9    
            sd_starch = [0.08 0.01 0.14]; %prr7/prr9
            sd_OA = [2.14 2.69 4.04]; %prr7/prr9

        case 4
            Measured_biomass = [0.02102	0.09238	0.10086	0.5743]; %prr7/prr9
            Measured_sd_biomass = [0.001413153	0.013399515	0.008260327	0.061144337]; %prr7/prr9
            Measured_starch = [9.93 67.31 23.77]; %prr7/prr9
            Measured_OA = [18.65 41.80 25.31]; %prr7/prr9    
            sd_starch = [0.08 0.01 0.14]; %prr7/prr9
            sd_OA = [2.14 2.69 4.04]; %prr7/prr9
    end        

    clock_parameters = P2011_parameter_call(clock_genotype);
    starch_parameters = starch_parameter_call(starch_genotype);

    [output,sim_data] = simulate_FM(hour,T,sunrise,sunset,CO2,PAR,Photoperiod,clock_parameters,starch_parameters,p,d,mf_use,run_phenology_model);
    
    Estimate_biomass = sim_data.Total_shoot(Measured_day*24)/d;
    cvRMSE = (sum((Estimate_biomass - Measured_biomass).^2)/numel(Measured_biomass))^0.5/...
                    (mean(Measured_biomass))

    nRMSE = (sum((Estimate_biomass - Measured_biomass).^2)/numel(Measured_biomass))^0.5/...
                    (max(Measured_biomass) - min(Measured_biomass));

    Rsquare = corrcoef(Measured_biomass,Estimate_biomass).^2;
    Rsquare = Rsquare(2);

    figure(1)
    hold on
    plot(0:24,sim_data.starchC6(28*24:29*24))
    errorbar(Measured_hour,Measured_starch,sd_starch,'o')
    axis([0 24 0 100])
    title('Starch profile')
    xlabel('Hour')
    ylabel('Starch level (micromol C6 per g FW)')
    hold off

    figure(2)
    hold on
    plot(0:24,sim_data.malatefumarateC4(28*24:29*24))
    errorbar(Measured_hour,Measured_OA,sd_OA,'o')
    axis([0 24 0 50])
    title('Malate + Fumarate profile')
    xlabel('Hour')
    ylabel('Malate + Fumarate level (micromol C4 per g FW)')
    hold off
    % 
    % figure
    % plot(0:24,sim_data.sugar(28*24:29*24))
    % title('Sucrose profile')
    % xlabel('Hour')
    % ylabel('Sucrose level (micromol per g FW)')
    % %axis([0 24 0 6])
    % 
    % figure
    % plot(1/24:1/24:40,sim_data.sugar)
    % title('Sucrose profile')
    % xlabel('DAS')
    % ylabel('Sucrose level (micromol per g FW)')
    % 
    figure(3)
    hold on
    errorbar(Measured_day,Measured_biomass,Measured_sd_biomass,'o')
    plot(1/24:1/24:40,sim_data.FW)
    title('Rosette FW')
    xlabel('DAS')
    ylabel('g FW')
    hold off
end