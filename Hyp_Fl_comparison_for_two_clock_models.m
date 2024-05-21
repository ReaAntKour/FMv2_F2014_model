% This script simulates and plots Flowering time and  Hypocotyl length  
% using for the clock model either the P2011 model or the F2014 model, with
% the Global or COP1-only YHB mutant version, in various photoperiod
% conditions. The script produces plots for each mutant and .csv files with
% the results from each of the four model variants.
% 
%   Copyright 2024 Rea L Antoniou-Kourounioti and The University of Glasgow
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
% 
% 

function Hyp_Fl_comparison_for_two_clock_models(YHB) % clear % 
run_phenology_model=1;

% Include model folders in path
addpath('PIF_CO_FT_model')
addpath('circadian_module')
addpath('parameters')
addpath('phenology')

%% input envirnmental conditions
temp = 22;
rise = 0;
load('weather.mat','weather')
hour=weather(:,1);
T=temp*weather(:,2);
sunrise=rise*weather(:,3);
options=struct();
options.entrain = 12; % entrain model at 12/12

% load parameters other than clock and PIF parameters
load('parameter.mat','parameter')
p=parameter;
load('Flowering_parameter_F2014','paramFl')
p(17) = paramFl;% CumPhenThrm flowering threshold for Col-0 is p(17)

% Define Hypocotyl parameters
options.hypocotyl_parameters=struct();
% further defined according to model

%Specify the genotypes to run
mutant_genotypes={{'wt'},{'YHB'},{'elf3'},{'elf3','YHB'}};
nG=length(mutant_genotypes);

% Name the models to run
Models = {'P2011_Red', 'P2011_COP1', 'F2014_Red', 'F2014_COP1'};

% set the parameter set for F2014 model
options.paramSet=1;
options.YHB=YHB; % 0.5 (33%)Y, 1 (50%), 3 (75%)Y, 9 (90%)Y or 999999 (99.9999%)

% set photoperiods to run
Phot=0:1:20;
nPh=length(Phot);
Hyp_length = zeros(nG,nPh,length(Models));
FT_area = zeros(nG,nPh,length(Models));
Days_to_flower = zeros(nG,nPh,length(Models));

Model_output_to_file=struct('wt',struct(),'YHB',struct(),'elf3',struct(),'elf3YHB',struct());
for ig=1:nG
	% set the genotype
	options.genotype = mutant_genotypes{ig};
	Model_output_to_file.(string(join(options.genotype,'')))=struct();
% 	figure('Name',string(join(options.genotype)))
	for clock_dynamics_model_i=1:length(Models)
		% set the clock model
		if clock_dynamics_model_i<3
			clock_parameters = load_P2011_parameters(options.genotype,options.YHB);
			if clock_dynamics_model_i==1
				clock_dynamics = @P2011_dynamics_Red;
			elseif clock_dynamics_model_i==2
				clock_dynamics = @P2011_dynamics_COP1;
			end
			options.y0=[1.0151 0.956 0.0755 0.0041 0.506 0.0977 0.0238 0.0731 0.0697 0.0196 0.0435 0.2505 0.0709 0.1017 0.0658 0.4016 0.1167 0.1012 0.207 0.0788 0.3102 0.0553 0.2991 0.1503 0.0286 0.65 0.2566 0.1012 0.576 0.3269];
			load('Hypocotyl_parameters_P2011','paramHy')
			options.hypocotyl_parameters.a1 = paramHy(1);
			options.hypocotyl_parameters.a2 = paramHy(2);
			options.hypocotyl_parameters.a3 = paramHy(3);
		else
			clock_parameters = load_F2014_parameters(options.genotype,options.paramSet,options.YHB);
			if clock_dynamics_model_i==3
				clock_dynamics = @F2014_dynamics_Red;
			elseif clock_dynamics_model_i==4
				clock_dynamics = @F2014_dynamics_COP1;
			end
			options.y0=0.1*ones(1,35);
			load('Hypocotyl_parameters_F2014','paramHy')
			options.hypocotyl_parameters.a1 = paramHy(1);
			options.hypocotyl_parameters.a2 = paramHy(2);
			options.hypocotyl_parameters.a3 = paramHy(3);
		end

		for ip=1:nPh
			% set the photoperiod
			set = Phot(ip); % photoperiod (Hr)
			options.photoperiod = set;
			sunset=set*weather(:,4); 

			[~,sim_data] = simulate_FM_simp(hour,T,sunrise,sunset,clock_parameters,clock_dynamics,clock_dynamics_model_i,p,run_phenology_model,options);

			Hyp_length(ig,ip,clock_dynamics_model_i)=sim_data.Hyp;
		
			% estimate days to flower based on area under FTm
			FT_area(ig,ip,clock_dynamics_model_i) = sim_data.PIF_output.FTarea;
			Days_to_flower(ig,ip,clock_dynamics_model_i) = sim_data.Fl;
		end
		photoperiod = Phot';
		Hypocotyl_length_model = Hyp_length(ig,:,clock_dynamics_model_i)';
		Days_to_flower_model = Days_to_flower(ig,:,clock_dynamics_model_i)';
		Model_output_to_file.(string(join(options.genotype,''))).(Models{clock_dynamics_model_i}) = table(photoperiod,Hypocotyl_length_model,Days_to_flower_model);
% 		make_and_format_plots(Phot,Hyp_length(ig,:,:),FT_area(ig,:,:),Days_to_flower(ig,:,:),Models,clock_dynamics_model_i)
	end
end

%% Save output
sz = [21,9];
varTypes = {'double','double','double','double','double','double','double','double','double'};
varNames=["photoperiod";"HypocotylLength_wt";"HypocotylLength_YHB";"HypocotylLength_elf3";"HypocotylLength_elf3YHB";"DaysToFlower_wt";"DaysToFlower_YHB";"DaysToFlower_elf3";"DaysToFlower_elf3YHB"];
for clock_dynamics_model_i=1:length(Models)
	ModelHypFlMut=table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
	ModelHypFlMut.photoperiod = photoperiod;
	for ig=1:nG
		% set the genotype
		options.genotype = mutant_genotypes{ig};
		ModelHypFlMut.(string(join([{'HypocotylLength'},join(options.genotype,'')],'_'))) = Model_output_to_file.(string(join(options.genotype,''))).(Models{clock_dynamics_model_i}).Hypocotyl_length_model;
		ModelHypFlMut.(string(join([{'DaysToFlower'},join(options.genotype,'')],'_'))) = Model_output_to_file.(string(join(options.genotype,''))).(Models{clock_dynamics_model_i}).Days_to_flower_model;
	end
	writetable(ModelHypFlMut,['output\',int2str(floor(100*YHB/(YHB+1))),'\ModelHypFlMut_',Models{clock_dynamics_model_i},'.csv'])
end

rmpath('PIF_CO_FT_model')
rmpath('circadian_module')
rmpath('parameters')
rmpath('phenology')
end

function make_and_format_plots(Phot,Hyp_length,FT_area,Days_to_flower,Models,clock_dynamics_model_i)
	% Plot Flowering time and Hypocotyl length for different photoperiods
	subplot(1,3,1)
	plot(Phot,Hyp_length(1,:,clock_dynamics_model_i),'DisplayName',Models{clock_dynamics_model_i})
    hold on
	box on
	xlim([0,20])
	ylim([0 7])
   	xlabel('Photoperiod (Hrs)')
   	ylabel('Hypocotyl length')

	subplot(1,3,2)
	plot(Phot,FT_area(1,:,clock_dynamics_model_i),'DisplayName',Models{clock_dynamics_model_i})
    hold on
	box on
	xlim([0,20])
	ylim([0 15])
   	xlabel('Photoperiod (Hrs)')
   	ylabel('Daily FT area')

	subplot(1,3,3)
	plot(Phot,Days_to_flower(1,:,clock_dynamics_model_i),'DisplayName',Models{clock_dynamics_model_i})
    hold on
	box on
	xlim([0,20])
	ylim([0 100])
   	xlabel('Photoperiod (Hrs)')
   	ylabel('Days to Flower')
end