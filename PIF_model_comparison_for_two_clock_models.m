% This script simulates and plots the PIF model outputs when using for the
% clock model either the P2011 model or the F2014 model in the given 
% photoperiod conditions. It runs one day at a time until flowering.

clear
run_phenology_model=1;

%% input envirnmental conditions
temp = 22;
rise = 0;
set = 12;
load('weather.mat')
hour=weather(:,1);
T=temp*weather(:,2);
sunrise=rise*weather(:,3);
sunset=set*weather(:,4); 
options=struct();
options.entrain = 12; % entrain model at 12/12
options.photoperiod=set-rise; % only used for circaplot, elsewhere the sunset-sunrise is used for each day.

% Define Hypocotyl parameters
options.hypocotyl_parameters=struct();
load('Hypocotyl_parameters')
options.hypocotyl_parameters.a1 = paramHy(1);%0.9;%10/24;
options.hypocotyl_parameters.a2 = paramHy(2);%0.03103;%0;
options.hypocotyl_parameters.a3 = paramHy(3);%0.8;%100;

% load parameters other than clock and PIF parameters
load('parameter.mat')
p=parameter;

%Specify the genotypes to run
mutant_genotypes={{'wt'},{'elf3'}};%,{'YHB'},{'elf3','YHB'}};
nG=length(mutant_genotypes);

% Name the models to run
Models = {'P2011', 'F2014'};

% Include model folders and plotting folders in path
addpath('C:\Users\ra134k\OneDrive - University of Glasgow\Projects\Matt\Seaton 2015\published_model\plotting_tools')
addpath('PIF_CO_FT_model')

for ig=1:nG
	% set the genotype
	options.genotype = mutant_genotypes{ig};
	figure('Name',string(join(options.genotype)))
	for clock_dynamics_model_i=1:2
		%% Clock model
		if clock_dynamics_model_i==1
			clock_parameters = load_P2011_parameters(options.genotype);
			clock_dynamics = @P2011_dynamics;
% 			clock_dynamics_wrapper = @wrap_P2011_model_dynamics;
			options.y0=[1.0151 0.956 0.0755 0.0041 0.506 0.0977 0.0238 0.0731 0.0697 0.0196 0.0435 0.2505 0.0709 0.1017 0.0658 0.4016 0.1167 0.1012 0.207 0.0788 0.3102 0.0553 0.2991 0.1503 0.0286 0.65 0.2566 0.1012 0.576 0.3269];
		else
			clock_parameters = load_F2014_parameters(options.genotype);
			clock_dynamics = @F2014_dynamics;
% 			clock_dynamics_wrapper = @wrap_F2014_model_dynamics;
			options.y0=0.1*ones(1,35);
		end

		[~,sim_data] = simulate_FM_simp(hour,T,sunrise,sunset,clock_parameters,clock_dynamics,clock_dynamics_model_i,p,run_phenology_model,options);

		make_and_format_plots(sim_data,Models{clock_dynamics_model_i},options)
	end
end

rmpath('C:\Users\ra134k\OneDrive - University of Glasgow\Projects\Matt\Seaton 2015\published_model\plotting_tools')
rmpath('PIF_CO_FT_model')

function make_and_format_plots(output,clock_dynamics_model,options)
	% PIF model - plot the results for each species on separate panels

	species_names_all = {'PIF4m', 'PIF5m', 'PIF', 'phyB', 'PR', 'INT', 'IAA29m', 'ATHB2m', 'CDF1m', 'FKF1m', 'FKF1', 'CDF1', 'COm', 'CO', 'FTm', 'InducedC1', 'InducedC2', 'RepressedC1'};
	Sidx = [3,4,7,8,13:15];
	ylimit_T=[150, 1.2, 2, 6, 1, 5, 2];
	nS = length(Sidx);
	species_names = [species_names_all(Sidx),{'Hypocotyl length'},{'Days to Flower'}];

	for j = 1:nS
    	subplot(1,nS,j)
		plot(output.PIF_output.T,output.PIF_output.Y(:,Sidx(j)),'DisplayName',clock_dynamics_model)
    	hold on
		box on
		xlim([0,24])
		ylim([0,ylimit_T(j)])
		xlabel('Time (ZT Hrs)')
		ylabel('Relative Expression')
		circaplot([],[],[0,options.photoperiod],['w','k'],24)
    	v = axis;
		text(v(1)+0.5,v(4)*0.8,species_names{j},'FontAngle','italic')
	end
end