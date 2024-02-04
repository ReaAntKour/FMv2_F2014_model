% This script simulates and plots Flowering time and  Hypocotyl length  
% using for the clock model either the P2011 model or the F2014 model in  
% various photoperiod conditions.

clear
run_phenology_model=1;

%% input envirnmental conditions
temp = 22;
rise = 0;
load('weather.mat')
hour=weather(:,1);
T=temp*weather(:,2);
sunrise=rise*weather(:,3);
options=struct();
options.entrain = 12; % entrain model at 12/12

% load parameters other than clock and PIF parameters
load('parameter.mat')
p=parameter;
load('Flowering_parameter_F2014')
p(17) = paramFl;% CumPhenThrm flowering threshold for Col-0 is p(17)

% Define Hypocotyl parameters
options.hypocotyl_parameters=struct();
% further defined according to model

%Specify the genotypes to run
mutant_genotypes={{'wt'},{'YHB'},{'elf3'},{'elf3','YHB'}};
nG=length(mutant_genotypes);

% Name the models to run
Models = {'P2011', 'F2014'};

% set photoperiods to run
Phot=0:1:20;
nPh=length(Phot);
Hyp_length = zeros(nG,nPh,2);
FT_area = zeros(nG,nPh,2);
Days_to_flower = zeros(nG,nPh,2);

% estimate hypocotyl parameters from our data
%paramHy=FitHypModelToData();
%save('Hypocotyl_parameters.mat','paramHy','paramHy')

% Include model folders and plotting folders in path
addpath('C:\Users\ra134k\OneDrive - University of Glasgow\Projects\Matt\Seaton 2015\published_model\plotting_tools')
addpath('PIF_CO_FT_model')

Model_output_to_file=struct('wt',struct(),'YHB',struct(),'elf3',struct(),'elf3YHB',struct());
for ig=1:nG
	% set the genotype
	options.genotype = mutant_genotypes{ig};
	Model_output_to_file.(string(join(options.genotype,'')))=struct();
	figure('Name',string(join(options.genotype)))
	for clock_dynamics_model_i=1:2
		% set the clock model
		if clock_dynamics_model_i==1
			clock_parameters = load_P2011_parameters(options.genotype);
			clock_dynamics = @P2011_dynamics;
% 			clock_dynamics_wrapper = @wrap_P2011_model_dynamics;
			options.y0=[1.0151 0.956 0.0755 0.0041 0.506 0.0977 0.0238 0.0731 0.0697 0.0196 0.0435 0.2505 0.0709 0.1017 0.0658 0.4016 0.1167 0.1012 0.207 0.0788 0.3102 0.0553 0.2991 0.1503 0.0286 0.65 0.2566 0.1012 0.576 0.3269];
			load('Hypocotyl_parameters_P2011')
			options.hypocotyl_parameters.a1 = paramHy(1);
			options.hypocotyl_parameters.a2 = paramHy(2);
			options.hypocotyl_parameters.a3 = paramHy(3);
		else
			clock_parameters = load_F2014_parameters(options.genotype);
			clock_dynamics = @F2014_dynamics;
% 			clock_dynamics_wrapper = @wrap_F2014_model_dynamics;
			options.y0=0.1*ones(1,35);
			load('Hypocotyl_parameters_F2014')
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
		make_and_format_plots(Phot,Hyp_length(ig,:,:),FT_area(ig,:,:),Days_to_flower(ig,:,:),Models,clock_dynamics_model_i)
	end
end

% %% Problematic code
% % combinations only since R2023a but in R2023b the writetable seems broken.
% % Fixed by running varNames in R2023b and saving output as .mat and then
% % running the rest of the code in R2022a where it works.
% sz = [21,9];
% varTypes = {'double','double','double','double','double','double','double','double','double'};
% genotypes = fieldnames(Model_output_to_file)';
% variables = {'HypocotylLength','DaysToFlower'};
% varNames = [{'photoperiod'};string(join(table2cell(combinations(variables,genotypes)),'_'))];
% for clock_dynamics_model_i=1:2
% 	ModelHypFlMut=table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
% 	ModelHypFlMut.photoperiod = photoperiod;
% 	for ig=1:nG
% 		% set the genotype
% 		options.genotype = mutant_genotypes{ig};
% 		ModelHypFlMut.(string(join([{'HypocotylLength'},join(options.genotype,'')],'_'))) = Model_output_to_file.(string(join(options.genotype,''))).(Models{clock_dynamics_model_i}).Hypocotyl_length_model;
% 		ModelHypFlMut.(string(join([{'DaysToFlower'},join(options.genotype,'')],'_'))) = Model_output_to_file.(string(join(options.genotype,''))).(Models{clock_dynamics_model_i}).Days_to_flower_model;
% 	end
% 	writetable(ModelHypFlMut,['ModelHypFlMut_',Models{clock_dynamics_model_i},'.csv'])
% end
% %% end of Problematic code

rmpath('C:\Users\ra134k\OneDrive - University of Glasgow\Projects\Matt\Seaton 2015\published_model\plotting_tools')
rmpath('PIF_CO_FT_model')

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