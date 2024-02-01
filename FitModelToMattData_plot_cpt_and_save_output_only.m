function [param]=FitModelToMattData_plot_cpt_and_save_output_only()
% This function fits flowering time and hypocotyl data

% Import experimental data
Fl=readtable('C:\Users\ra134k\OneDrive - University of Glasgow\Projects\Matt\Seaton 2015\Matt_flowering_data.csv');
Hy=readtable('C:\Users\ra134k\OneDrive - University of Glasgow\Projects\Matt\Seaton 2015\Matt_hypocotyl_data.csv');

% Include model folders and plotting folders in path
addpath('C:\Users\ra134k\OneDrive - University of Glasgow\Projects\Matt\Seaton 2015\published_model\plotting_tools')
addpath('PIF_CO_FT_model')

% input envirnmental conditions
temp = 22;
rise = 0;
load('weather.mat','weather')
hour=weather(:,1);
T=temp*weather(:,2);
sunrise=rise*weather(:,3);

% Set photoperiod sets to run
PhotAll=[0,8,12,16];
PhotCell=cell(1,2);%{[2,3,4],[1,2,4]};% for plotting: [5:1:24,16];
[~,PhotCell{1}]=ismember(Fl.Photoperiod,PhotAll);
[~,PhotCell{2}]=ismember(Hy.Photoperiod,PhotAll);

% Set model related parameters
run_phenology_model=1;
options = struct();
options.entrain = 12; % entrain model at 12/12
options.genotype = {'wt'}; % not mutant
options.temperature = 22; % temperature (oC) (only 22 or 27oC accepted)
clock_parameters = load_F2014_parameters(options.genotype);
clock_dynamics = @F2014_dynamics;
clock_dynamics_model_i=2;
options.y0=0.1*ones(1,35);

% load parameters other than clock and PIF parameters
load('parameter.mat','parameter')
p=parameter;
load('Flowering_parameter_F2014','paramFl')
p(17) = paramFl;% CumPhenThrm flowering threshold for Col-0 is p(17)
options.hypocotyl_parameters=struct();
load('Hypocotyl_parameters_P2011','paramHy')
options.hypocotyl_parameters.a1 = paramHy(1);
options.hypocotyl_parameters.a2 = paramHy(2);
options.hypocotyl_parameters.a3 = paramHy(3);

% Initialise arrays
Phot = PhotAll;
nPh=length(Phot);
CumPhenThrm = cell(nPh,1);

% run for each photoperiod
figure
for ip=1:length(PhotAll)
	% set the conditions in this case
	set = Phot(ip); % photoperiod (Hr)
	options.photoperiod = set;
	sunset=set*weather(:,4); 

	[~,sim_data] = simulate_FM_simp(hour,T,sunrise,sunset,clock_parameters,clock_dynamics,clock_dynamics_model_i,p,run_phenology_model,options);

	CumPhenThrm{ip} = sim_data.CumPhenThrm;
	
	plot(1:90,sim_data.CumPhenThrm)
	hold on
end


load('Flowering_parameter_F2014','paramFl')
load('Hypocotyl_parameters_F2014','paramHy')
param=[paramHy,paramFl];
%sse

ModelFitting = Plotting(param, PhotAll, options, Fl, Hy, hour,T,sunrise,weather,clock_parameters,clock_dynamics,clock_dynamics_model_i,p,run_phenology_model);
writetable(ModelFitting,'ModelFitting.csv')
end

%% Plotting
function out=Plotting(param, PhotAll, options, Fl, Hy, hour,T,sunrise,weather,clock_parameters,clock_dynamics,clock_dynamics_model_i,p,run_phenology_model)
	figure('Name','Physiology parameter fitting')

	% Read parameters
	options.hypocotyl_parameters.a1 = param(1);
	options.hypocotyl_parameters.a2 = param(2);
	options.hypocotyl_parameters.a3 = param(3);
	p(17) = param(4);
	
	%% Run core model
	% Initialise arrays
	Phot = PhotAll;
	nPh=length(Phot);
	Hyp_length = zeros(nPh,1);
% 	FT_area = zeros(nPh,1);
	Days_to_flower = zeros(nPh,1);

	% run for each photoperiod
	for ip=1:length(PhotAll)
		% set the conditions in this case
		set = Phot(ip); % photoperiod (Hr)
		options.photoperiod = set;
		sunset=set*weather(:,4); 
	
		[~,sim_data] = simulate_FM_simp(hour,T,sunrise,sunset,clock_parameters,clock_dynamics,clock_dynamics_model_i,p,run_phenology_model,options);
	
		Hyp_length(ip)=sim_data.Hyp;
	
		% estimate days to flower based on area under FTm
% 		FT_area(ip) = sim_data.PIF_output.FTarea;
		Days_to_flower(ip) = sim_data.Fl;
	end

	% Hypocotyl plot
	subplot(1,2,1)
	plot(PhotAll,Hyp_length,'--k')
	hold on
	% Plot Hypocotyl data
% 	plot([0,8,16],[mean(Hy.Hypocotyl_length(Hy.Photoperiod==0)),mean(Hy.Hypocotyl_length(Hy.Photoperiod==8)),mean(Hy.Hypocotyl_length(Hy.Photoperiod==16))])
	errorbar([0,8,16],[mean(Hy.Hypocotyl_length(Hy.Photoperiod==0)),mean(Hy.Hypocotyl_length(Hy.Photoperiod==8)),mean(Hy.Hypocotyl_length(Hy.Photoperiod==16))],[std(Hy.Hypocotyl_length(Hy.Photoperiod==0))/sqrt(sum(Hy.Photoperiod==0)),std(Hy.Hypocotyl_length(Hy.Photoperiod==8))/sqrt(sum(Hy.Photoperiod==8)),std(Hy.Hypocotyl_length(Hy.Photoperiod==16))/sqrt(sum(Hy.Photoperiod==16))],'sk','MarkerFaceColor','k')
	title('Hypocotyl Length')
	xlabel('Photoperiod')
	
	% Flowering plot
	subplot(1,2,2)
	plot(PhotAll,Days_to_flower,'--k')
	hold on
	% Plot Flowering data
% 	plot([8,12,16],[mean(Fl.Bolting(Fl.Photoperiod==8)),mean(Fl.Bolting(Fl.Photoperiod==12)),mean(Fl.Bolting(Fl.Photoperiod==16))])
	errorbar([8,12,16],[mean(Fl.Bolting(Fl.Photoperiod==8)),mean(Fl.Bolting(Fl.Photoperiod==12)),mean(Fl.Bolting(Fl.Photoperiod==16))],[std(Fl.Bolting(Fl.Photoperiod==8))/sqrt(sum(Fl.Photoperiod==8)),std(Fl.Bolting(Fl.Photoperiod==12))/sqrt(sum(Fl.Photoperiod==12)),std(Fl.Bolting(Fl.Photoperiod==16))/sqrt(sum(Fl.Photoperiod==16))],'sk','MarkerFaceColor','k')
	title('Flowering time')
	xlabel('Photoperiod')
	
	% Format plots
	species_names = {'Hypocotyl length','Days to Flower'};
	for j=1:2
		subplot(1,2,j)
		box on
		if j==1
			ylim([0 7])
			xlim([-0.5,24])
		else
			ylim([0 100])
			xlim([5.5,24])
		end
		xlabel('Photoperiod (Hrs)')
		ylabel(species_names{j})
	end

	
	photoperiod = PhotAll';
	Hypocotyl_length_model = Hyp_length;
	Days_to_flower_model = Days_to_flower;
	Hypocotyl_length_data_mean = [mean(Hy.Hypocotyl_length(Hy.Photoperiod==0)),mean(Hy.Hypocotyl_length(Hy.Photoperiod==8)),NaN,mean(Hy.Hypocotyl_length(Hy.Photoperiod==16))]';
	Days_to_flower_data_mean = [NaN,mean(Fl.Bolting(Fl.Photoperiod==8)),mean(Fl.Bolting(Fl.Photoperiod==12)),mean(Fl.Bolting(Fl.Photoperiod==16))]';
	out = table(photoperiod,Hypocotyl_length_model,Days_to_flower_model,Hypocotyl_length_data_mean,Days_to_flower_data_mean);
end
