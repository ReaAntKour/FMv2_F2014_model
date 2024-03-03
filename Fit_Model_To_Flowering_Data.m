function [param]=FitModelToMattData_Fl_only(clock_dynamics_model_i)
% This function fits the model to the lab's flowering time data

if nargin<1
	clock_dynamics_model_i=4; % F2014 is default
end

% Include model folders folders in path
addpath('PIF_CO_FT_model')
addpath('circadian_module')
addpath('parameters')
addpath('phenology')
addpath('Fitting_Data')

% Import experimental data
Fl=readtable('Fitting_data_flowering.csv');

% input envirnmental conditions
temp = 22;
rise = 0;
load('weather.mat','weather')
hour=weather(:,1);
T=temp*weather(:,2);
sunrise=rise*weather(:,3);

% Set photoperiod sets to run
PhotAll=[8,12,16];
PhotCell=cell(1,2);
[~,PhotCell{1}]=ismember(Fl.Photoperiod,PhotAll);

% Set model related parameters
run_phenology_model=1;
options = struct();
options.entrain = 12; % entrain model at 12/12
options.genotype = {'wt'}; % not mutant
options.temperature = 22; % temperature (oC) (only 22 or 27oC accepted)

if clock_dynamics_model_i==2
	clock_parameters = load_P2011_parameters(options.genotype);
	clock_dynamics = @P2011_dynamics_COP1;
	options.y0=[1.0151 0.956 0.0755 0.0041 0.506 0.0977 0.0238 0.0731 0.0697 0.0196 0.0435 0.2505 0.0709 0.1017 0.0658 0.4016 0.1167 0.1012 0.207 0.0788 0.3102 0.0553 0.2991 0.1503 0.0286 0.65 0.2566 0.1012 0.576 0.3269];
elseif clock_dynamics_model_i==4
	clock_parameters = load_F2014_parameters(options.genotype);
	clock_dynamics = @F2014_dynamics_COP1;
	options.y0=0.1*ones(1,35);
end

% load parameters other than clock and PIF parameters
load('parameter.mat','parameter')
p=parameter; 
%specify Col accession for phenology model threshold
flowering_thresh_geno = 2; % threshold for Col-0 is p(17)
options.hypocotyl_parameters=struct();
load('Hypocotyl_parameters_P2011','paramHy')
options.hypocotyl_parameters.a1 = paramHy(1);
options.hypocotyl_parameters.a2 = paramHy(2);
options.hypocotyl_parameters.a3 = paramHy(3);

% Define initial parameters
x0=[p(17)+700];

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

	% set flowering threshold very high to stop simulation due to N_max_days
	% and so to have record of CumPhenThrm for all 90 days
	p(17)=1e7;

	[~,sim_data] = simulate_FM_simp(hour,T,sunrise,sunset,clock_parameters,clock_dynamics,clock_dynamics_model_i,p,run_phenology_model,options);

	CumPhenThrm{ip} = sim_data.CumPhenThrm;
	
	plot(1:90,sim_data.CumPhenThrm)
	hold on
end

%% RUN optimiser
[param, sse] = fminsearch(@(param) costfunction(param, CumPhenThrm, PhotAll, PhotCell{1}, Fl, p,flowering_thresh_geno), x0);

param
sse

Plotting(param, options, Fl, hour,T,sunrise,weather,clock_parameters,clock_dynamics,clock_dynamics_model_i,p,run_phenology_model)

rmpath('PIF_CO_FT_model')
rmpath('circadian_module')
rmpath('parameters')
rmpath('phenology')
rmpath('Fitting_Data')
end

%% Plotting
function Plotting(param, options, Fl, hour,T,sunrise,weather,clock_parameters,clock_dynamics,clock_dynamics_model_i,p,run_phenology_model)
	figure('Name','Physiology parameter fitting')

	% Read parameters
	p(17) = param;

	%% Run core model
	% Initialise arrays
	Phot = 5:1:24;
	nPh=length(Phot);
	Days_to_flower = zeros(nPh,1);

	% run for each photoperiod
	for ip=1:length(Phot)
		% set the conditions in this case
		set = Phot(ip); % photoperiod (Hr)
		options.photoperiod = set;
		sunset=set*weather(:,4); 
	
		[~,sim_data] = simulate_FM_simp(hour,T,sunrise,sunset,clock_parameters,clock_dynamics,clock_dynamics_model_i,p,run_phenology_model,options);
		
		% estimate days to flower
		Days_to_flower(ip) = sim_data.Fl;
	end
	
	% Flowering plot
	plot(Phot,Days_to_flower,'--k')
	hold on
	% Plot Flowering data
	errorbar([8,12,16],[mean(Fl.Bolting(Fl.Photoperiod==8)),mean(Fl.Bolting(Fl.Photoperiod==12)),mean(Fl.Bolting(Fl.Photoperiod==16))],[std(Fl.Bolting(Fl.Photoperiod==8))/sqrt(sum(Fl.Photoperiod==8)),std(Fl.Bolting(Fl.Photoperiod==12))/sqrt(sum(Fl.Photoperiod==12)),std(Fl.Bolting(Fl.Photoperiod==16))/sqrt(sum(Fl.Photoperiod==16))],'sk','MarkerFaceColor','k')
	title('Flowering time')
	xlabel('Photoperiod')
	box on
	ylim([0 100])
	xlim([5.5,24])
	xlabel('Photoperiod (Hrs)')
	ylabel('Days to Flower')
end


%% Cost calculation Fuction
function cost=costfunction(param, CumPhenThrm, PhotAll, PhotCell, Fl, p, flowering_thresh_geno)
	% Set param values
	p(17) = param;
	
	% Initialise arrays
	Phot = PhotAll;
	nPh=length(Phot);
	Days_to_flower = zeros(nPh,1);

	% run for each photoperiod
	for ip=1:length(PhotAll)
		% reset initial parameters
		day_idx = 1;
		has_flowered = false;
		N_max_days = 90;
		while day_idx <= N_max_days && ~(has_flowered)
    		% check if plant has flowered
			has_flowered = flowering_threshold_test(CumPhenThrm{ip}(day_idx),flowering_thresh_geno,p);
    		%Update day index
    		day_idx = day_idx+1;
		end
		Days_to_flower(ip) = day_idx - 1;
	end

	% Cost Flowering
	cost=sum((Fl.Bolting-Days_to_flower(PhotCell)).^2)
end