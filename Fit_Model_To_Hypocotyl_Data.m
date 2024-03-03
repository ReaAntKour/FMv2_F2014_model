function [param]=FitModelToMattData_Hy_only(clock_dynamics_model_i)
% This function fits the model to the lab's hypocotyl data

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
Hy=readtable('Fitting_data_hypocotyl.csv');

% input envirnmental conditions
temp = 22;
rise = 0;
load('weather.mat','weather')
hour=weather(:,1);
T=temp*weather(:,2);
sunrise=rise*weather(:,3);

% Set photoperiod sets to run
PhotAll=[0,8,16];
PhotCell=cell(1,2);
[~,PhotCell{2}]=ismember(Hy.Photoperiod,PhotAll);

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

% load initial Hypocotyl parameters
options.hypocotyl_parameters=struct();
load('Hypocotyl_parameters_P2011','paramHy')
options.hypocotyl_parameters.a1 = paramHy(1);
options.hypocotyl_parameters.a2 = paramHy(2);
options.hypocotyl_parameters.a3 = paramHy(3);

% Define initial parameters
x0=paramHy;

% Initialise arrays
Phot = PhotAll;
nPh=length(Phot);
PIF_Output = cell(nPh,1);

% run for each photoperiod
for ip=1:length(PhotAll)
	% set the conditions in this case
	set = Phot(ip); % photoperiod (Hr)
	options.photoperiod = set;
	sunset=set*weather(:,4); 

	[~,sim_data] = simulate_FM_simp(hour,T,sunrise,sunset,clock_parameters,clock_dynamics,clock_dynamics_model_i,p,run_phenology_model,options);

	PIF_Output{ip}=sim_data.PIF_output;
end

%% RUN optimiser
[param, sse] = fminsearch(@(param) costfunction(param, PIF_Output, PhotCell{2}, Hy), x0);

param % P2011: 0.2527   -0.3820    0.6436
sse

Plotting(param, options, Hy, hour,T,sunrise,weather,clock_parameters,clock_dynamics,clock_dynamics_model_i,p,run_phenology_model)

rmpath('PIF_CO_FT_model')
rmpath('circadian_module')
rmpath('parameters')
rmpath('phenology')
rmpath('Fitting_Data')
end

%% Plotting
function Plotting(param, options, Hy, hour,T,sunrise,weather,clock_parameters,clock_dynamics,clock_dynamics_model_i,p,run_phenology_model)
	figure('Name','Physiology parameter fitting')

	% Read parameters
	options.hypocotyl_parameters.a1 = param(1);
	options.hypocotyl_parameters.a2 = param(2);
	options.hypocotyl_parameters.a3 = param(3);
	
	%% Run core model
	% Initialise arrays
	Phot = 0:1:24;
	nPh=length(Phot);
	Hyp_length = zeros(nPh,1);

	% run for each photoperiod
	for ip=1:nPh
		% set the conditions in this case
		set = Phot(ip); % photoperiod (Hr)
		options.photoperiod = set;
		sunset=set*weather(:,4); 
		% run model
		[~,sim_data] = simulate_FM_simp(hour,T,sunrise,sunset,clock_parameters,clock_dynamics,clock_dynamics_model_i,p,run_phenology_model,options);
		% save hypocotyl length for each photoperiod
		Hyp_length(ip)=sim_data.Hyp;
	end

	% Hypocotyl plot
	plot(Phot,Hyp_length,'--k')
	hold on
	% Plot Hypocotyl data
	errorbar([0,8,16],[mean(Hy.Hypocotyl_length(Hy.Photoperiod==0)),mean(Hy.Hypocotyl_length(Hy.Photoperiod==8)),mean(Hy.Hypocotyl_length(Hy.Photoperiod==16))],[std(Hy.Hypocotyl_length(Hy.Photoperiod==0))/sqrt(sum(Hy.Photoperiod==0)),std(Hy.Hypocotyl_length(Hy.Photoperiod==8))/sqrt(sum(Hy.Photoperiod==8)),std(Hy.Hypocotyl_length(Hy.Photoperiod==16))/sqrt(sum(Hy.Photoperiod==16))],'sk','MarkerFaceColor','k')
	title('Hypocotyl Length')
	xlabel('Photoperiod')
	box on
	ylim([0 7])
	xlim([-0.5,24])
	xlabel('Photoperiod (Hrs)')
	ylabel('Hypocotyl length')
end


%% Cost calculation Fuction
function cost=costfunction(param, PIF_Output, PhotCell, Hy)
	nPh=length(PhotCell);
	
	% Set param values
	a1 = param(1);
	a2 = param(2);
	a3 = param(3);

	
	% Initialise arrays
	Hyp_length = zeros(nPh,1);
	
	for ip=1:nPh
		T=PIF_Output{PhotCell(ip)}.T;
		Y=PIF_Output{PhotCell(ip)}.Y;
		
		% calculate area under ATHB2m as proxy for hypocotyl length
		Hyp_length(ip)=a1*trapz(T,min(Y(:,8),a3)-a2);
	end

	param
	cost=sum((Hy.Hypocotyl_length-Hyp_length).^2)
end