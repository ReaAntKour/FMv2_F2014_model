% This script simulates and plots the clock model outputs for the P2011
% and F2014 models, with the Global or COP1-only YHB mutant version, in the
% given photoperiod conditions. The script produces plots for each mutant
% and .csv files with the results from each of the four model variants.
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

% start fresh
clear

% Include model folders in path
addpath('PIF_CO_FT_model')
addpath('circadian_module')
addpath('parameters')
addpath('phenology')

% Indices of the clock species of interest (i.e. CCA1/LHY mRNA)
clock_species_names_all = {'LHYm','P','GIZ','GIE','LHYp','TOCm','PRR9p','PRR5m','PRR5p','GIp','TOCp','ZTL','EC','GIm','PRR9m','PRR7m','PRR7p','ELF4m','ELF4p','LHYpm','HY5','HFR1','ELF3m','ELF3pc','ELF3pn','COP1pnn','COP1pnd','LUXm','LUXp','COP1pc'};
clock_Sidx1 = [1,14,15,26]; % LHYm, GIm, PRR9m, COP1pnn
clock_nS = length(clock_Sidx1);
clock_species_names = [clock_species_names_all(clock_Sidx1)];
clock_Sidx2 = [1,29,6,25]; % LHYm, GIm, PRR9m, COP1pnn % Replace index 1 with 3 for CCA1
clock_ylimit_T=[1.5,4,1.2,1];

mutant_genotypes={{'wt'},{'YHB'},{'elf3'},{'elf3','YHB'}};
nG=length(mutant_genotypes);

% Name the models to run
Models = {'P2011_Red', 'P2011_COP1', 'F2014_Red', 'F2014_COP1'};

% set the parameter set for F2014 model
paramSet=1;
YHB=3; % 0.5 (33%), 1 (50%), 3 (75%), 9 (90%) or 999999(99.9999%)

% set the common conditions in this case
options = struct();
options.temperature = 22; % temperature (oC) (only 22 or 27oC accepted)
options.entrain = 12; % entrain model at 12/12

% Load light conditions into 'c' for common light function
c.period = 24;
c.phase = 0;
c.dawn = 0;
twilightPer=0.00005; %The duration of time between value of force in dark and value of force in light.
LightAmp = 0.3; %The light intensity.

Phot_All={12, 24, 0, 18};
Time_Phot_All={'Time_Since_Dawn','Time_In_Constant_Light','Time_In_Constant_Dark','Time_Since_Dawn'};

for phot_i=1:length(Phot_All)
	options.photoperiod = Phot_All{phot_i}; % run model in constant dark
	Time_Phot = Time_Phot_All{phot_i};
	Phot_name = [int2str(options.photoperiod),'h'];
	
	Model_output_to_file=struct();
	for ig=1:nG
		% set the genotype
		options.genotype = mutant_genotypes{ig};
% 		figure('Name',string(join([options.genotype,{Phot_name}])))
		for clock_dynamics_model_i=1:length(Models)
			%% Clock model
			if clock_dynamics_model_i<3
				clock_species_colNames = {'LHYm','P','GIZ','GIE','LHYp','TOCm','PRR9p','PRR5m','PRR5p','GIp','TOCp','ZTL','EC','GIm','PRR9m','PRR7m','PRR7p','ELF4m','ELF4p','LHYpm','HY5','HFR1','ELF3m','ELF3pc','ELF3pn','COP1pnn','COP1pnd','LUXm','LUXp','COP1pc'};
				parameters = load_P2011_parameters(options.genotype,YHB);
				if clock_dynamics_model_i==1
					clock_dynamics = @P2011_dynamics_Red;
				elseif clock_dynamics_model_i==2
					clock_dynamics = @P2011_dynamics_COP1;
				end
				clock_dynamics_wrapper = @wrap_P2011_model_dynamics;
				y0=[1.0151 0.956 0.0755 0.0041 0.506 0.0977 0.0238 0.0731 0.0697 0.0196 0.0435 0.2505 0.0709 0.1017 0.0658 0.4016 0.1167 0.1012 0.207 0.0788 0.3102 0.0553 0.2991 0.1503 0.0286 0.65 0.2566 0.1012 0.576 0.3269];
				clock_Sidx = clock_Sidx1;
			else
				clock_species_colNames = {'LHYm','LHYp','CCA1m','CCA1p','P','PRR9m','PRR9p','PRR7m','PRR7p','PRR5m','PRR5c','PRR5n','TOC1m','TOC1n','TOC1c','ELF4m','ELF4p','ELF4d','ELF3m','ELF3p','ELF34','LUXm','LUXp','COP1c','COP1n','COP1d','ZTL','ZG','GIm','GIc','GIn','NOXm','NOXp','RVE8m','RVE8p'};
				parameters = load_F2014_parameters(options.genotype,paramSet,YHB);
				if clock_dynamics_model_i==3
					clock_dynamics = @F2014_dynamics_Red;
				elseif clock_dynamics_model_i==4
					clock_dynamics = @F2014_dynamics_COP1;
				end
				clock_dynamics_wrapper = @wrap_F2014_model_dynamics;
				y0=0.1*ones(1,35);
				clock_Sidx = clock_Sidx2;
			end
	
			[u,Tc,Yc]=entrain_and_run_clock_model(parameters,clock_dynamics,clock_dynamics_wrapper,c,options,y0);
	
% 			make_and_format_plots(Tc,Yc,clock_dynamics_model_i,options,clock_species_names,clock_nS,clock_Sidx,clock_ylimit_T)

			Model_output_to_file.(Models{clock_dynamics_model_i})(ig).Time = Tc;
			for ic=1:length(clock_species_colNames)
				Model_output_to_file.(Models{clock_dynamics_model_i})(ig).(clock_species_colNames{ic}) = Yc(:,ic);
			end
		end
	end
	
	for clock_dynamics_model_i=1:length(Models)
		if clock_dynamics_model_i<3
			clock_species_colNames = {'LHYm','P','GIZ','GIE','LHYp','TOCm','PRR9p','PRR5m','PRR5p','GIp','TOCp','ZTL','EC','GIm','PRR9m','PRR7m','PRR7p','ELF4m','ELF4p','LHYpm','HY5','HFR1','ELF3m','ELF3pc','ELF3pn','COP1pnn','COP1pnd','LUXm','LUXp','COP1pc'};
		else
			clock_species_colNames = {'LHYm','LHYp','CCA1m','CCA1p','P','PRR9m','PRR9p','PRR7m','PRR7p','PRR5m','PRR5c','PRR5n','TOC1m','TOC1n','TOC1c','ELF4m','ELF4p','ELF4d','ELF3m','ELF3p','ELF34','LUXm','LUXp','COP1c','COP1n','COP1d','ZTL','ZG','GIm','GIc','GIn','NOXm','NOXp','RVE8m','RVE8p'};
		end
		ModelHypFlMut=table();
	
		if options.photoperiod==0
			ModelHypFlMut.(Time_Phot) = Model_output_to_file.(Models{clock_dynamics_model_i})(1).Time+12;
		else
			ModelHypFlMut.(Time_Phot) = Model_output_to_file.(Models{clock_dynamics_model_i})(1).Time;
		end
		for ic=1:length(clock_species_colNames)
			for ig=1:nG
				% set the genotype
				options.genotype = mutant_genotypes{ig};
				ModelHypFlMut.(string(join(options.genotype,''))) = Model_output_to_file.(Models{clock_dynamics_model_i})(ig).(clock_species_colNames{ic});
			end
			writetable(ModelHypFlMut,['output\',int2str(floor(100*YHB/(YHB+1))),'\ModelClockMut_phot',Phot_name,'_',Models{clock_dynamics_model_i},'.xlsx'],'Sheet',clock_species_colNames{ic})
		end
	end
end

rmpath('PIF_CO_FT_model')
rmpath('circadian_module')
rmpath('parameters')
rmpath('phenology')

function [u,Tc,Yc]=entrain_and_run_clock_model(parameters,clock_dynamics,clock_dynamics_wrapper,c,options,y0)
	u=struct;
	
	% Initialise clock model (12 days)
	c.photoperiod = options.entrain;
	[~,Ycinit] = ode15s(@(t,y) clock_dynamics(t,y,parameters,c.dawn,c.dawn+c.photoperiod),0:0.1:(12*c.period),y0);
	y0 = Ycinit(end,:)';
	
	n_e=3;
	% Simulate for n_e days in entrainment conditions
	[Tce,Yce] = ode15s(@(t,y) clock_dynamics(t,y,parameters,c.dawn,c.dawn+c.photoperiod),0:0.1:(n_e*c.period),y0);
	y0 = Yce(end,:)';
	% Convert dynamics into input struct for PIF_CO_FT model
	u_e = clock_dynamics_wrapper(Tce,Yce,parameters);
		
	if options.entrain==options.photoperiod
		Tc = Tce;
		Yc = Yce;
		u = u_e;
	else
		c.photoperiod = options.photoperiod;
		% Simulate for five days in selected photoperiod conditions
		[Tcp,Ycp] = ode15s(@(t,y) clock_dynamics(t,y,parameters,c.dawn,c.dawn+c.photoperiod),0:0.1:(5*c.period),y0);
		% Convert dynamics into input struct for PIF_CO_FT model
		u_p = clock_dynamics_wrapper(Tcp,Ycp,parameters);

		clock_species_names = fieldnames(u_p);
		clock_nS = length(clock_species_names);
		for clock_j = 1:clock_nS
			if ismember(clock_species_names(clock_j),{'T'})
				u.(clock_species_names{clock_j}) = [u_e.(clock_species_names{clock_j})-24*n_e;u_p.(clock_species_names{clock_j})];
			else
				u.(clock_species_names{clock_j}) = [u_e.(clock_species_names{clock_j});u_p.(clock_species_names{clock_j})];
			end
		end
		Tc=[Tce-24*n_e;Tcp];
		Yc=[Yce;Ycp];
	end
end
	


function make_and_format_plots(Tc,Yc,clock_dynamics_model_i,options,clock_species_names,clock_nS,clock_Sidx,clock_ylimit_T)
	% Clock model - plot the results for each species on separate panels
	for clock_j = 1:clock_nS
		subplot(1,clock_nS,clock_j)
		plot(Tc,Yc(:,clock_Sidx(clock_j)),'DisplayName',int2str(clock_dynamics_model_i),'linewidth',2)
		hold on
		box on
		ylim([0,clock_ylimit_T(clock_j)])
		xlim([Tc(1),Tc(end)])
		xlabel('Time (ZT Hrs)')
		ylabel('Relative Expression')
		title(clock_species_names{clock_j})
		v = axis;

		if clock_dynamics_model_i==4
			if options.entrain==options.photoperiod
				for plot_i=0:ceil((Tc(end)-Tc(1))/24)
					a=fill([options.entrain,24,24,options.entrain]+(Tc(1)-mod(Tc(1),24)+plot_i*24),[0,0,v(4),v(4)],'k','FaceAlpha',0.3,'LineStyle','none');
					uistack(a,'bottom')
				end
			else
				for plot_i=0:(ceil((0-Tc(1))/24)-1)
					a=fill([options.entrain,24,24,options.entrain]+(Tc(1)-mod(Tc(1),24)+plot_i*24),[0,0,v(4),v(4)],'k','FaceAlpha',0.3,'LineStyle','none');
					uistack(a,'bottom')
				end
	
				if options.photoperiod==0
					for plot_i=0:ceil((Tc(end)-0)/24)
						a=fill([0,options.entrain,options.entrain,0]+(plot_i*24),[0,0,v(4),v(4)],'k','FaceAlpha',0.1,'LineStyle','none');
						uistack(a,'bottom')
						a=fill([options.entrain,24,24,options.entrain]+(plot_i*24),[0,0,v(4),v(4)],'k','FaceAlpha',0.3,'LineStyle','none');
						uistack(a,'bottom')
					end
				elseif options.photoperiod==24
					for plot_i=0:ceil((Tc(end)-0)/24)
						a=fill([options.entrain,24,24,options.entrain]+(plot_i*24),[0,0,v(4),v(4)],'k','FaceAlpha',0.1,'LineStyle','none');
						uistack(a,'bottom')
					end
				else
					for plot_i=0:ceil((Tc(end)-0)/24)
						a=fill([options.photoperiod,24,24,options.photoperiod]+(plot_i*24),[0,0,v(4),v(4)],'k','FaceAlpha',0.3,'LineStyle','none');
						uistack(a,'bottom')
					end
				end
			end
		end
	end
end
