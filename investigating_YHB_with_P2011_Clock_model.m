% This script simulates and plots the clock model outputs for the P2011
% model and the F2014 model in the given photoperiod conditions. It can
% also run (transition to) constant dark or constant light due to a change
% in "light_conditions"


% Include model folders and plotting folders in path
if contains(pwd,'nenya')
	addpath('C:\Users\nenya\OneDrive - University of Glasgow\Projects\Matt\Seaton 2015\published_model\plotting_tools')
else
	addpath('C:\Users\ra134k\OneDrive - University of Glasgow\Projects\Matt\Seaton 2015\published_model\plotting_tools')
end
addpath('PIF_CO_FT_model')

%% Clock model
clock_dynamics_model_i=1;
clock_species_names_all = {'LHYm','P','GIZ','GIE','LHYp','TOCm','PRR9p','PRR5m','PRR5p','GIp','TOCp','ZTL','EC','GIm','PRR9m','PRR7m','PRR7p','ELF4m','ELF4p','LHYpm','HY5','HFR1','ELF3m','ELF3pc','ELF3pn','COP1pnn','COP1pnd','LUXm','LUXp','COP1pc'};
clock_dynamics = @P2011_dynamics_COP1;
clock_dynamics_wrapper = @wrap_P2011_model_dynamics;
y0=[1.0151 0.956 0.0755 0.0041 0.506 0.0977 0.0238 0.0731 0.0697 0.0196 0.0435 0.2505 0.0709 0.1017 0.0658 0.4016 0.1167 0.1012 0.207 0.0788 0.3102 0.0553 0.2991 0.1503 0.0286 0.65 0.2566 0.1012 0.576 0.3269];
clock_Sidx = 1:length(clock_species_names_all);% [1, 3, 5, 6, 9, 11, 29, 25]; % {'LHYm'}    {'CCA1m'}    {'GIm'}    {'PRR9m'}    {'COP1pnn'}
clock_nS = length(clock_Sidx);
clock_species_names = [clock_species_names_all(clock_Sidx)];

clock_ylimit_T=[1,1,1,1,1];
clock_species_names_u = {'LHY','EC','PRR9','PRR7','PRR5','TOC1','cP','COP1n_n','GIn'};
%clock_ylimit_T_u_gens={[1.6,0.14,1.1,4,1.5,0.25,1,1,4.5],[1.6,0.14,1.1,4,1.5,0.86,1,1,8.5]};
clock_ylimit_T_u_gens={[1.6,0.14,1.1,4,1.5,0.25,1,1,4.5],[1.6,0.14,1.1,4,1.5,0.86,1,1,8.5],[1.6,0.14,1.1,4,1.5,0.86,1,1,8.5],[1.6,0.14,1.1,4,1.5,0.86,1,1,8.5]};

mutant_genotypes={{'wt'},{'elf3'},{'YHB'},{'elf3','YHB'}};
nG=length(mutant_genotypes);

Models = {'P2011', 'F2014'};

% set the parameter set for F2014 model
paramSet=1;

% set the common conditions in this case
options = struct();
options.temperature = 22; % temperature (oC) (only 22 or 27oC accepted)
options.entrain = 12; % entrain model at 12/12
options.photoperiod = 0;

% Load light conditions into 'c' for common light function
c.period = 24;
c.phase = 0;
c.dawn = 0;

figure
set(gcf,'windowstyle', 'docked')
for ig=1:nG
	% set the genotype
	options.genotype = mutant_genotypes{ig};
	clock_ylimit_T_u=clock_ylimit_T_u_gens{ig};
% 	figure('Name',string(join(options.genotype)))
	parameters = load_P2011_parameters(options.genotype);

	[u,Tc,Yc,Lc]=entrain_and_run_clock_model(parameters,clock_dynamics,clock_dynamics_wrapper,c,options,y0);

	make_and_format_plots(Tc,Yc,Lc,clock_dynamics_model_i,options,clock_species_names,clock_nS,clock_Sidx,clock_ylimit_T)
	legend('show')
% 		figure('Name',[char(string(join(options.genotype))),' u model: ',int2str(clock_dynamics_model_i)])
% 		make_and_format_plots_from_u(u,Lc,Models{clock_dynamics_model_i},options,clock_species_names_u,clock_ylimit_T_u)
	
end

if contains(pwd,'nenya')
	rmpath('C:\Users\nenya\OneDrive - University of Glasgow\Projects\Matt\Seaton 2015\published_model\plotting_tools')
else
	rmpath('C:\Users\ra134k\OneDrive - University of Glasgow\Projects\Matt\Seaton 2015\published_model\plotting_tools')
end
rmpath('PIF_CO_FT_model')

function [u,Tc,Yc,Lc]=entrain_and_run_clock_model(parameters,clock_dynamics,clock_dynamics_wrapper,c,options,y0)
	u=struct;
	
	% Initialise clock model (12 days)
	c.photoperiod = options.entrain;
	[~,Ycinit] = ode15s(@(t,y) clock_dynamics(t,y,parameters,c.dawn,c.dawn+c.photoperiod),0:0.1:(12*c.period),y0);
	y0 = Ycinit(end,:)';
	
	% Simulate for one day in entrainment conditions
	[Tce,Yce] = ode15s(@(t,y) clock_dynamics(t,y,parameters,c.dawn,c.dawn+c.photoperiod),0:0.1:c.period,y0);
	y0 = Yce(end,:)';
	% Convert dynamics into input struct for PIF_CO_FT model
	u_e = clock_dynamics_wrapper(Tce,Yce,parameters);
	% Store light conditions
	Lce = light_conditions(Tce,c);
	
	if options.entrain==options.photoperiod
		Tc = Tce;
		Yc = Yce;
		Lc = Lce;
		u = u_e;
	else
		c.photoperiod = options.photoperiod;
		% Simulate for three days
		[Tcp,Ycp] = ode15s(@(t,y) clock_dynamics(t,y,parameters,c.dawn,c.dawn+c.photoperiod),0:0.1:(3*c.period),y0);
		% Convert dynamics into input struct for PIF_CO_FT model
		u_p = clock_dynamics_wrapper(Tcp,Ycp,parameters);
		% Store light conditions
		Lcp = light_conditions(Tcp,c);

		clock_species_names = fieldnames(u_p);
		clock_nS = length(clock_species_names);
		for clock_j = 1:clock_nS
			if ismember(clock_species_names(clock_j),{'T'})
				u.(clock_species_names{clock_j}) = [u_e.(clock_species_names{clock_j})-24;u_p.(clock_species_names{clock_j})];
			else
				u.(clock_species_names{clock_j}) = [u_e.(clock_species_names{clock_j});u_p.(clock_species_names{clock_j})];
			end
		end
		Tc=[Tce-24;Tcp];
		Yc=[Yce;Ycp];
		Lc=[Lce;Lcp];
	end
end
	

function make_and_format_plots_from_u(u,Lc,clock_dynamics_model,options,clock_species_names_u,clock_ylimit_T)
	clock_species_names = fieldnames(u);
	clock_nS = length(clock_species_names);

	beforeT=1;
	% Clock model - plot the results for each species on separate panels
	for clock_j = 1:clock_nS
		if ismember(clock_species_names(clock_j),{'T'})
			beforeT=0;
			continue
		end
		subplot(1,clock_nS,clock_j+beforeT)
		plot(u.T,u.(clock_species_names{clock_j}),'DisplayName',clock_dynamics_model)
		hold on
		box on
		ylim([0,clock_ylimit_T(ismember(clock_species_names_u,clock_species_names(clock_j)))])
		xlabel('Time (ZT Hrs)')
		ylabel('Relative Expression')

		circaplotlocal(options)

		v = axis;
		text(v(1)+0.5,v(4)*0.8,clock_species_names{clock_j},'FontAngle','italic')
	end
	
	% Plot the light conditions
	subplot(1,clock_nS+1,1)
	plot(u.T,Lc,'b--')
	hold on
	box on
	ylim([0,1.1])
	xlabel('Time (ZT Hrs)')
	ylabel('Light')

	circaplotlocal(options)
end

function make_and_format_plots(Tc,Yc,Lc,clock_dynamics_model_i,options,clock_species_names,clock_nS,clock_Sidx,clock_ylimit_T)
% 	figure('Name',[char(string(join(options.genotype))),' model: ',int2str(clock_dynamics_model_i)])
	
	% Clock model - plot the results for each species on separate panels
	for clock_j = 1:clock_nS
		if clock_nS>6
			subplot(4,ceil(clock_nS/4),clock_j)
		else
			subplot(1,clock_nS,clock_j)
		end
		plot(Tc,Yc(:,clock_Sidx(clock_j)),'DisplayName',string(join(options.genotype)))
		hold on
		box on
% 		ylim([0,clock_ylimit_T(clock_j)])
		xlabel('Time (ZT Hrs)')
		ylabel('Relative Expression')
		title(clock_species_names{clock_j})
		circaplotlocal(options)

% 		v = axis;
% 		text(v(1)+0.5,v(4)*0.8,clock_species_names{clock_j},'FontAngle','italic')
	end
end

function circaplotlocal(options)
	if options.entrain==options.photoperiod
		circaplot([],[],[0,options.entrain],{[1 1 1],[0 0 0]},24)
		xlim([0,24])
	else
		circaplot([],[],[0,options.entrain]-24,['w','k'],24)

		if options.photoperiod==0
			circaplot([],[],[0,options.entrain],{[0.5 0.5 0.5],[0 0 0]},24)
			circaplot([],[],[0,options.entrain]+24,{[0.5 0.5 0.5],[0 0 0]},48)
			circaplot([],[],[0,options.entrain]+48,{[0.5 0.5 0.5],[0 0 0]},72)
			xlabel('Time in constant darkness (hours)')
			set(gca,'xtick',[-24   -12     0    12    24    36    48    60    72],'xticklabel',{'-12','0','12','24','36','48','60','72','84'})
		elseif options.photoperiod==24
			circaplot([],[],[0,options.entrain],{[1 1 1],[0.5 0.5 0.5]},24)
			circaplot([],[],[0,options.entrain]+24,{[1 1 1],[0.5 0.5 0.5]},48)
			circaplot([],[],[0,options.entrain]+48,{[1 1 1],[0.5 0.5 0.5]},72)
		else
			circaplot([],[],[0,options.photoperiod],{[1 1 1],[0 0 0]},24)
			circaplot([],[],[0,options.photoperiod]+24,{[1 1 1],[0 0 0]},48)
			circaplot([],[],[0,options.photoperiod]+48,{[1 1 1],[0 0 0]},72)
		end
		
		xlim([-24,3*24])
	end
end