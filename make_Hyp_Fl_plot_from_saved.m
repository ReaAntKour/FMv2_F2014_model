figure
for ig=1:nG
	for clock_dynamics_model_i=1:2
		subplot(1,3,1)
		plot(Phot,Hyp_length(ig,:,clock_dynamics_model_i),'DisplayName',string(join([mutant_genotypes{ig},Models{clock_dynamics_model_i}])))
    	hold on
		box on
		xlim([0,20])
		ylim([0 7])
		xlabel('Photoperiod (Hrs)')
		ylabel('Hypocotyl length')
	
		subplot(1,3,2)
		plot(Phot,FT_area(ig,:,clock_dynamics_model_i),'DisplayName',string(join([mutant_genotypes{ig},Models{clock_dynamics_model_i}])))
		hold on
		box on
		xlim([0,20])
		ylim([0 15])
		xlabel('Photoperiod (Hrs)')
		ylabel('Daily FT area')
	
		subplot(1,3,3)
		plot(Phot,Days_to_flower(ig,:,clock_dynamics_model_i),'DisplayName',string(join([mutant_genotypes{ig},Models{clock_dynamics_model_i}])))
		hold on
		box on
		xlim([0,20])
		ylim([0 100])
		xlabel('Photoperiod (Hrs)')
		ylabel('Days to Flower')
	end
end