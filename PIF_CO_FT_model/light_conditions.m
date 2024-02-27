function L = light_conditions(t,c)
%
%   Copyright 
	
	if isfield(c,'LightAmp') %%%%%%%%% Not currently passed because c is not original c!
		LightAmp = c.LightAmp; %The light intensity. %%%%%%%%% Not currently passed because c is not original c!
	else
		LightAmp = 1; %The light intensity.
	end
	twilightPer=0.00005; %The duration of time between value of force in dark and value of force in light.
	
	if length(c.photoperiod)==1
		if (c.photoperiod == 0)
			L = 0*t;
		elseif (c.photoperiod == c.period)
    		L = LightAmp+0*t;
		else
	    	L = interp1([c.dawn, c.dawn+twilightPer, c.dawn+c.photoperiod-twilightPer, c.dawn+c.photoperiod, c.period],[0, LightAmp, LightAmp, 0, 0],mod(t,c.period));
		end
	else
    	L = interp1(c.photoperiod(:,1),c.photoperiod(:,2),mod(t,c.period));
	end

end
