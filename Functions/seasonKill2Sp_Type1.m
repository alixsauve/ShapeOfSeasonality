function xdot = seasonKill2Sp_Type1(tt, xt, scale)
	% seasonKill2Sp_Type1() is the system of differential equation for the endogeneous model: a generalised Lotka-Volterra model for 2 species (type I functional responses)
	% only attack rate is seasonal

	% xt is a vector of species abundance at time t
	% tt is a scalar for time
	% scale is a logical input
	if islogical(scale) == 0
		error('Error in seasonKill2Sp(): `scale` must be a logical.');
	end


	global a K r d c
	global epsilon theta

	CorrSignal = NaN;
	if scale == true
		T = [1/52:1/52:1]; % a time vector over one period
		SignalT_theta = 0;
		SignalT_1 = arrayfun(@(X) sin(2 * pi * X), T);
		for i1 = T
			SignalT_theta = [SignalT_theta sin(2 * pi * i1) / abs(sin(2 * pi * i1)) * abs(sin(2 * pi * i1))^theta];
		end
		CorrSignal = std(SignalT_1) / std(SignalT_theta);
	end

	if epsilon ~= 0
		if sin(2 * pi * tt) == 0
			fAR = 1; % avoid division by 0
		else
			SignalAR = sin(2 * pi * tt) / abs(sin(2 * pi * tt)) * abs(sin(2 * pi * tt))^theta;
			if scale
				SignalAR = SignalAR * CorrSignal ;
			end
			fAR = (1 + epsilon * SignalAR); % forcing on the attack rate
		end
	else
		fAR = 1;
	end

	x1dot = xt(1) * (c * a * fAR * xt(2) - d); % predator's growth rate
	x2dot = xt(2) * (r - r * xt(2) / K - a * fAR * xt(1)); %prey's growth rate
  	xdot = [x1dot; x2dot];
end
