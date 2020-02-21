function [TSol, YSol] = odeIntegSeason(Y0, TFinal, TSave, FuncGrowth, Options)
	% odeIntegSeason() runs an ode integrator which is stopped and reinitialised each time we change season.
	% It returns TSol and YSol, respectively the time steps of the numerical solution of the ODE, and the numerical solution.

	% Y0 is a vector for the initial conditions
	% TFinal is a numeric for the final time of integration
	% FuncGrowth is a function handle for the ode
	if ~isa(FuncGrowth, 'function_handle')
		error('FuncGrowth must be a function handle for the ode.');
	end
	
	TStep = 1/52; % sample every month
	TInit = 0; % Initial time
	YSol = [];
	TSol = [];
	Y0i = Y0;
	while TInit < TFinal
		[ti, Yi] = ode15s(FuncGrowth, [TInit:TStep:(TInit+0.5)], Y0i, Options);
		if (TInit + 0.5) >= TSave(1)
			YSol = [YSol; Yi]; % paste together the solutions to get the whole time serie
			TSol = [TSol; ti];
		end
		Y0i = Yi(end, :); % save the last state
		TInit = ti(end);
	end

end
