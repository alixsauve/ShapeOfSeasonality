clc
clear all
close all

global a K r d c h
global epsilon theta

epsilon = 0; theta = 1;

 % parameter values used in Taylor et al. 2013, J. Math. Biol.
d = 2*pi; % predator's mortality
c = 1; % conversion efficiency
r = d; % intrinsic-growth rate of the prey
K = 1; % carrying capacity of the prey
h = 1/(4*pi); % handling time

rng(2019);

TFinal = 2000;
TInit = 0;
TStep = 1/52;
TSpan = [TInit:TStep:TFinal];
TRecord = [(TFinal - 50):TStep:TFinal];

TimeVect = 0:1/52:2;
TimeVectPeriod = 0:1/52:1;

addpath('../Functions');

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-6, 'NonNegative', [1:2]);

ThetaCol = [59 154 178; 83 165 185; 107 177 193; 143 187 165; 189 195 103; 235 204 42; 230 192 25; 227 180 8; 228 145 0; 235 85 0; 242 26 0]/255;

% the joint effect of theta and epsilon
% 3d plots (one for each value of a)
IC = [0.3; 0.3];
aValues = [8 12 16];
for aIndex = 1:3
	a = aValues(aIndex)*pi;
	for epsilon = 0.1:0.1:1
		thetaCount = 1;
		for theta = 0:0.1:1
			disp([aIndex epsilon theta]);

			% first, we need to verify that |epsilon * ĝ(t)| <= 1, where g(t) is the forcing signal
			ForcingSignal_epsilon_theta = arrayfun(@(X) forcingSignal(X, epsilon, theta), TimeVect);
			SdPeriod = std(arrayfun(@(X) forcingSignal(X, epsilon, theta), TimeVectPeriod));
			% type 1
			disp('Simulating predator-prey community with type I functional response...');
			[T, Y] = odeIntegSeason(IC, TFinal, TRecord, @(t, y) seasonKill2Sp_Type1(t, y, true), options);
			Index10Y = (T >= (TFinal - 10)); % index of the last 10 years
			EValues = epsilon * ones(sum(Index10Y), 1);
			Prey10Y = Y(Index10Y, 2); Pred10Y = Y(Index10Y, 1);
	
			subplot(3, 2, 2*aIndex - 1); hold on;
			if (epsilon == 0.1) & (theta == 0)
				set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 10], 'PaperUnits', 'Inches', 'PaperSize', [5, 10])
				set(gca, 'View', [-25 10]);
				xlim([0 1]); ylim([0 1]); zlim([0 0.8]);
				if aIndex == 1
					xlim([0 1]); ylim([0 1]); zlim([0 0.5]);
				end
				text(min(xlim) + 0.05*diff(xlim), max(ylim) - 0.05*diff(ylim), max(zlim) - 0.05*diff(zlim), char(2*aIndex - 1 + 64));
				xlabel('\epsilon'); ylabel('Prey density'); zlabel('Predator density');
			  	title(strcat('a = ',  int2str(aValues(aIndex)), '\pi'), 'FontSize', 8);
				set(gca, 'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', 'FontSize', 8);
			end
			plot3(EValues, Prey10Y, Pred10Y, 'Color', ThetaCol(thetaCount, :)); hold on;

			% type 2
			disp('Simulating predator-prey community with type II functional response...');
			[T, Y] = odeIntegSeason(IC, TFinal, TRecord, @(t, y) seasonKill2Sp_Type2(t, y, true), options);
			Index10Y = (T >= (TFinal - 10)); % index of the last 10 years
			Prey10Y = Y(Index10Y, 2); Pred10Y = Y(Index10Y, 1);
	
			subplot(3, 2, 2*aIndex); hold on;
			if (epsilon == 0.1) & (theta == 0)
				set(gca, 'View', [-25 10]);
				xlim([0 1]); ylim([0 1]); zlim([0 0.8]);
				if aIndex == 1
					xlim([0 1]); ylim([0 1]); zlim([0 0.5]);
				end
				text(min(xlim) + 0.05*diff(xlim), max(ylim) - 0.05*diff(ylim), max(zlim) - 0.05*diff(zlim), char(2*aIndex + 64));
				xlabel('\epsilon'); ylabel('Prey density'); zlabel('Predator density');
			  	title(strcat('a = ',  int2str(aValues(aIndex)), '\pi'), 'FontSize', 8);
				set(gca, 'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', 'FontSize', 8);
			end
			plot3(EValues, Prey10Y, Pred10Y, 'Color', ThetaCol(thetaCount, :)); hold on;
			thetaCount = thetaCount + 1;
		end
	end
end
subplot(3, 2, 5);
legend( {'\theta = 0', '\theta = 0.1', '\theta = 0.2', '\theta = 0.3', '\theta = 0.4', '\theta = 0.5', '\theta = 0.6', '\theta = 0.7', '\theta = 0.8', '\theta = 0.9', '\theta = 1'}, 'Location', 'West');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9, 10], 'PaperUnits', 'Inches', 'PaperSize', [9, 10]);
savefig('../Outputs/fig_phaseportraits_controlvar.fig');
print('../Figures/fig_phaseportraits_controlvar', '-dpng');
