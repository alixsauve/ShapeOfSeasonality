clc
clear all
close all

global a K r d c h
global epsilon theta

epsilon = 1; theta = 1;

 % parameter values used in Taylor et al. 2013, J. Math. Biol.
d = 2*pi; % predator's mortality
c = 1; % conversion efficiency
r = d; % intrinsic-growth rate of the prey
K = 1; % carrying capacity of the prey
h = 1/(4*pi); % handling time

rng(2019);

TFinal = 5000;
TInit = 0;
TStep = 1/52;
TSpan = [TInit:TStep:TFinal];
TRecord = [(TFinal - 4000):TStep:TFinal];

addpath('../Functions');

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-6, 'NonNegative', [1:2]);

aValues = [8 12 16];
ICRef = [0.3; 0.3];

ThetaCol = [59 154 178; 235 204 42; 242 26 0]/255;
a = aValues(3) * pi;

SubPlotCount = 1;
for epsilon = [0.3 0.36 0.38 0.41 0.44 0.47 0.48 0.5]
	ThetaCount = 1;
	for theta = [0 0.5 1]
		disp([epsilon theta]);
		subplot(2, 4, SubPlotCount); hold on;
		[T, Y] = odeIntegSeason(ICRef, TFinal, TRecord, @(t, y) seasonKill2Sp_Type2(t, y, false), options);
		Y = Y((T > (TFinal - 4000)), :); T = T(T > (TFinal - 4000));
		IndexBegY = find(ismember(T, unique(floor(T))));
		scatter(Y(IndexBegY, 2), Y(IndexBegY, 1), 10, ThetaCol(ThetaCount, :), '.'); xlim([0 1]); ylim([0 0.4]); hold on;
		ThetaCount = ThetaCount + 1;
	end
	xlim([0 1]); ylim([0 0.4]);
	text(min(xlim) + 0.05 * diff(xlim), max(ylim)*(1 - 0.1), char(SubPlotCount + 64));
	xlabel("Prey density", 'FontSize', 8); ylabel("Predator density", 'FontSize', 8);
  set(gca, 'FontSize', 8);
	SubPlotCount = SubPlotCount + 1;
end
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 5], 'PaperUnits', 'Inches', 'PaperSize', [10, 5]);
savefig('../Outputs/fig_stroboscmaps_t2.fig');
print('../Figures/fig_stroboscmaps_t2', '-dpng');

