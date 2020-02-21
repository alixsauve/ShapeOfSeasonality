% In this script, I draw attractors and their bassins of attraction when we observe coexistence in the bifurcation diagrams.

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
aValues = [8 12 16];

rng(2019);

TFinal = 2000;
TInit = 0;
TStep = 1/52;
TSpan = [TInit:TStep:TFinal];
TRecord = [(TFinal - 100):TStep:TFinal];

addpath('../Functions');

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-6, 'NonNegative', [1:2]);

PeriodCol = [228 26 28; 55 126 184; 77 175 74; 152 78 163; 255 127 0; 255 255 51; 166 86 40; 247 129 191; 153 153 153]/255;

% type I, a = 12 * pi, epsilon = 0.8
a = aValues(2) * pi; epsilon = 0.8;
subplot(2, 2, 1); hold on;
for Rep = 1:2000
	IC = rand(2, 1); % pick initial conditions between 0.1 and 1
	[T, Y] = ode15s(@(t, y) seasonKill2Sp_Type1(t, y, false), [0 TRecord], IC, options);
	Y = Y((T > (TFinal - 100)), :); T = T(T > (TFinal - 100));
	
	% pick colour depending on periodicity
	Col = 'k';
	PEst = periodAnalysis(T, Y);
	IndexBegY = find(ismember(T, unique(floor(T) + 0.25)));
	if ~isnan(PEst)
		Col = PeriodCol(PEst, :);
		scatter(Y(IndexBegY(end-PEst-1:end), 2), Y(IndexBegY(end-PEst-1:end), 1), 20, Col, 'filled'); hold on; % plot PEst points to reduce the size of the figure
	else
		scatter(Y(IndexBegY, 2), Y(IndexBegY, 1), 20, Col, 'filled'); hold on;
	end
	scatter(IC(2), IC(1), 5, Col, 'filled'); hold on;
	disp([1 2 Rep PEst]);
end
xlim([0 1]); ylim([0 1]);
text(min(xlim) + 0.05 * diff(xlim), max(ylim) - 0.1 * diff(ylim), char(1 + 64));
xlabel('Prey density', 'FontSize', 8); ylabel('Predator density', 'FontSize', 8);
set(gca, 'FontSize', 8);

% type II, a = 12 * pi, epsilon = 0.25
a = aValues(2) * pi; epsilon = 0.25;
subplot(2, 2, 2); hold on;
for Rep = 1:2000
	IC = rand(2, 1); % pick initial conditions between 0.1 and 1
	[T, Y] = ode15s(@(t, y) seasonKill2Sp_Type2(t, y, false), [0 TRecord], IC, options);
	Y = Y((T > (TFinal - 100)), :); T = T(T > (TFinal - 100));
	
	% pick colour depending on periodicity
	Col = 'k';
	PEst = periodAnalysis(T, Y);
	IndexBegY = find(ismember(T, unique(floor(T) + 0.25)));
	if ~isnan(PEst)
		Col = PeriodCol(PEst, :);
		scatter(Y(IndexBegY(end-PEst-1:end), 2), Y(IndexBegY(end-PEst-1:end), 1), 20, Col, 'filled'); hold on; % plot PEst points to reduce the size of the figure
	else
		scatter(Y(IndexBegY, 2), Y(IndexBegY, 1), 20, Col, 'filled'); hold on;
	end
	scatter(IC(2), IC(1), 5, Col, 'filled'); hold on;
	disp([2 2 Rep PEst]);
end
xlim([0 1]); ylim([0 1]);
text(min(xlim) + 0.05 * diff(xlim), max(ylim) - 0.1 * diff(ylim), char(2 + 64));
xlabel('Prey density', 'FontSize', 8); ylabel('Predator density', 'FontSize', 8);
set(gca, 'FontSize', 8);

% type I, a = 16 * pi, epsilon = 0.86 / was 0.82
a = aValues(3) * pi; epsilon = 0.86;
subplot(2, 2, 3); hold on;
for Rep = 1:2000
	IC = rand(2, 1); % pick initial conditions between 0.1 and 1
	[T, Y] = ode15s(@(t, y) seasonKill2Sp_Type1(t, y, false), [0 TRecord], IC, options);
	Y = Y((T > (TFinal - 100)), :); T = T(T > (TFinal - 100));
	
	% pick colour depending on periodicity
	Col = 'k';
	PEst = periodAnalysis(T, Y);
	IndexBegY = find(ismember(T, unique(floor(T) + 0.25)));
	if ~isnan(PEst)
		Col = PeriodCol(PEst, :);
		scatter(Y(IndexBegY(end-PEst-1:end), 2), Y(IndexBegY(end-PEst-1:end), 1), 20, Col, 'filled'); hold on; % plot PEst points to reduce the size of the figure
	else
		scatter(Y(IndexBegY, 2), Y(IndexBegY, 1), 20, Col, 'filled'); hold on;
	end
	scatter(IC(2), IC(1), 5, Col, 'filled'); hold on;
	disp([1 3 Rep PEst]);
end
xlim([0 1]); ylim([0 1]);
text(min(xlim) + 0.05 * diff(xlim), max(ylim) - 0.1 * diff(ylim), char(3 + 64));
xlabel('Prey density', 'FontSize', 8); ylabel('Predator density', 'FontSize', 8);
set(gca, 'FontSize', 8);

% type II, a = 16 * pi, epsilon = 0.59 / 0.68: 6T + chaos
a = aValues(3) * pi; epsilon = 0.59;
subplot(2, 2, 4); hold on;
for Rep = 1:2000
	IC = rand(2, 1); % pick initial conditions between 0.1 and 1
	[T, Y] = ode15s(@(t, y) seasonKill2Sp_Type2(t, y, false), [0 TRecord], IC, options);
	Y = Y((T > (TFinal - 100)), :); T = T(T > (TFinal - 100));
	
	% pick colour depending on periodicity
	Col = 'k';
	PEst = periodAnalysis(T, Y);
	IndexBegY = find(ismember(T, unique(floor(T) + 0.25)));
	if ~isnan(PEst)
		Col = PeriodCol(PEst, :);
		scatter(Y(IndexBegY(end-PEst-1:end), 2), Y(IndexBegY(end-PEst-1:end), 1), 20, Col, 'filled'); hold on; % plot PEst points to reduce the size of the figure
	else
		scatter(Y(IndexBegY, 2), Y(IndexBegY, 1), 20, Col, 'filled'); hold on;
	end
	scatter(IC(2), IC(1), 5, Col, 'filled'); hold on;
	disp([2 3 Rep PEst]);
end
xlim([0 1]); ylim([0 1]);
text(min(xlim) + 0.05 * diff(xlim), max(ylim) - 0.1 * diff(ylim), char(4 + 64));
xlabel('Prey density', 'FontSize', 8); ylabel('Predator density', 'FontSize', 8);
set(gca, 'FontSize', 8);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8, 8], 'PaperUnits', 'Inches', 'PaperSize', [8, 8]);
savefig('../Outputs/fig_basins_attraction.fig');
print('../Figures/fig_basins_attraction', '-dpng');

