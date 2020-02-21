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

% type 1, epsilon = 0.8, a = 12 * pi
subplot(2, 2, 1); hold on;
a = aValues(2) * pi; epsilon = 0.8;
IC = [0.3523 0.3522; 0.7190 0.4315];
YPrey = []; YPred = []; PEst = [];
for ICIndex = 1:size(IC, 1)
	[T, Y] = ode15s(@(t, y) seasonKill2Sp_Type1(t, y, false), [0 TRecord], IC(ICIndex, :), options);
	% estimating the period
	Y = Y((T > (TFinal - 100)), :); T = T(T > (TFinal - 100));
	PEst = [PEst, max([periodAnalysis(T, Y(:, 1)) periodAnalysis(T, Y(:, 2))])];
	if isempty(periodAnalysis(T, Y))
		PEst = [PEst nan];
	end
	% saving time series to display
	Y = Y((T > (TFinal - 10)), :); T = T(T > (TFinal - 10));
	YPrey = [YPrey, Y(:, 2)]; YPred = [YPred, Y(:, 1)];
end
% plot time series
CeilPred = ceil(max(max(YPred))*10)/10;
CeilPreyPred = CeilPred + ceil(max(max(YPrey))*10)/10 + 0.1;
line([min(T) max(T)], [CeilPred + 0.1 CeilPred + 0.1], 'Color', 'k'); hold on;
for ICIndex = 1:size(IC, 1)
	Col = 'k';
	if ~isnan(PEst(ICIndex))
		Col = PeriodCol(PEst(ICIndex), :);
	end
	plot(T, YPred(:, ICIndex), '-.', 'Color', Col); hold on;
	plot(T, YPrey(:, ICIndex) + CeilPred + 0.1, 'Color', Col); hold on;
end
ylim([0 CeilPreyPred]);
set(gca, 'YTick', [0:0.2:CeilPreyPred]);
set(gca, 'YTickLabel', [0:0.2:CeilPred, 0:0.2:CeilPreyPred]);
text(min(xlim) + 0.05 * diff(xlim), max(ylim) - 0.1 * diff(ylim), char(1 + 64));
xlabel('Time', 'FontSize', 8); ylabel('Species density', 'FontSize', 8);
set(gca, 'FontSize', 8);

% type 2, epsilon = 0.25, a = 12 * pi
subplot(2, 2, 2); hold on;
a = aValues(2) * pi; epsilon = 0.25;
IC = [0.1604 0.4053; 0.3985 0.7742];
YPrey = []; YPred = []; PEst = [];
for ICIndex = 1:size(IC, 1)
	[T, Y] = ode15s(@(t, y) seasonKill2Sp_Type2(t, y, false), [0 TRecord], IC(ICIndex, :), options);
	% estimating the period
	Y = Y((T > (TFinal - 100)), :); T = T(T > (TFinal - 100));
	PEst = [PEst, max([periodAnalysis(T, Y(:, 1)) periodAnalysis(T, Y(:, 2))])];
	if isempty(periodAnalysis(T, Y))
		PEst = [PEst nan];
	end
	% saving time series to display
	Y = Y((T > (TFinal - 10)), :); T = T(T > (TFinal - 10));
	YPrey = [YPrey, Y(:, 2)]; YPred = [YPred, Y(:, 1)];
end
% plot time series
CeilPred = ceil(max(max(YPred))*10)/10;
CeilPreyPred = CeilPred + ceil(max(max(YPrey))*10)/10 + 0.1;
line([min(T) max(T)], [CeilPred + 0.1 CeilPred + 0.1], 'Color', 'k'); hold on;
for ICIndex = 1:size(IC, 1)
	Col = 'k';
	if ~isnan(PEst(ICIndex))
		Col = PeriodCol(PEst(ICIndex), :);
	end
	plot(T, YPred(:, ICIndex), '-.', 'Color', Col); hold on;
	plot(T, YPrey(:, ICIndex) + CeilPred + 0.1, 'Color', Col); hold on;
end
ylim([0 CeilPreyPred]);
set(gca, 'YTick', [0:0.1:CeilPreyPred]);
set(gca, 'YTickLabel', [0:0.2:CeilPred, 0:0.2:CeilPreyPred]);
text(min(xlim) + 0.05 * diff(xlim), max(ylim) - 0.1 * diff(ylim), char(2 + 64));
xlabel('Time', 'FontSize', 8); ylabel('Species density', 'FontSize', 8);
set(gca, 'FontSize', 8);

% type 1, a = 16 * pi, epsilon = 0.82
subplot(2, 2, 3); hold on;
a = aValues(3) * pi; epsilon = 0.82;
IC = [0.8452 0.7915; 0.3233 0.2756];
YPrey = []; YPred = []; PEst = [];
for ICIndex = 1:size(IC, 1)
	[T, Y] = ode15s(@(t, y) seasonKill2Sp_Type1(t, y, false), [0 TRecord], IC(ICIndex, :), options);
	% estimating the period
	Y = Y((T > (TFinal - 100)), :); T = T(T > (TFinal - 100));
	PEst = [PEst, max([periodAnalysis(T, Y(:, 1)) periodAnalysis(T, Y(:, 2))])];
	if isempty(periodAnalysis(T, Y))
		PEst = [PEst nan];
	end
	% saving time series to display
	Y = Y((T > (TFinal - 10)), :); T = T(T > (TFinal - 10));
	YPrey = [YPrey, Y(:, 2)]; YPred = [YPred, Y(:, 1)];
end
% plot time series
CeilPred = ceil(max(max(YPred))*10)/10;
CeilPreyPred = CeilPred + ceil(max(max(YPrey))*10)/10 + 0.1;
line([min(T) max(T)], [CeilPred + 0.1 CeilPred + 0.1], 'Color', 'k'); hold on;
for ICIndex = 1:size(IC, 1)
	Col = 'k';
	if ~isnan(PEst(ICIndex))
		Col = PeriodCol(PEst(ICIndex), :);
	end
	plot(T, YPred(:, ICIndex), '-.', 'Color', Col); hold on;
	semilogy(T, YPrey(:, ICIndex) + CeilPred + 0.1, 'Color', Col); hold on;
end
ylim([0 CeilPreyPred]);
set(gca, 'YTick', [0:0.2:CeilPreyPred]);
set(gca, 'YTickLabel', [0:0.2:CeilPred, 0:0.2:CeilPreyPred]);
text(min(xlim) + 0.05 * diff(xlim), max(ylim) - 0.1 * diff(ylim), char(3 + 64));
xlabel('Time', 'FontSize', 8); ylabel('Species density', 'FontSize', 8);
set(gca, 'FontSize', 8);

% type 2, a = 16 * pi, epsilon = 0.59
subplot(2, 2, 4); hold on;
a = aValues(3) * pi; epsilon = 0.59;
IC = [0.7554 0.8586; 0.6036 0.9601];
YPrey = []; YPred = []; PEst = [];
for ICIndex = 1:size(IC, 1)
	[T, Y] = ode15s(@(t, y) seasonKill2Sp_Type2(t, y, false), [0 TRecord], IC(ICIndex, :), options);
	% estimating the period
	Y = Y((T > (TFinal - 100)), :); T = T(T > (TFinal - 100));
	PEst = [PEst, max([periodAnalysis(T, Y(:, 1)) periodAnalysis(T, Y(:, 2))])];
	if isempty(periodAnalysis(T, Y))
		PEst = [PEst nan];
	end
	% saving time series to display
	Y = Y((T > (TFinal - 10)), :); T = T(T > (TFinal - 10));
	YPrey = [YPrey, Y(:, 2)]; YPred = [YPred, Y(:, 1)];
end
% plot time series
CeilPred = ceil(max(max(YPred))*10)/10;
CeilPreyPred = CeilPred + ceil(max(max(YPrey))*10)/10 + 0.1;
line([min(T) max(T)], [CeilPred + 0.1 CeilPred + 0.1], 'Color', 'k'); hold on;
for ICIndex = 1:size(IC, 1)
	Col = 'k';
	if ~isnan(PEst(ICIndex))
		Col = PeriodCol(PEst(ICIndex), :);
	end
	plot(T, YPred(:, ICIndex), '-.', 'Color', Col); hold on;
	plot(T, YPrey(:, ICIndex) + CeilPred + 0.1, 'Color', Col); hold on;
end
ylim([0 CeilPreyPred]);
set(gca, 'YTick', [0:0.2:CeilPreyPred]);
set(gca, 'YTickLabel', [0:0.2:CeilPred, 0:0.2:CeilPreyPred]);
text(min(xlim) + 0.05 * diff(xlim), max(ylim) - 0.1 * diff(ylim), char(4 + 64));
xlabel('Time', 'FontSize', 8); ylabel('Species density', 'FontSize', 8);
set(gca, 'FontSize', 8);

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 6, 6], 'PaperUnits', 'Inches', 'PaperSize', [6, 86]);
savefig('../Outputs/fig_time_series.fig');
print('../Figures/fig_time_series', '-dpng');
