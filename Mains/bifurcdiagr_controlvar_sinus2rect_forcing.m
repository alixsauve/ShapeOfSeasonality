% Because the simulations are long to run, we recommend to segment the following script into smaller subsets so as to generate each subplot separately.
% Figures can be gathered within the same multi-panel figure afterwards with the help of 'SubPlot':
% Farhad Sedaghati (2020). SubPlot (https://www.mathworks.com/matlabcentral/fileexchange/51236-subplot), MATLAB Central File Exchange. Retrieved February 20, 2020. 

% This script requires functions from:
% Fred Gruber (2020). Annotation pinned to axes (https://www.mathworks.com/matlabcentral/fileexchange/32536-annotation-pinned-to-axes), MATLAB Central File Exchange. Retrieved February 21, 2020. 

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

TFinal = 5000;
TInit = 0;
TStep = 1/52;
TSpan = [TInit:TStep:TFinal];
TRecord = [(TFinal - 100):TStep:TFinal];

TimeVect = 0:1/52:2;

addpath('../Functions');
prompt = 'Enter path to functions to annotate axes (by Fred Gruber): '
path_annotation_pinned = input(prompt, 's');
addpath(path_annotation_pinned);

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-6, 'NonNegative', [1:2]);

PeriodCol = [228 26 28; 55 126 184; 77 175 74; 152 78 163; 255 127 0; 255 255 51; 166 86 40; 247 129 191; 153 153 153]/255;

% type I functional response
a = 16*pi;
SubPlotCount = 1;
for theta = [0 0.5 1]
  for epsilon = 0.5:0.001:1
		for Rep = 1:10
			disp([1 theta epsilon Rep]);
			IC = rand(2, 1); % pick initial conditions between 0.1 and 1
			[T, Y] = odeIntegSeason(IC, TFinal, TRecord, @(t, y) seasonKill2Sp_Type1(t, y, true), options);
			Y = Y((T >= TFinal - 100), :); T = T(T>= TFinal -100);
			Index10Y = find(ismember(T, [(TFinal - 9):TFinal]));

			% add points to graph
			subplot(3, 3, 1 + 3 * (SubPlotCount - 1));
			Col = 'k';
			PEst = max(arrayfun(@(X) periodAnalysis(T, Y(:, X)), [1 2]), [], 'includenan');
			if ~isnan(PEst)
				Col = PeriodCol(PEst, :);
			end
			NPts = size(unique(Y(Index10Y, 2)), 1);
	 	 	scatter(epsilon * ones(NPts, 1), unique(Y(Index10Y, 2)), 7.5, Col, '.'); hold on;
			clear T Y
	  end
	end

  SubPlotCount = SubPlotCount + 1;
	set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 8], 'PaperUnits', 'Inches', 'PaperSize', [10, 8]);
	print('../Figures/fig_bifurcdiagr_controlvar_sinus2rect_forcing', '-dpng');
	savefig('../Outputs/fig_bifurcdiagr_controlvar_sinus2rect_forcing.fig');
end

% type II functional response
a = 12*pi;
SubPlotCount = 1;
for theta = [0 0.5 1]
  for epsilon = 0.01:0.001:0.5
		for Rep = 1:10
			disp([2 theta epsilon Rep]);
			IC = rand(2, 1); % pick initial conditions between 0.1 and 1
			[T, Y] = odeIntegSeason(IC, TFinal, TRecord, @(t, y) seasonKill2Sp_Type2(t, y, true), options);
			Y = Y((T >= TFinal - 100), :); T = T(T>= TFinal -100);
			Index10Y = find(ismember(T, [(TFinal - 9):TFinal]));

			% add points to graph
			subplot(3, 3, 2 + 3 * (SubPlotCount - 1));
			Col = 'k';
			PEst = max(arrayfun(@(X) periodAnalysis(T, Y(:, X)), [1 2]), [], 'includenan');
			if ~isnan(PEst)
				Col = PeriodCol(PEst, :);
			end
			NPts = size(unique(Y(Index10Y, 2)), 1);
	 	 	scatter(epsilon * ones(NPts, 1), unique(Y(Index10Y, 2)), 7.5, Col, '.'); hold on;
			clear T Y
    end
	end

  SubPlotCount = SubPlotCount + 1;

	set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 8], 'PaperUnits', 'Inches', 'PaperSize', [10, 8]);
	print('../Figures/fig_bifurcdiagr_controlvar_sinus2rect_forcing', '-dpng');
	savefig('../Outputs/fig_bifurcdiagr_controlvar_sinus2rect_forcing.fig');
end

% type II functional response
a = 16*pi;
SubPlotCount = 1;
for theta = [0 0.5 1]
	for epsilon = 0.2:0.01:0.7
		for Rep = 1:10
			disp([2 a theta epsilon Rep]);
			IC = rand(2, 1); % pick initial conditions between 0.1 and 1
			[T, Y] = odeIntegSeason(IC, TFinal, TRecord, @(t, y) seasonKill2Sp_Type2(t, y, true), options);
			Y = Y((T >= TFinal - 100), :); T = T(T>= TFinal -100);
			Index10Y = find(ismember(T, [(TFinal - 9):TFinal]));

			% add points to graph
			subplot(3, 3, 3 * SubPlotCount);
			Col = 'k';
			PEst = max(arrayfun(@(X) periodAnalysis(T, Y(:, X)), [1 2]), [], 'includenan');
			if ~isnan(PEst)
				Col = PeriodCol(PEst, :);
			end

			NPts = size(unique(Y(Index10Y, 2)), 1);
	 	 	scatter(epsilon * ones(NPts, 1), unique(Y(Index10Y, 2)), 7.5, Col, '.'); hold on;
			clear T Y
	  end
	end
  SubPlotCount = SubPlotCount + 1;

	set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 8], 'PaperUnits', 'Inches', 'PaperSize', [10, 8]);
	print('../Figures/fig_bifurcdiagr_controlvar_sinus2rect_forcing', '-dpng');
	savefig('../Outputs/fig_bifurcdiagr_controlvar_sinus2rect_forcing.fig');
end

for ii = 1:3:7
	subplot(3, 3, ii); hold on; xlim([0.5 1]); ylim([0 1]);
	xlabel("\epsilon", 'FontSize', 8); ylabel("Prey densities", 'FontSize', 8);
	text(min(xlim) + 0.1 * diff(xlim), max(ylim) - 0.1*diff(ylim), char(ii + 64));
	set(gca, 'FontSize', 8);
end
for ii = 2:3:8
	subplot(3, 3, ii); hold on; xlim([0 0.5]); ylim([0 0.8]);
	xlabel("\epsilon", 'FontSize', 8); ylabel("Prey densities", 'FontSize', 8);
	text(min(xlim) + 0.1 * diff(xlim), max(ylim) - 0.1*diff(ylim), char(ii + 64));
	set(gca, 'FontSize', 8);
end
for ii = 3:3:9
	subplot(3, 3, ii); hold on; xlim([0.2 0.7]); ylim([0 1]);
	xlabel("\epsilon", 'FontSize', 8); ylabel("Prey densities", 'FontSize', 8);
	text(min(xlim) + 0.1 * diff(xlim), max(ylim) - 0.1*diff(ylim), char(ii + 64));
	set(gca, 'FontSize', 8);
end
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 8], 'PaperUnits', 'Inches', 'PaperSize', [10, 8]);

subplot(3, 3, 2); hold on;
for Col = 1:size(PeriodCol, 1)
	p(Col) = plot([nan nan], '.', 'color', PeriodCol(Col, :)); hold on;
end
p(Col + 1) = plot([nan nan], '.', 'color', 'k'); hold on;
legend(p(:), {'1', '2', '3', '4', '5', '6', '7', '8', '9', '>9'}, 'Location', 'southeast');

subplot(3, 3, 1); hold on;
title('Type I, $\bar{a} = 16\pi$', 'Interpreter','Latex');
subplot(3, 3, 2); hold on;
title('Type II, $\bar{a} = 12\pi$', 'Interpreter','Latex');
subplot(3, 3, 3); hold on;
title('Type II, $\bar{a} = 16\pi$', 'Interpreter','Latex');

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 8], 'PaperUnits', 'Inches', 'PaperSize', [10, 8]);
print('../Figures/fig_bifurcdiagr_controlvar_sinus2rect_forcing', '-dpng');
savefig('../Outputs/fig_bifurcdiagr_controlvar_sinus2rect_forcing.fig');

% add vertical lines with the net forcing amplitude = EE * StdQ_TT1 / StdQ_TT
ForcingSignal_TT1 = arrayfun(@(X) forcingSignal(X, 1, 1), TimeVect);
StdQ_TT1 = std(ForcingSignal_TT1);
ForcingSignal_TT0 = arrayfun(@(X) forcingSignal(X, 1, 0), TimeVect);
StdQ_TT0 = std(ForcingSignal_TT0);
ForcingSignal_TT05 = arrayfun(@(X) forcingSignal(X, 1, 0.5), TimeVect);
StdQ_TT05 = std(ForcingSignal_TT05);

subplot(3, 3, 7); hold on;
NetEE1_TT0 = StdQ_TT1 / StdQ_TT0;
NetEE1_TT05 = StdQ_TT1 / StdQ_TT05;
line([NetEE1_TT0 NetEE1_TT0], [0 1], 'Color', 'k');
annotation_pinned('textarrow', [NetEE1_TT0-0.05 NetEE1_TT0], [0.8 0.9], 'String', '$$\tilde{\epsilon}(0, 1)$$', 'Interpreter', 'LaTeX');
line([NetEE1_TT05 NetEE1_TT05], [0 1], 'Color', 'k', 'LineStyle', '--');
annotation_pinned('textarrow', [NetEE1_TT05-0.05 NetEE1_TT05], [0.9 0.8], 'String', '$$\tilde{\epsilon}(0.5, 1)$$', 'Interpreter', 'LaTeX');

subplot(3, 3, 8); hold on;
NetEE05_TT0 = 0.5 * StdQ_TT1 / StdQ_TT0;
NetEE05_TT05 = 0.5 * StdQ_TT1 / StdQ_TT05;
line([NetEE05_TT0 NetEE05_TT0], [0 1], 'Color', 'k');
annotation_pinned('textarrow', [NetEE05_TT0-0.05 NetEE05_TT0], [0.1 0.2], 'String', '$$\tilde{\epsilon}(0, 0.5)$$', 'Interpreter', 'LaTeX');
line([NetEE05_TT05 NetEE05_TT05], [0 1], 'Color', 'k', 'LineStyle', '--');
annotation_pinned('textarrow', [NetEE05_TT05-0.025 NetEE05_TT05], [0.5 0.3], 'String', '$$\tilde{\epsilon}(0.5, 0.5)$$', 'Interpreter', 'LaTeX');

subplot(3, 3, 9); hold on;
NetEE07_TT0 = 0.7 * StdQ_TT1 / StdQ_TT0;
NetEE07_TT05 = 0.7 * StdQ_TT1 / StdQ_TT05;
line([NetEE07_TT0 NetEE07_TT0], [0 1], 'Color', 'k');
annotation_pinned('textarrow', [NetEE07_TT0-0.05 NetEE07_TT0], [0.4 0.3], 'String', '$$\tilde{\epsilon}(0, 0.7)$$', 'Interpreter', 'LaTeX');
line([NetEE07_TT05 NetEE07_TT05], [0 1], 'Color', 'k', 'LineStyle', '--');
annotation_pinned('textarrow', [NetEE07_TT05+0.03 NetEE07_TT05], [0.3 0.5], 'String', '$$\tilde{\epsilon}(0.5, 0.7)$$', 'Interpreter', 'LaTeX');

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 8], 'PaperUnits', 'Inches', 'PaperSize', [10, 8]);
print('../Figures/fig_bifurcdiagr_controlvar_sinus2rect_forcing', '-dpng');
savefig('../Outputs/fig_bifurcdiagr_controlvar_sinus2rect_forcing.fig');

