% Because the simulations are long to run, we recommend to segment the following script into smaller subsets so as to generate each subplot separately.
% Figures can be gathered within the same multi-panel figure afterwards with the help of 'SubPlot':
% Farhad Sedaghati (2020). SubPlot (https://www.mathworks.com/matlabcentral/fileexchange/51236-subplot), MATLAB Central File Exchange. Retrieved February 20, 2020. 

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

EXTINCT_THRS = 1e-5; % extinction threshold

rng(2019);

TFinal = 5000;
TInit = 0;
TStep = 1/52;
TSpan = [TInit:TStep:TFinal];
TRecord = [(TFinal - 50):TStep:TFinal];

addpath('../Functions');

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-6, 'NonNegative', [1:2]);

PeriodCol = [228 26 28; 55 126 184; 77 175 74; 152 78 163; 255 127 0; 255 255 51; 166 86 40; 247 129 191; 153 153 153]/255;

% type I functional response
a = 16*pi;
SubPlotCount = 1;
for theta = [0 0.5 1]
  for epsilon = 0.5:0.001:1
		for Rep = 1:10
			disp([1 a theta epsilon Rep]);

			IC = rand(2, 1); % pick initial conditions between 0.1 and 1
			[T, Y] = odeIntegSeason(IC, TFinal, TRecord, @(t, y) seasonKill2Sp_Type1(t, y, false), options);
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
	print('../Figures/fig_bifurcdiagr_sinus2rect_forcing', '-dpng');
	savefig('../Outputs/fig_bifurcdiagr_sinus2rect_forcing.fig');
end

% type II functional response
a = 12*pi;
SubPlotCount = 1;
for theta = [0 0.5 1]
  for epsilon = 0:0.001:0.5
		for Rep = 1:10
			disp([2 a theta epsilon Rep]);

			IC = rand(2, 1); % pick initial conditions between 0.1 and 1
			[T, Y] = odeIntegSeason(IC, TFinal, TRecord, @(t, y) seasonKill2Sp_Type2(t, y, false), options);
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
	print('../Figures/fig_bifurcdiagr_sinus2rect_forcing', '-dpng');
	savefig('../Outputs/fig_bifurcdiagr_sinus2rect_forcing.fig');
end

% type II functional response
a = 16*pi;
SubPlotCount = 1;
for theta = [0 0.5 1]
	for epsilon = 0.2:0.01:0.7
		for Rep = 1:10
			disp([2 a theta epsilon Rep]);
			IC = rand(2, 1); % pick initial conditions between 0.1 and 1
			[T, Y] = odeIntegSeason(IC, TFinal, TRecord, @(t, y) seasonKill2Sp_Type2(t, y, false), options);
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
	print('../Figures/fig_bifurcdiagr_sinus2rect_forcing', '-dpng');
	savefig('../Outputs/fig_bifurcdiagr_sinus2rect_forcing.fig');
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
savefig('../Outputs/fig_bifurcdiagr_sinus2rect_forcing.fig');
print('../Figures/fig_bifurcdiagr_sinus2rect_forcing', '-dpng');

