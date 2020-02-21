% In this script, I run a bifurcation analysis to find a Neimark-Sacker bifurcation if it exists when the parameterisation is the same as in Taylor et al. (2013, J. Math. Biol.) but forcing on the discovery rate.

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
a = 4*pi/0.3; % discovery rate

rng(2019);

TFinal = 2000;
TInit = 0;
TStep = 1/52;
TSpan = [TInit:TStep:TFinal];
TRecord = [(TFinal - 100):TStep:TFinal];

addpath('../Functions');

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-6, 'NonNegative', [1:2]);
PeriodCol = [228 26 28; 55 126 184; 77 175 74; 152 78 163; 255 127 0; 255 255 51; 166 86 40; 247 129 191; 153 153 153]/255;

for epsilon = 0.01:0.01:1
	for Rep = 1:10
		IC = rand(2, 1); % pick initial conditions between 0.1 and 1
		[T, Y] = ode15s(@(t, y) seasonKill2Sp_Type2(t, y, false), [0 TRecord], IC, options);
		Index10Y = find(ismember(T, [(TFinal - 99):TFinal]));
		disp([epsilon Rep]);

		Col = 'k';
		PEst = periodAnalysis(T(T >= TFinal - 100), Y((T >= TFinal - 100), 2));
		if ~isempty(PEst)
			Col = PeriodCol(PEst, :);
		end
		subplot(1, 2, 1); hold on;
		scatter(epsilon * ones(10, 1), unique(Y(Index10Y, 2)), 10, Col, 'filled'); hold on;
		subplot(1, 2, 2); hold on;
		scatter(epsilon * ones(10, 1), unique(Y(Index10Y, 1)), 10, Col, 'filled'); hold on;
	end
end
subplot(1, 2, 1); hold on;
ylim([0 1]); xlim([0 1]);
xlabel("\epsilon", 'FontSize', 8); ylabel("Prey density", 'FontSize', 8);
set(gca, 'FontSize', 8);
subplot(1, 2, 2); hold on;
ylim([0 0.4]); xlim([0 1]);
xlabel("\epsilon", 'FontSize', 8); ylabel("Predator density", 'FontSize', 8);
set(gca, 'FontSize', 8);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 4], 'PaperUnits', 'Inches', 'PaperSize', [10, 4]);
print('../Figures/fig_bifurcdiagr_zoom_sinus_forcing', '-dpng');
savefig('../Outputs/fig_bifurcdiagr_zoom_sinus_forcing.fig');


for epsilon = 0.01:0.01:1
	for Rep = 1:10
		IC = rand(2, 1); % pick initial conditions between 0.1 and 1
		[T, Y] = ode15s(@(t, y) seasonRepro2Sp_Type2(t, y, false), [0 TRecord], IC, options);
		Index10Y = find(ismember(T, [(TFinal - 9):TFinal]));
		disp([epsilon Rep]);

		Col = 'k';
		PEst = periodAnalysis(T(T >= TFinal - 100), Y((T >= TFinal - 100), 2));
		if ~isempty(PEst)
			Col = PeriodCol(PEst, :);
		end
		subplot(1, 2, 1); hold on;
		scatter(epsilon * ones(10, 1), unique(Y(Index10Y, 2)), 10, Col, 'filled'); hold on;
		subplot(1, 2, 2); hold on;
		scatter(epsilon * ones(10, 1), unique(Y(Index10Y, 1)), 10, Col, 'filled'); hold on;
	end
end
subplot(1, 2, 1); hold on;
ylim([0 1]); xlim([0 1]);
xlabel("\epsilon", 'FontSize', 8); ylabel("Prey density", 'FontSize', 8);
set(gca, 'FontSize', 8);
subplot(1, 2, 2); hold on;
ylim([0 0.4]); xlim([0 1]);
xlabel("\epsilon", 'FontSize', 8); ylabel("Predator density", 'FontSize', 8);
set(gca, 'FontSize', 8);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 4], 'PaperUnits', 'Inches', 'PaperSize', [10, 4]);
print('../Figures/BifurcDiagEpsilon_SeasonalRepro_20190620', '-dpng');
savefig('../Outputs/BifurcDiagEpsilon_SeasonalRepro_20190620.fig');

% zoom on epsilon = {0.3:0.5}, a = 16*pi
figure();
a = 16*pi;
for epsilon = 0.3:0.001:0.5
	for Rep = 1:20
		IC = rand(2, 1); % pick initial conditions between 0.1 and 1
		[T, Y] = ode15s(@(t, y) seasonKill2Sp_Type2(t, y, false), [0 TRecord], IC, options);
		Index10Y = find(ismember(T, [(TFinal - 9):TFinal]));
		disp([epsilon Rep]);
		Col = 'k';
		PEst_Prey = periodAnalysis(T(T >= TFinal - 100), Y((T >= TFinal - 100), 2));
		PEst_Pred = periodAnalysis(T(T >= TFinal - 100), Y((T >= TFinal - 100), 1));
		PEst = max([PEst_Prey PEst_Pred]);
		if ~isempty(PEst)
			Col = PeriodCol(PEst, :);
		end
		subplot(1, 2, 1); hold on;
		scatter(epsilon * ones(10, 1), unique(Y(Index10Y, 2)), 10, Col, 'filled'); hold on;
		subplot(1, 2, 2); hold on;
		scatter(epsilon * ones(10, 1), unique(Y(Index10Y, 1)), 10, Col, 'filled'); hold on;
	end
end
subplot(1, 2, 2); hold on;
xlabel("\epsilon", 'FontSize', 8); ylabel("Predator density", 'FontSize', 8);
set(gca, 'FontSize', 8);
xlim([0.3 0.5]); ylim([0 0.4]); text(min(xlim) + 0.1 * diff(xlim), max(ylim) - 0.1*diff(ylim), char(2 + 64));
subplot(1, 2, 1); hold on;
xlabel("\epsilon", 'FontSize', 8); ylabel("Prey density", 'FontSize', 8);                                    
set(gca, 'FontSize', 8);                                                                                     
xlim([0.3 0.5]); ylim([0 1]); text(min(xlim) + 0.1 * diff(xlim), max(ylim) - 0.1*diff(ylim), char(1 + 64));  
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8, 4], 'PaperUnits', 'Inches', 'PaperSize', [8, 4]);
print('../Figures/BifurcDiagEpsilon_Type2_a16pi_ZoomE03_20190621', '-dpng');
savefig('../Outputs/BifurcDiagEpsilon_Type2_a16pi_ZoomE03_20190621.fig');

% zoom on epsilon = {0.64:0.75}, a = 16*pi
a = 16*pi;
for epsilon = 0.64:0.001:0.75
	for Rep = 1:20
		IC = rand(2, 1); % pick initial conditions between 0.1 and 1
		[T, Y] = ode15s(@(t, y) seasonKill2Sp_Type2(t, y, false), [0 TRecord], IC, options);
		Index10Y = find(ismember(T, [(TFinal - 9):TFinal]));
		disp([epsilon Rep]);
		Col = 'k';
		PEst = max(arrayfun(@(X) periodAnalysis(T, Y(:, X)), [1 2]), [], 'includenan');
		if ~isnan(PEst)
			Col = PeriodCol(PEst, :);
		end
		subplot(1, 2, 1); hold on;
		scatter(epsilon * ones(10, 1), unique(Y(Index10Y, 2)), 7.5, Col, '.'); hold on;
		subplot(1, 2, 2); hold on;
		scatter(epsilon * ones(10, 1), unique(Y(Index10Y, 1)), 7.5, Col, '.'); hold on;
	end
end
subplot(1, 2, 2); hold on;
xlabel("\epsilon", 'FontSize', 8); ylabel("Predator density", 'FontSize', 8);
set(gca, 'FontSize', 8);
xlim([0.64 0.75]); ylim([0 0.4]); text(min(xlim) + 0.1 * diff(xlim), max(ylim) - 0.1*diff(ylim), char(2 + 64));
subplot(1, 2, 1); hold on;
xlabel("\epsilon", 'FontSize', 8); ylabel("Prey density", 'FontSize', 8);                                    
set(gca, 'FontSize', 8);                                                                                     
xlim([0.64 0.75]); ylim([0 1]); text(min(xlim) + 0.1 * diff(xlim), max(ylim) - 0.1*diff(ylim), char(1 + 64));  
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8, 4], 'PaperUnits', 'Inches', 'PaperSize', [8, 4]);
print('../Figures/BifurcDiagEpsilon_Type2_a16pi_ZoomE064_20190621', '-dpng');
savefig('../Outputs/BifurcDiagEpsilon_Type2_a16pi_ZoomE064_20190621.fig');

% zoom on epsilon = {0:0.1}, a = 12*pi
figure();
a = 12*pi;
for epsilon = 0:0.001:0.1
	for Rep = 1:20
		IC = rand(2, 1); % pick initial conditions between 0.1 and 1
		[T, Y] = ode15s(@(t, y) seasonKill2Sp_Type2(t, y, false), [0 TRecord], IC, options);
		Index10Y = find(ismember(T, [(TFinal - 9):TFinal]));
		disp([epsilon Rep]);
		Col = 'k';
		PEst_Prey = periodAnalysis(T(T >= TFinal - 100), Y((T >= TFinal - 100), 2));
		PEst_Pred = periodAnalysis(T(T >= TFinal - 100), Y((T >= TFinal - 100), 1));
		PEst = max([PEst_Prey PEst_Pred]);
		if ~isempty(PEst)
			Col = PeriodCol(PEst, :);
		end
		subplot(1, 2, 1); hold on;
		scatter(epsilon * ones(10, 1), unique(Y(Index10Y, 2)), 10, Col, 'filled'); hold on;
		subplot(1, 2, 2); hold on;
		scatter(epsilon * ones(10, 1), unique(Y(Index10Y, 1)), 10, Col, 'filled'); hold on;
	end
end

a = 16*pi; epsilon = 0.64;
IC = rand(2, 1); % pick initial conditions between 0.1 and 1
[T, Y] = ode15s(@(t, y) seasonKill2Sp_Type2(t, y, false), [0 TRecord], IC, options);
Y = Y((T > (TFinal - 100)), :); T = T(T > (TFinal - 100));
PEst = periodAnalysis(T, Y(:, 2)); PEst = periodAnalysis(T, Y(:, 1));

epsilon = 0.65;
[T, Y] = ode15s(@(t, y) seasonKill2Sp_Type2(t, y, false), [0 TRecord], IC, options);
Y = Y((T > (TFinal - 100)), :); T = T(T > (TFinal - 100));
PEst = periodAnalysis(T, Y(:, 2)); PEst = periodAnalysis(T, Y(:, 1));
