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

TFinal = 2000;
TInit = 0;
TStep = 1/52;
TSpan = [TInit:TStep:TFinal];
TRecord = [(TFinal - 50):TStep:TFinal];

addpath('../Functions');

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-6, 'NonNegative', [1:2]);

PeriodCol = [228 26 28; 55 126 184; 77 175 74; 152 78 163; 255 127 0; 255 255 51; 166 86 40; 247 129 191; 153 153 153]/255;

% testing the effect of epsilon with a sinusoidal signal
aValues = [8 12 16];
for aIndex = 1:3
	a = aValues(aIndex)*pi;
  for epsilon = 0.01:0.01:1
		for Rep = 1:10
			IC = 0.1 + 0.9.*rand(2, 1); % pick initial conditions between 0.1 and 1
	
			% type 1
			[T, Y] = ode15s(@(t, y) seasonKill2Sp_Type1(t, y, false), [0 TRecord], IC, options);
			Y = Y((T >= TFinal - 100), :); T = T(T>= TFinal -100);
			Index10Y = find(ismember(T, [(TFinal - 9):TFinal]));
			disp([1 aIndex epsilon Rep]);
			
			% estimating period with FFT
			PEstFFT = periodFFT(T, Y);
			dlmwrite('../Outputs/period_fft_t1.csv', [aIndex epsilon Rep PEstFFT], 'delimiter', ',', '-append');

			% add points to graph
			figure(1);
			subplot(2, 3, aIndex); hold on;
			if (epsilon == 0.01) & (Rep == 1)
				xlim([0 1]); ylim([0 1]);
			  xlabel("\epsilon", 'FontSize', 8); ylabel("Prey density", 'FontSize', 8);
			  title(strcat('a = ',  int2str(aValues(aIndex)), '\pi'), 'FontSize', 8);
				text(0.05*1, max(ylim)*(1 - 0.1), char(aIndex + 64));
			  set(gca, 'FontSize', 8);
			end

			Col = 'k';
			PEst = max(arrayfun(@(X) periodAnalysis(T, Y(:, X)), [1 2]), [], 'includenan');
			if ~isnan(PEst)
				Col = PeriodCol(PEst, :);
			end
			NPts = size(unique(Y(Index10Y, 2)), 1);
	 	 	scatter(epsilon * ones(NPts, 1), unique(Y(Index10Y, 2)), 7.5, Col, '.'); hold on;


			% add predator's points to another figure
			figure(2);
			subplot(2, 3, aIndex); hold on;
			if (epsilon == 0.01) & (Rep == 1)
				xlim([0 1]); ylim([0 0.45]);
			  xlabel("\epsilon", 'FontSize', 8); ylabel("Predator density", 'FontSize', 8);
			  title(strcat('a = ',  int2str(aValues(aIndex)), '\pi'), 'FontSize', 8);
				text(min(xlim) + 0.05*diff(xlim), max(ylim) - 0.1*diff(ylim), char(aIndex + 64));
			  set(gca, 'FontSize', 8);
			end
			NPts = size(unique(Y(Index10Y, 1)), 1);
	 	 	scatter(epsilon * ones(NPts, 1), unique(Y(Index10Y, 1)), 7.5, Col, '.'); hold on;
			clear T Y

			% type 2
			[T, Y] = ode15s(@(t, y) seasonKill2Sp_Type2(t, y, false), [0 TRecord], IC, options);
			Y = Y((T >= TFinal - 100), :); T = T(T>= TFinal -100);
			Index10Y = find(ismember(T, [(TFinal - 9):TFinal]));
			disp([2 aIndex epsilon Rep]);

			% estimating period with FFT
			PEstFFT = periodFFT(T, Y);
			dlmwrite('../Outputs/period_fft_t2.csv', [aIndex epsilon Rep PEstFFT], 'delimiter', ',', '-append');

			% add points to graph
			figure(1);
			subplot(2, 3, aIndex + 3); hold on;
			if (epsilon == 0.01) & (Rep == 1)
				xlim([0 1]); ylim([0 1]);
			  xlabel("\epsilon", 'FontSize', 8); ylabel("Prey density", 'FontSize', 8);
				text(0.05*1, max(ylim)*(1 - 0.1), char(aIndex + 3 + 64));
			  set(gca, 'FontSize', 8);
			end
			Col = 'k';
			PEst = max(arrayfun(@(X) periodAnalysis(T, Y(:, X)), [1 2]), [], 'includenan');
			if ~isnan(PEst)
				Col = PeriodCol(PEst, :);
			end
			NPts = size(unique(Y(Index10Y, 2)), 1);
	 	 	scatter(epsilon * ones(NPts, 1), unique(Y(Index10Y, 2)), 7.5, Col, '.'); hold on;


			% add predator's points to another figure
			figure(2);
			subplot(2, 3, aIndex + 3); hold on;
			if (epsilon == 0.01) & (Rep == 1)
				xlim([0 1]); ylim([0 0.45]);
			  xlabel("\epsilon", 'FontSize', 8); ylabel("Predator density", 'FontSize', 8);
				text(min(xlim) + 0.05*diff(xlim), max(ylim) - 0.1*diff(ylim), char(aIndex + 3 + 64));
			  set(gca, 'FontSize', 8);
			end
			NPts = size(unique(Y(Index10Y, 1)), 1);
	 	 	scatter(epsilon * ones(NPts, 1), unique(Y(Index10Y, 1)), 7.5, Col, '.'); hold on;
			clear T Y
    end
	end
	figure(1);	
	set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 5], 'PaperUnits', 'Inches', 'PaperSize', [10, 5]);
	print('../Figures/fig_prey_bifurcdiagr_sinus_forcing', '-dpng');
	savefig('../Outputs/fig_prey_bifurcdiagr_sinus_forcing.fig');

	figure(2);
	set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 5], 'PaperUnits', 'Inches', 'PaperSize', [10, 5]);
	print('../Figures/fig_pred_bifurcdiagr_sinus_forcing', '-dpng');
	savefig('../Outputs/fig_pred_bifurcdiagr_sinus_forcing.fig');
end

figure(1)
subplot(2, 3, 4); hold on;
for Col = 1:size(PeriodCol, 1)
	p(Col) = plot([nan nan], '.', 'color', PeriodCol(Col, :)); hold on;
end
p(Col + 1) = plot([nan nan], '.', 'color', 'k'); hold on;
legend(p(:), {'1', '2', '3', '4', '5', '6', '7', '8', '9', '>9'}, 'Location', 'southeast');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8, 4], 'PaperUnits', 'Inches', 'PaperSize', [8, 4]);
print('../Figures/fig_prey_bifurcdiagr_sinus_forcing', '-dpng');
savefig('../Outputs/fig_prey_bifurcdiagr_sinus_forcing.fig');
close(1);

figure(2)
subplot(2, 3, 4); hold on;
for Col = 1:size(PeriodCol, 1)
	p(Col) = plot([nan nan], '.', 'color', PeriodCol(Col, :)); hold on;
end
p(Col + 1) = plot([nan nan], '.', 'color', 'k'); hold on;
legend(p(:), {'1', '2', '3', '4', '5', '6', '7', '8', '9', '>9'}, 'Location', 'northeast');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 5], 'PaperUnits', 'Inches', 'PaperSize', [10, 5]);
print('../Figures/fig_pred_bifurcdiagr_sinus_forcing', '-dpng');
savefig('../Outputs/fig_pred_bifurcdiagr_sinus_forcing.fig');
close(2);
