% this script draws a 2 by 2 figure panel with 1) a sinusoidal forcing signal and a rectangular one, with each recorded point highlighted with a symbol; 2) the amount of forcing exerted during the predation season; 3) the variance of the signal for different levels of rectangularity; and 4) the amount of forcing once the forcing signal is corrected.

% This script requires functions from:
% Fred Gruber (2020). Annotation pinned to axes (https://www.mathworks.com/matlabcentral/fileexchange/32536-annotation-pinned-to-axes), MATLAB Central File Exchange. Retrieved February 21, 2020. 

clc
clear all
close all

addpath('../Functions');

prompt = 'Enter path to functions to annotate axes (by Fred Gruber): '
path_annotation_pinned = input(prompt, 's');
addpath(path_annotation_pinned);

ThetaCol = [59 154 178; 242 26 0]/255;

TimeVect = 0:1/52:2;

subplot(2, 3, 1);
set(gca, 'FontSize', 8);
OPSize = get(gca, 'OuterPosition');
TISize = get(gca, 'TightInset');
SPSize = get(gca, 'Position'); %close all;
subplot('Position', [0.5-OPSize(3)+(SPSize(1)-OPSize(1)) SPSize(2) SPSize(3) SPSize(4)]);
EE = 1;
TT = 1; ForcingSignal_T1 = arrayfun(@(X) forcingSignal(X, EE, TT), TimeVect);
TT = 0;ForcingSignal_T0 = arrayfun(@(X) forcingSignal(X, EE, TT), TimeVect);
line([0 2], [0 0], 'LineStyle', '--', 'Color', 'k'); hold on;
plot(TimeVect, ForcingSignal_T1, '-o', 'color', ThetaCol(2, :), 'MarkerFaceColor', ThetaCol(2, :), 'MarkerSize', 3); hold on;
plot(TimeVect, ForcingSignal_T0, '-o', 'color', ThetaCol(1, :), 'MarkerFaceColor', ThetaCol(1, :), 'MarkerSize', 3);
ylim([-1.2 1.2]);
set(gca, 'FontSize', 8);
xlabel('Time', 'FontSize', 8); ylabel('s(t)', 'FontSize', 8);
text(min(xlim) + 0.05 * diff(xlim), max(ylim) - 0.1*diff(ylim), char(1 + 64));
	
subplot('Position', [0.5-OPSize(3)+(SPSize(1)-OPSize(1))+0.2808+OPSize(3)-SPSize(3)-SPSize(1)+OPSize(1) SPSize(2) SPSize(3) SPSize(4)]);
ThetaCol = [59 154 178; 83 165 185; 107 177 193; 143 187 165; 189 195 103; 235 204 42; 230 192 25; 227 180 8; 228 145 0; 235 85 0; 242 26 0]/255;
ThetaCount = 1;
for TT = 0:0.1:1
	VarForcing = [];
	for EE = 0.1:0.1:1
		ForcingSignal_EE_TT = arrayfun(@(X) forcingSignal(X, EE, TT), TimeVect);
		VarForcing_EE_TT = var(ForcingSignal_EE_TT);
		VarForcing = [VarForcing VarForcing_EE_TT];
	end
	plot([0.1:0.1:1], VarForcing, '-o', 'color', ThetaCol(ThetaCount, :), 'MarkerFaceColor', ThetaCol(ThetaCount, :), 'MarkerSize', 3); hold on;
	ThetaCount = ThetaCount + 1;
end
box off;
set(gca, 'FontSize', 8);
xlabel('\epsilon', 'FontSize', 8); ylabel('Variance of s(t)', 'FontSize', 8);
hold on;
annotation_pinned('textarrow', [0.6 0.7], [0.6 var(arrayfun(@(X) forcingSignal(X, 0.7, 0), TimeVect))], 'String', '\theta = 0');
annotation_pinned('textarrow', [0.7 0.6], [0.1 var(arrayfun(@(X) forcingSignal(X, 0.6, 1), TimeVect))], 'String', '\theta = 1');
text(min(xlim) + 0.05 * diff(xlim), max(ylim) - 0.1*diff(ylim), char(2 + 64));
set(gca, 'FontSize', 8);

subplot(2, 3, 4);
ThetaCol = [59 154 178; 242 26 0]/255;
EE = 1;
TT = 1; ForcingSignal_T1 = arrayfun(@(X) forcingSignal(X, EE, TT), TimeVect);
TT = 0;ForcingSignal_T0 = arrayfun(@(X) forcingSignal(X, EE, TT), TimeVect);
ForcingSignal_1 = arrayfun(@(X) forcingSignal(X, 1, 1), TimeVect);
StdQ_1 = std(ForcingSignal_1);
ForcingSignal_0 = arrayfun(@(X) forcingSignal(X, 1, 0), TimeVect);
StdQ_0 = std(ForcingSignal_0);
line([0 2], [0 0], 'LineStyle', '--', 'Color', 'k'); hold on;
plot(TimeVect, ForcingSignal_T1, '-o', 'color', ThetaCol(2, :), 'MarkerFaceColor', ThetaCol(2, :), 'MarkerSize', 3); hold on;
plot(TimeVect, ForcingSignal_T0* StdQ_1 / StdQ_0, '-o', 'color', ThetaCol(1, :), 'MarkerFaceColor', ThetaCol(1, :), 'MarkerSize', 3);
ylim([-1.2 1.2]);
set(gca, 'FontSize', 8);
xlabel('Time', 'FontSize', 8); ylabel('Corrected s(t)', 'FontSize', 8);
text(min(xlim) + 0.05 * diff(xlim), max(ylim) - 0.1*diff(ylim), char(3 + 64));

subplot(2, 3, 5);
ThetaCol = [59 154 178; 83 165 185; 107 177 193; 143 187 165; 189 195 103; 235 204 42; 230 192 25; 227 180 8; 228 145 0; 235 85 0; 242 26 0]/255;
ThetaCount = 1;
for TT = 0:0.1:1
	VarForcing = [];
	ForcingSignal_TT = arrayfun(@(X) forcingSignal(X, 1, TT), TimeVect);
	StdQ_TT = std(ForcingSignal_TT);
	for EE = 0.1:0.1:1
		ForcingSignal_EE_TT = arrayfun(@(X) forcingSignal(X, EE, TT), TimeVect);
		VarForcing_EE_TT = var(ForcingSignal_EE_TT * StdQ_1 / StdQ_TT);
		VarForcing = [VarForcing VarForcing_EE_TT];
	end
	plot([0.1:0.1:1], VarForcing, '-o', 'color', ThetaCol(ThetaCount, :), 'MarkerFaceColor', ThetaCol(ThetaCount, :), 'MarkerSize', 3); hold on;
	ThetaCount = ThetaCount + 1;
end
box off;
set(gca, 'FontSize', 8);
xlabel('\epsilon', 'FontSize', 8); ylabel('Variance of corrected s(t)', 'FontSize', 8);
text(min(xlim) + 0.05 * diff(xlim), max(ylim) - 0.1*diff(ylim), char(4 + 64));

subplot(2, 3, 6);
ThetaCount = 1;
for TT = 0:0.1:1
	EffForcing = [];
	ForcingSignal_TT = arrayfun(@(X) forcingSignal(X, 1, TT), TimeVect);
	StdQ_TT = std(ForcingSignal_TT);
	for EE = 0.1:0.1:1
		EffForcing_EE_TT = EE * StdQ_1 / StdQ_TT;
		EffForcing = [EffForcing EffForcing_EE_TT];
	end
	plot([0.1:0.1:1], EffForcing, '-o', 'color', ThetaCol(ThetaCount, :), 'MarkerFaceColor', ThetaCol(ThetaCount, :), 'MarkerSize', 3); hold on;
	set(gca, 'FontSize', 8);
	xlabel('\epsilon', 'FontSize', 8); ylabel('Net \epsilon', 'FontSize', 8);
	ThetaCount = ThetaCount + 1;
end
box off;
text(min(xlim) + 0.05 * diff(xlim), max(ylim) - 0.1*diff(ylim), char(5 + 64));
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8, 5], 'PaperUnits', 'Inches', 'PaperSize', [8, 5]);
savefig('../Outputs/fig_describe_forcing.fig');
print('../Figures/fig_describe_forcing', '-dpng');
