function PeriodEst = periodFFT(TVect, YVect)
	N = length(TVect);
	if rem(N, 2)
		TVect = TVect(2:end); YVect = YVect(2:end, :);
		N = length(TVect);
	end
	YCorr = YVect - mean(YVect);
	HannWindow = hann(N);
	Yfft = fft(YCorr.*[HannWindow HannWindow]);
	Yfft = abs(Yfft.^2); % raw power spectrum density
	Yfft = Yfft(1:(1 + N/2));
	Fs = 52;
	FreqScale = (0:N/2) * Fs/N; % frequency scale
	[Vfft, k] = max(Yfft); % find maximum
	PeriodEst = 1/FreqScale(k);
end
