% same as periodAnalysis2(), but return 'nan' if no period is detected
function PeriodEst = periodAnalysis2(TVect, YVect)
	TCandidates = [1:9]; IsPeriodic = zeros(size(TCandidates));
	for P = TCandidates
		PVect = find(ismember(TVect, sort([max(TVect):-P:(max(TVect)-100)])));
		MinY = min(YVect(PVect)); MaxY = max(YVect(PVect));
		MeanY = mean(YVect(PVect));
		if (MaxY - MinY) < (2.5/100 * MeanY)
			IsPeriodic(TCandidates == P) = true;
		end
	end
	PeriodEst = min(find(IsPeriodic));
	if isempty(PeriodEst)
		PeriodEst = nan;
	end
end
