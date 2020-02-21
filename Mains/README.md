In this folder, we share several scripts to simulate predator-prey dynamics with a seasonal forcing signal varying in shape and magnitude.

## Simulating predator-prey dynamics in a seasonal environment

* To draw **bifurcation diagrams** of species densities {x(t), y(t)} at the beginning of the predation season, with the magnitude of the forcing signal &#949; being the bifurcation parameter:
	- `bifurcdiagr_sinus_forcing.m` when the forcing signal is fully sinusoidal (i.e., &#952; = 1). The resulting bifurcation diagrams correspond to Figs. 3 and A1 in the *JTB* manuscript.

	- `bifurcdiagr_sinus2rect_forcing.m` when considering different shapes of the forcing signal. With &#952; = {0, 0.5, 1}, the forcing signal is either rectangular, intermediate, or sinusoidal. The resulting bifurcation diagrams correspond to Fig. 5 in the *JTB* manuscript.

	- `bifurcdiagr_controlvar_sinus2rect_forcing.m` when considering different shapes of the forcing signal and controlling the signal variance. For these simulations, signal variance thus remains the same despite changing the shape of the forcing signal. The resulting bifurcation diagrams correspond to Fig. 8 in the *JTB* manuscript.

In these three scripts, bifurcation diagrams are drawn for type I and type II functional responses, and for three values of mean discovery rate: {8&#960;, 12&#960;, 16&#960;}.

* `bifurcdiagr_zoom_sinus_forcing.m`: This script zooms on specific areas of the bifurcation diagram drawed with `bifurcdiagr_sinus_forcing.m`.
The resulting bifurcation diagrams correspond to Figs. A2 and A3 in the *JTB* manuscript.

* `basins_attraction.m`: This script estimates the periodicity of the seasonal predator-prey system for various initial conditions to identify coexisting bassins of attraction.
The resulting bassins of attractions are represented in Fig. A4 in the *JTB* manuscript.

* `time_series.m`: This script produces time series illustrating the bassins of attraction identified with `basins_attraction.m`.
The resulting time series are displayed in Fig. A5 in the *JTB* manuscript.

* To draw **phase portraits** for predator-prey dynamics undergoing various magnitude (&#949; &#8712; ]0, 1] with increment 0.1) and shape (&#952; &#8712; ]0, 1] with increment 0.1) of seasonal forcing. Both type I and type II are simulated within the following scripts. For each mean value of discovery rate considered, phase portraits are displayed in a 3D-space ((&#949;, x(t), y(t))):

	- `phaseportraits.m` when signal variance **is not** controlled. The resulting phase portraits correspond to Fig. 4 in the *JTB* manuscript.

	- `phaseportraits_controlvar.m` when signal variance is controlled. The resulting phase portraits correspond to Fig. 7 in the *JTB* manuscript.

* To draw **stroboscopic maps** of the system for various magnitudes and shapes of seasonal forcings:

	- `stroboscmaps_t1.m` when predation is ruled by a type I functional response for these simulations. The resulting stroboscopic maps are displayed in Fig. B1 in the *JTB* manuscript.

	- `stroboscmaps_t2.m` when predation is ruled by a type II functional response for these simulations. The resulting stroboscopic maps are displayed in Fig. B2 in the *JTB* manuscript.

	- `stroboscmaps_t1_controlvar.m` when predation is ruled by a type I functional response for these simulations, and signal variance is controlled so as to remain the same for various shapes of forcing signals. The resulting stroboscopic maps are displayed in Fig. B3 in the *JTB* manuscript.

	- `stroboscmaps_t2_controlvar.m` when predation is ruled by a type II functional response for these simulations, and signal variance is controlled so as to remain the same for various shapes of forcing signals. The resulting stroboscopic maps are displayed in Fig. B4 in the *JTB* manuscript.

NB: These scripts call the following functions: `seasonKill2Sp_Type1()`, `seasonKill2Sp_Type2()`, `odeIntegSeason()`, `periodAnalysis()`, and `periodFFT()`.

## Other scripts

* `describe_forcing.m`: This script describes the forcing signal depending on &#949; and &#952;, its variance over one period, and its net magnitude when the variance is controlled.
The resulting figures correspond to Fig. 6 in the *JTB* manuscript.

