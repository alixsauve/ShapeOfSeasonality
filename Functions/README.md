## Functions

These functions are called by the MATLAB scripts shared in `../Mains`.

# Species growth rates

* `seasonKill2Sp_Type1.m`: This function describes the population growth rates with a generalised Lotka-Volterra model of predation. The predator functional response is thus of type I.

* `seasonKill2Sp_Type2.m`: This function describes the population growth rates with a generalised Lotka-Volterra model of predation with a type II functional response for the predator.

Both functions integrate a seasonal forcing exerted on the predator attack rate $$a$$ with magnitude controlled by parameter $$\epsilon$$ and shape controlled by parameter $$\theta$$. All parameters are defined as global variables. The functions take three input arguments: `xt` which is a vector of species densities at time $$t$$, `tt` which is a scalar for time $$t$$, and `scale` which is a logical (`TRUE` if the signal variance is controlled, and `FALSE` otherwise).

# Period detection

* `periodAnalysis2.m`: This function estimates the system periodicity following [Taylor et al.'s (2013)](https://doi.org/10.1098/rspb.2012.2714) algorithm. Non-integer periods and periods above 9 are classified similarly (ouput value being `nan`).

* `periodFFT.m`: This function returns the system periodicity based on the discrete fast Fourier transform (using `fft()` of MATLAB). These outputs are not used in the *JTB* manuscript.

# Numerical integration of discontinuous functions

* `odeIntegSeason.m`: This function runs MATLAB `ode15s()`, but interrupts the numerical integration at each discontinuity of the ODE (i.e., when changing season) to avoid errors from the ODE integrator. The functions takes initial conditions (`Y0`), the number of years simulated (`TFinal`), the time range to record (`TSave`), the function to integrate (`FuncGrowth`), and numerical integration options (`Options`). Just like any MATLAB ODE integrator, it returns the state variables at each recorded time step.
Note that this function is specifically designed to simulated dynamics between 0 and `TFinal`, recording the state variables every $$\nicefrac{1}{52}$$ time units (i.e., every week), with interruption every season change, i.e., when *t = Y* or *t = Y + 0.5* (*Y* being the year of the simulation).
