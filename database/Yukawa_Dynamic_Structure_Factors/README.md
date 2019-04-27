
# Yukawa Dynamic Structure Factors

### What?

We have computed the dynamic structure factor (DSF) $S(k,\omega)$ for the Yukawa model for a wide range of $\kappa$ and $\Gamma$: 

* kappa = 0; Gamma = [10, 50, 150]
* kappa = 1; Gamma = [14, 72, 217]
* kappa = 2; Gamma = [31, 158, 476]
* kappa = 3; Gamma = [100, 503, 1510]

These values were chosen to roughly follow contours of constant effective coupling. 


The data contains the DSF for only $0.0 < \omega < 3$. If you need higher frequency data, please contact Prof. Murillo murillom@msu.edu. 

* The data is written in Python Numpy format.
* The data size is (166x763) and covers $0.18 < q < 30$, and $0 < \omega < 3$.
* The 0th row is the frequency
* The 0th column is the wavenumber

### How? 

The Sarkas MD code was used.

Yukawa reduced unit is used (time is normalized by the ion plasma frequency, length is normalized by the ion sphere
radius)

Simulation parameters:
* $dt = 0.01$
* Number of timesteps = 80000
* Number of particles = 10000


### Citation

Yongjun Choi, Gautham Dharuman, and Michael S. Murillo, High-Frequency Response of Classical Strongly Coupled Plasmas,
Physical Review E, currently under review.
