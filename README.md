# Non-reciprocal interactions on a membrane
This project simulates the 1D dynamics of a membrane (height $h$) coupled to two non-reciprocally interacting enyzmes with concentrations $\psi_1$ and $\psi_2$. We use a pseudo-spectral method to solve the equations
```math
\partial_t\psi_1 = \partial_x^2\bigg[a_1\psi_1 + (\chi + \alpha)\psi_2 - K\partial_x^2\psi_1 + \kappa s_1\partial_x^2 h\bigg] + N_{\psi_1},
```
```math
\partial_t\psi_2 = \partial_x^2\bigg[a_2\psi_2 + (\chi - \alpha)\psi_1 - K\partial_x^2\psi_2 + \kappa s_2\partial_x^2 h\bigg] + N_{\psi_2},
```
```math
\zeta\partial_t h = -\kappa\partial^4_xh - \kappa s_1\partial_x^2\psi_1 - \kappa s_2\partial_x^2\psi_2 + \sigma\partial_x^2h + N_h,
```
where $N_{\psi_1}$, $N_{\psi_2}$, and $N_h$ are nonlinear terms that will be detailed in the upcoming paper. The parameters $s_1$ and $s_2$ are called spontaneous curvatutures and represent the local curvature induced by an enzyme being attached to the membrane. The linear part of the coupled eqautions is eveloved with a [Crank-Nicolson](https://en.wikipedia.org/wiki/Crank%E2%80%93Nicolson_method) scheme, and the non-linear part is treated with a 2nd order [Adams-Bashforth](https://en.wikipedia.org/wiki/Linear_multistep_method) scheme. 

## How to use
- Input parameters in ```params.txt``` and run the simulation using ```bash run_which.sh -script pbc_main.py```. 
- Plotting can be done with ```bash run_which.sh -script pbc_mov.py```.
- Various options and parameters can be changed in params.txt and new options can be added by updating the ```get_args()``` function in ```pbc_utils.py```

## Dependencies and Packages
Outside of the packages in the standard python library, you will need these:
- [ffmpeg](https://ffmpeg.org/) for turning a series of plots into a movie
- [Matplotlib](https://matplotlib.org/) is a comprehensive library for creating static, animated, and interactive visualizations in Python.
- [NumPy](https://numpy.org/) is the standard package for scientific computing with Python
- [SciPy](https://scipy.org/) provides higher level scientific computing - this project uses the [fast fourier transform](https://docs.scipy.org/doc/scipy/tutorial/fft.html)
