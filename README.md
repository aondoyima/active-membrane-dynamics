# Non-reciprocal interactions on a membrane
This project simulates the 1D dynamics of a membrane (height $h$) coupled to two non-reciprocally interacting enyzmes with concentrations $\psi_1$ and $\psi_2$. We use a pseudo-spectral method to solve the equations
```math
\partial_t\psi_1 = \partial_x^2\bigg[\tilde{a}_1\psi_1 + (\tilde{\chi} + \alpha)\psi_2 - K\partial_x^2\psi_1 + \kappa s_1\partial_x^2 h\bigg] + N_{\psi_1},
```
where $N_{\psi_1}$, $N_{\psi_2}$, and $N_h$ are nonlinear terms that will be detailed in the upcoming paper. The parameters $s_1$ and $s_2$ are called spontaneous curvatutures and represent the local curvature induced by an enzyme being attached to the membrane. 

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
