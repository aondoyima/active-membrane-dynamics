# Non-reciprocal interactions on a membrane
This project simulates the 1D dynamics of a membrane (height $h$) coupled to two non-reciprocally interacting enyzmes with concentrations $\psi_1$ and $\psi_2$.

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
