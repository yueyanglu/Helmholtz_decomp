# Helmholtz_decomp

Helmholtz decomposition of a vector field, i.e., computing the streamfuction and potential by solving a minimization problem proposed by [Li et al. (2006)](https://journals.ametsoc.org/doi/full/10.1175/MWR3249.1), without the need of specifying boundary conditions.

The optimization function is [minFunc](https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html) instead of MATLAB's [fminunc](https://www.mathworks.com/help/optim/ug/fminunc.html) for computational efficiency. The quasi-Newton algorithm with limited-memory BFGS updates is used by default.

The code is inspired by the [Python version](https://github.com/iuryt/vector_fields) but with a more general configuration of grid and an accurate approach to calculating the Jacobian of the Tikhonov’s functional.

This code has been used in the studies of [Kamenkovich et al. (2021)](https://doi.org/10.1029/2020GL091719), [Kamenkovich and Garraffo (2022)](https://doi.org/10.1175/JPO-D-21-0244.1), and [Lu et al. (2022)](https://doi.org/10.1175/JPO-D-22-0108.1).
