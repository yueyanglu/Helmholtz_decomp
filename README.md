# Helmholtz_decomp

Code of Helmholtz decomposition of a vector field, i.e., computing the streamfuction and potential by solving a minimization problem proposed by [Li et al. (2006)](https://journals.ametsoc.org/doi/full/10.1175/MWR3249.1), without specifying boundary conditions.

The optimization function is [minFunc](https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html) instead of MATLAB's [fminunc](https://www.mathworks.com/help/optim/ug/fminunc.html) for computational efficiency. The quasi-Newton algorithm with limited-memory BFGS updates is used by default.

The code is inspired by the [Python version](https://github.com/iuryt/vector_fields) with a more general configuration of grid and a different way of dealing with the objective functional and Jacobian.
