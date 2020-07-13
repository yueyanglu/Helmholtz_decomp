# Helmholtz_decomp

Code of Helmholtz decomposition of a vector field, i.e., computing the streamfuction and potential by solving a minimization problem proposed by [Li et al. (2006)](https://journals.ametsoc.org/doi/full/10.1175/MWR3249.1).

See the [Example](https://github.com/iuryt/vector_fields/blob/master/Example.ipynb) to understand the usage.

## Adaptations from the original
* uv2psiphi() function in which waives the input of the first guess and resolution fields
* periodify() function is used in uv2psiphi() mirroring the data to avoid boundary problems
* General adaptations for Python 3
