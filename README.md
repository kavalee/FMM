# FMM
An implementation of the Fast Multipole Method in C, as decribed in [this](https://math.berkeley.edu/~strain/128b.S20/fmm1.pdf) paper by [John Strain](https://math.berkeley.edu/~strain/).
# Using FMM
Create a `CauchyMultiplier` struct using 
```C
newCauchyMultiplier(double* sources, double* targets, double* input, int n, int precision)
```
 `double* fastMultiply(CauchyMultiplier cm)` uses the FMM to apply the Cauchy matrix determined by the sources and targets to the input with precision p. When finished with the caluculations free the CauchyMultiplier struct with `freeCauchyMultiplier(CauchyMultiplier cm)`.
# Testing
Run `make test` to generate error, speed, and flops data which will be outputed to `output` in `.csv` format.
# Performance
![errorplot](https://github.com/kavalee/FMM/tree/master/images/error.png)


