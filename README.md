# FMM
An implementation of the Fast Multipole Method in C, as decribed in [this](https://math.berkeley.edu/~strain/128b.S20/fmm1.pdf) paper by [John Strain](https://math.berkeley.edu/~strain/).
# Using FMM
Create a `CauchyMultiplier` struct using 
```C
newCauchyMultiplier(double* sources, double* targets, double* input, int n, int precision)
```
 `double* fastMultiply(CauchyMultiplier cm)` uses the FMM to apply the Cauchy matrix determined by the sources and targets to the input with precision p. When finished with the calculations, free the CauchyMultiplier struct with `freeCauchyMultiplier(CauchyMultiplier cm)`.
# Testing
Run `make test` to generate error, speed, and flops data which will be outputed to `output` in `.csv` format.
# Performance
This implementation uses max(log2(n), 2) levels. Also included are some graphs of number of levels being log2(n), just because it is interesting. For matrices of size about 200 and larger this implementation of the fast multipole method uses fewer floating point operations than the naive method (direct multiplication). Interestingly, as far as speed when tested on my computer goes, the degree to which increasing the precision increases the runtime of the program is smaller than expected, while the numbers of flops more properly reflects this increase.
## Number of Levels = max(log2(n), 2)
![errorplot](https://github.com/kavalee/FMM/raw/master/images/error-6.png)
![flopsplot](https://github.com/kavalee/FMM/raw/master/images/flops-6.png)
![speedplot](https://github.com/kavalee/FMM/raw/master/images/speed-6.png)
## Number of Levels = log2(n)
![errorplot](https://github.com/kavalee/FMM/raw/master/images/error.png)
![flopsplot](https://github.com/kavalee/FMM/raw/master/images/flops.png)
![speedplot](https://github.com/kavalee/FMM/raw/master/images/speed.png)
