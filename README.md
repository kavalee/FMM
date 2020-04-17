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
This implementation uses max(log2(n), 2) levels. Also included are some graphs log2(n) levels, just because it is interesting. 
For matrices of size about 200 and larger this implementation of the fast multipole method uses fewer floating point operations than the naive method (direct multiplication). Interestingly, as far as speed when tested on my computer goes, the degree to which increasing the precision increases the runtime of the program is smaller than expected; the numbers of flops more properly reflects this increase. The main purpose of the speed plot is to confirm that the flop plot can be trusted.

This algorithm runs in O(nlogn) time, however it seems to run close to linear time as shown by the plot of flops/n and flops/(nlogn) below (e.g flops/(nlogn) seems to approach a constant while flops/n increases but at a very slow rate). This is most likely due to how in this implementation the source moments are only computed at the finest level, but the Taylor coeffiecients are computed at every level. Since there is O(logn) levels, and the Taylor coeffiecient evaluation is called about 3n times per level, the resulting O(nlogn) behaviour is expected.

As far as programming goes, binomial coefficients are cached from n = 0 to 2 * p, and powers are accumulated instead of using pow(). Both of these save quite a few of flops. The array of source moments, and Taylor coefficients are reallocated on each level which is inefficent, however it probably doesn't have a huge effect of performance. As of now this implementation only works for Cauchy matricies with entries on (0, 1) however making it work for arbitrary intervals would only add a few flops per cell lookup, (if anyone wants this implemented contact me and I can add it).
## Number of Levels = max(log2(n), 2)
![errorplot](https://github.com/kavalee/FMM/raw/master/images/error-6.png)
![flopsplot](https://github.com/kavalee/FMM/raw/master/images/flops-6.png)
![speedplot](https://github.com/kavalee/FMM/raw/master/images/speed-6.png)
![accuracyp](https://github.com/kavalee/FMM/raw/master/images/asymptotic.png)
## Number of Levels = log2(n)
![errorplot](https://github.com/kavalee/FMM/raw/master/images/error.png)
![flopsplot](https://github.com/kavalee/FMM/raw/master/images/flops.png)
![speedplot](https://github.com/kavalee/FMM/raw/master/images/speed.png)

