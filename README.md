# Benchmark Results Overview

This document presents benchmark results for scalar multiplication on three elliptic curves: MNT-992, Lollipop-956-451, The benchmarks compare the performance of the optimized GLV (Gallant–Lambert–Vanstone) method against the standard double-and-add scalar multiplication. For more details, refer to the paper and especially to its Section 3.

## Usage

```shell
cd python
sage -python bench.py
================================================================================
                                    MNT6-992                                    
================================================================================
scalar multiplication by 2^l':                   5.800ms
endomorphismϕ:           3.095ms (87% faster)

standard scalar multiplication:                  19.007ms
GLV based on ϕ:            14.392ms (32% faster)
================================================================================
                                    MNT4-992                                    
================================================================================
scalar multiplication by 2^l':                   5.758ms
endomorphism ϕ:           3.110ms (85% faster)

standard scalar multiplication:                 19.105ms
GLV based on ϕ:           14.620ms (31% faster)
================================================================================
                                lollipop-956451                                 
================================================================================
scalar multiplication by 2^l':                  2.593ms
endomorphism ϕ:           1.407ms (84% faster)

standard scalar multiplication:                 8.782ms
GLV based on ϕ:            6.384ms (38% faster)
```

**N.B.** The benchmarks are averaged over 1000 instances of the routine. Besides, the endomorphism $\phi$ is not evaluated via (the homogeneous version of) Horner's scheme as proposed in the paper, but more elementarily. Therefore, (GLV based on) $\phi$ is even faster when implemented more properly.

## How to Interpret Results

The benchmark results provide valuable insights into the efficiency of scalar multiplication optimizations. Below are guidelines to help interpret the results:

1. **Execution Time**
  * The results are displayed in milliseconds (ms), indicating the time taken for each operation.
  * Lower values signify better performance.
2. **Percentage Improvement**
 * The percentage improvements demonstrate how much more efficient the optimized solutions (the endomorphism $\phi$ and GLV based on $\phi$) are compared to the double(-and-add) scalar multiplication.
3. **Comparative Context**
 * In the first two lines of each curve’s results, the time of evaluating the endomorphism $\phi$ is compared against that taken to double a point $\ell^\prime = \lceil \ell/2 \rceil$ times, where the value $\ell$ is indicated in the third column of Table 1 from the paper.
 * The second entry compares the standard scalar multiplication against the GLV optimization using $\phi$.
