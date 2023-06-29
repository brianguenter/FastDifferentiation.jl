# Benchmarks
See [Benchmarks.jl](https://github.com/brianguenter/Benchmarks) for the benchmark code used to generate these results.

The benchmarks test the speed of gradients, Jacobians, Hessians, and the ability to exploit sparsity in the derivative. The last problem, `ODE`, also compares the AD algorithms to a hand optimized Jacobian. Take these with a grain of salt; they may be useful for order of magnitude comparisons but not much more.

I believe the benchmark code reflects the best way to use each package. However, I am not an expert in any of these packages. For some of the benchmarks I have not yet figured out how to correctly and efficiently compute all the derivatives.

A notable case is Zygote which has unusually slow timings. It is possible it is not being used as efficiently as possible. 

If you are expert in any of these packages please submit a PR to fill in, improve, or correct a benchmark.

When determining which AD algorithm to use keep in mind the limitations of **FD**: operation count and conditionals. The total operation count of your expression should be less than 10⁵. You may get reasonable performance for expressions as large as 10⁶ operations but expect very long compile times. FD does not support conditionals which involve the differentiation variables (yet). The other algorithms do not have these limitations.

These timings are just for evaluating the derivative function. They do not include preprocessing time to generate either the function or auxiliary data structures that make the evaluation more efficient.

The times in each row are normalized to the shortest time in that row. The fastest algorithm will have a relative time of 1.0 and all other algorithms will have a time ≥ 1.0. Smaller numbers are better.


| Function | FD sparse | FD dense | ForwardDiff | ReverseDiff | Enzyme | Zygote |
|---------|-----------|----------|-------------|-------------|--------|--------|
| Rosenbrock Hessian | **1.00** | 75.60 | 571669.52 | 423058.61 | [^notes] | 1015635.96 |
| Rosenbrock gradient | [^notes] | 1.28 | 682.41 | 306.27 | **1.00** | 4726.62 |
| Simple matrix Jacobian | [^notes] | **1.00** | 42.61 | 54.60 | [^notes] | 130.13 |
| Spherical harmonics Jacobian | [^notes] | **1.00** | 36.00 | [^notes] | [^notes] | [^notes] |


 ## Comparison of AD algorithms with a hand optimized Jacobian for an ODE problem
| FD sparse | FD Dense | ForwardDiff | ReverseDiff | Enzyme | Zygote | Hand optimized|
|-----------|----------|-------------|-------------|--------|--------|---------------|
 **1.00** | 1.81 | 29.45 | [^notes] | [^notes] | 556889.67 | 2.47 |


It is worth nothing that both FD sparse and FD dense are faster than the hand optimized Jacobian.


It is also intersting to note the ratio of the number of operations of the FD Jacobian of a function to the number of operations in the original function. 

Problem sizes in approximately the ratio 1:10:100 were computed for several of the benchmarks. The parameters which give these ratios were: ((10,4,2),(100,11,4),(1000,35,9)) for (Rosenbrock Jacobian, Spherical harmonics Jacobian, Simple matrix ops Jacobian), respectively. 

The ratio (jacobian operations)/(original function operations) stays close to a constant over 2 orders of magnitude of problem size. In the case of the simple matrix ops Jacobian both the input and output dimensions scale from 4 to 81, a 20x increase. But the Jacobian computation relative to the original function computation only scales by 1.5x.

|Relative problem size | Rosenbrock Jacobian | Spherical harmonics Jacobian | Simple matrix ops Jacobian |
|-------|---------------------|------------------------------|------------------------|
|  1x     | 1.13                | 2.2                          |          2.6           |
|  10x     | 1.13                | 2.34                          |          3.5          |
|  100x     | 1.13                | 2.4                          |          3.8          |

This is a very small sample of functions but it will be interesting to see if this generalizes to all functions or only applies to functions with special graph structure.

[^notes]: For the FD sparse column, FD sparse was slower than FD dense so times are not listed for this column. For all other columns either the benchmark code crashes or I haven't yet figured out how to make it work correctly.
