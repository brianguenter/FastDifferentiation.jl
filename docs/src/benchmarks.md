# Comparison of FD with other AD algorithms
See [Benchmarks.jl](https://github.com/brianguenter/Benchmarks) for the benchmark code used to generate this table.
| Function | FD sparse | FD dense | ForwardDiff | ReverseDiff | Enzyme | Zygote |
|---------|-----------|----------|-------------|-------------|--------|--------|
| Rosenbrock Hessian | **1.00** | 75.60 | 571669.52 | 423058.61 | [^notes] | 1015635.96 |
| Rosenbrock gradient | [^notes] | 1.28 | 682.41 | 306.27 | **1.00** | 4726.62 |
| Simple matrix Jacobian | [^notes] | **1.00** | 42.61 | 54.60 | [^notes] | 130.13 |
| Spherical harmonics Jacobian | [^notes] | **1.00** | 36.00 | [^notes] | [^notes] | [^notes] |


 ## Comparison of AD algorithms with a hand optimized Jacobian
| FD sparse | FD Dense | ForwardDiff | ReverseDiff | Enzyme | Zygote | Hand optimized|
|-----------|----------|-------------|-------------|--------|--------|---------------|
 **1.00** | 1.81 | 29.45 | [^notes] | [^notes] | 556889.67 | 2.47 |


It is worth nothing that both FD sparse and FD dense are faster than the hand optimized Jacobian.

It is also intersting to note the ratio of the number of operations of the Jacobian or Hessian of a function to the number of operations in the original function. For FD these ratios are:

|Ratios | Rosenbrock Jacobian | Rosenbrock Hessian | Spherical harmonics Jacobian | Simple matrix function | hand optimized function |
|-------|---------------------|--------------------|------------------------------|------------------------|-------------------------|
|       | 1.13                | 1.13               | 2.4                          |          4.0           |     .59                 |

Contrary to expectation in most of these benchmarks the computation count of the Jacobian or Hessian is a small constant times larger than the computation count of the original function. This is a very small sample of functions but it will be interesting to see if this generalizes to all functions or just functions with special graph structure.

[^notes]: For the FD sparse column, FD sparse was slower than FD dense so times are not listed for this column. For all other columns either the benchmark code crashes or I haven't yet figured out how to make it work correctly.
