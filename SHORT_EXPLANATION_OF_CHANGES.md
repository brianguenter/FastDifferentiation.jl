# Summary of Changes: D* Derivative Graph Factorization for Multi-Output Functions

This pull request resolves long-standing issues in the $D^*$ derivative graph factorization algorithm for functions $f: \mathbb{R}^n \to \mathbb{R}^m$ with densely interconnected expression DAGs (specifically GitHub issues [#23](https://github.com/brianguenter/FastDifferentiation.jl/issues/23), [#65](https://github.com/brianguenter/FastDifferentiation.jl/issues/65), and [#108](https://github.com/brianguenter/FastDifferentiation.jl/issues/108)).

---

## 1. Algorithmic Context & Problem Formulation

In $D^*$, the derivative of a multi-output function $f: \mathbb{R}^n \to \mathbb{R}^m$ is represented by a directed acyclic derivative graph $\mathcal{G} = (\mathcal{V}, \mathcal{E})$. Each directed edge $e = (u, v)$ represents the local partial derivative $\frac{\partial v_u}{\partial v_v}$. Each of the $n \times m$ constituent functions
$$f_{ij} = \frac{\partial f_i}{\partial x_j}$$
is evaluated as the **sum of path products** over all directed paths from root/output node $i$ ($1 \le i \le m$) to leaf/variable node $j$ ($1 \le j \le n$):
$$f_{ij} = \sum_{\pi \in \Pi(i, j)} \prod_{e \in \pi} \text{value}(e)$$

Factorization reduces the number of multiplications and additions by identifying **factor subgraphs** $[b, c]$:
- In a **dominator subgraph** ($b \text{ dom } c$), factor node $b$ dominates factor base $c$ (every path from $c$ to root $i$ passes through $b$).
- In a **postdominator subgraph** ($b \text{ pdom } c$), factor node $b$ postdominates factor base $c$ (every path from $c$ to variable leaf $j$ passes through $b$).

Evaluating $[b, c]$ yields an equivalent factored subgraph edge $S$ from $b$ to $c$. In $f: \mathbb{R}^1 \to \mathbb{R}^1$, all paths in $[b, c]$ belong to the single derivative function. However, in $f: \mathbb{R}^n \to \mathbb{R}^m$ (Section 4.1 of the paper), edges inside $[b, c]$ are shared across distinct constituent functions $f_{ij}$.

---

## 2. Root Causes of Failures in `main`

Prior to this branch, factorization on expressions with shared multi-output/multi-variable structure suffered from four interrelated failure modes:

1. **Path Reachability Mismatch & Erasure (Issue #23)**:
   In `old_edge_path`, candidate forward edges were filtered using `subset(reachable_mask, reachable_variables(x))` (or roots), where `reachable_mask` was the union of *all* reachable variables across the entire subgraph. If a path reached only a single variable $x_j$, `subset` evaluated to `false`, causing path discovery to fail and dropping valid derivative paths to zero (`0.0`).
2. **Lossy Reachability Collapsing (Issues #65 & #108)**:
   When evaluating factor subgraph $[b, c]$, `evaluate_subgraph` summed all paths into a single scalar expression and created a single replacement edge $S$ with the union mask. When paths within $[b, c]$ reached different subsets of variables or roots, this cross-contaminated reachability, assigning derivative path products to wrong $f_{ij}$ components and causing subsequent downstream factorizations to assert:
   `AssertionError: Should only be one path from root R to variable V`.
3. **Flawed Branching Fallback**:
   `is_branching` incorrectly flagged subgraphs as "branching" whenever two paths crossed a shared trunk edge, even when those paths belonged to disjoint variables. It routed evaluation to `evaluate_branching_subgraph`, which discarded reachability masks entirely.
4. **Path Traversal Crossover & Non-Dominance Splitting**:
   `next_valid_edge` did not enforce reachability mask continuity along traversed paths. Additionally, walking interior paths via `edge_path` during non-dominance splitting failed when internal branching existed.

---

## 3. Summary of Key Fixes

- **Path Groups Keyed by Both Masks** ([`src/Factoring.jl`](src/Factoring.jl)):
  `subgraph_path_groups` uses dynamic programming to group the paths through $[b, c]$ by the pair of masks, dominance and non-dominance, for which every edge of the path is live.  Tracking both masks keeps visible the paths through edges that an earlier factorization split (`on_subgraph_path`).
- **Topological Edge Mask Initialization** ([`src/DerivativeGraph.jl`](src/DerivativeGraph.jl)):
  Added `initialize_edge_masks!` to rigorously propagate variable reachability upward from leaves and root reachability downward from roots, ensuring consistent initial bitmasks for all derivative edges.
- **Pruning Stale Reachability After Factoring** ([`src/Factoring.jl`](src/Factoring.jl)):
  `prune_stale_masks!` removes the bits that factoring leaves on the edges inside $[b, c]$ and deletes edges left empty, so that no root–variable pair keeps more than one path and no dead-end edges remain.
- **Reachability-Partitioned Subgraph Factoring** ([`src/Factoring.jl`](src/Factoring.jl)):
  `factored_edges` partitions the root–variable pairs of $[b, c]$ into blocks whose pairs lie on the same set of paths, and it replaces $[b, c]$ with one edge per block, preserving distinct $f_{ij}$ reachabilities.  The partition visits only the indices that occur on some path, so factorization scales linearly in the number of inputs and outputs.  A subgraph is factored only if some pair lies on at least two of its paths.
- **Boundary-Cut Non-Dominance Edge Splitting & Reset** ([`src/FactorableSubgraph.jl`](src/FactorableSubgraph.jl)):
  In `add_non_dom_edges!` and `reset_edge_masks!`, operations are performed directly on the boundary edges incident to `dominated_node(subgraph)` (`forward_edges(subgraph, dominated_node(subgraph))`). This cleanly eliminates the fragile internal path-walking (`edge_path`) while leaving interior bypass paths entering from outside $[b, c]$ intact.
- **Robust Path Evaluation for Factored Graphs** ([`src/Factoring.jl`](src/Factoring.jl)):
  `follow_path` walks the single remaining path from root $i$ to variable $j$ and multiplies its edges in decreasing order of use, so that products shared between derivatives are formed once.  If it meets a residual branch, `sum_all_paths` sums all valid path products without triggering assertions.
- **Idempotent DerivativeGraph Construction** ([`src/DerivativeGraph.jl`](src/DerivativeGraph.jl)):
  `create_NoOp` now unwraps pre-existing `NoOp` wrappers (`is_NoOp(root) ? children(root)[1] : root`), ensuring graph rebuilding preserves exact postorder node numbering.
- **Multi-Output Hessian Construction** ([`src/Jacobian.jl`](src/Jacobian.jl)):
  Updated `hessian` to differentiate the gradient vector through `jacobian(vec(gradient), variable_order)`, enabling full subexpression reuse across the entire Hessian.
- **Removal of Unused Code**:
  The functions that the new factorization no longer calls, such as the path-walking functions, `subgraph_exists`, `is_branching`, and `evaluate_branching_subgraph`, have been removed, as has code that was already unused on `main`.  Sections 3.2, 3.3, and 3.5 of [`EXPLANATION_OF_CHANGES.md`](EXPLANATION_OF_CHANGES.md) list them.

---

## 4. Test Suite & Verification

All 59 test items pass on Julia 1.12 and 1.13:
- **Existing Suite** ([`test/runtests.jl`](test/runtests.jl)): Fixed test typos and corrected analytical derivative reference in `jacobian` test item ($\frac{\partial}{\partial y}(x^2 y^3) = 3 x^2 y^2$ instead of $4 x^2 y^2$).  The tests of removed functions were removed, and two were rewritten to test `subgraph_region` and `factored_edges`.
- **Complex Expressions Regression Suite** ([`test/test-complex-expressions.jl`](test/test-complex-expressions.jl)): Added rigorous numerical validation for issues [#23](https://github.com/brianguenter/FastDifferentiation.jl/issues/23), [#65](https://github.com/brianguenter/FastDifferentiation.jl/issues/65), and [#108](https://github.com/brianguenter/FastDifferentiation.jl/issues/108) against finite differences.
- **DifferentiationInterface Standard Scenarios** ([`test/test-differentiation-interface.jl`](test/test-differentiation-interface.jl)): All 170 scenarios of `DifferentiationInterfaceTest.default_scenarios()` pass (7,278 assertions).
- **OptimizationProblems Test Suite** ([`test/test-optimization-problems.jl`](test/test-optimization-problems.jl)): 52 unconstrained and constrained optimization problems from `OptimizationProblems.ADNLPProblems` pass 100%. (On `main`, this suite fails on `:helical` with a multi-path assertion error).
- **Random Expression Graphs**: The Jacobians of 3,900 random expression graphs agree with ForwardDiff.  On the same graphs, `main` gives 584 wrong Jacobians and raises 213 errors.  Section 5 of [`EXPLANATION_OF_CHANGES.md`](EXPLANATION_OF_CHANGES.md) has the details, along with timings and operation counts.
