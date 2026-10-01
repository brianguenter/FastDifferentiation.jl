# Comprehensive Technical Review of Factorization Changes

This document provides a detailed, file-by-file and function-by-function technical review of all modifications on the `fix-complex-expression-factorization` branch compared to `main` in `FastDifferentiation.jl`.  It describes the code as it will be merged, including the follow-up changes from schrimpf/FastDifferentiation.jl#1.

---

## 1. Algorithmic Background & Mathematical Notation

### 1.1 The Derivative Graph and Sum of Path Products
Let $f: \mathbb{R}^n \to \mathbb{R}^m$ be a differentiable function defined by an expression directed acyclic graph (DAG). The derivative graph $\mathcal{G} = (\mathcal{V}, \mathcal{E})$ is a DAG where:
- Each vertex $v \in \mathcal{V}$ corresponds to an intermediate subexpression or variable.
- Each directed edge $e = (u, v) \in \mathcal{E}$ represents the local partial derivative $\frac{\partial v_u}{\partial v_v}$ of parent node $u$ with respect to child argument $v$. Nodes have no operational function; they serve only to connect edges.
- Leaf nodes correspond to input variables $x_1, \dots, x_n$ (domain dimension $n$).
- Root nodes correspond to output components $f_1, \dots, f_m$ (codomain dimension $m$).

The function $f$ has $n \times m$ scalar constituent derivative functions:
$$f_{ij} = \frac{\partial f_i}{\partial x_j}, \quad 1 \le i \le m, \; 1 \le j \le n$$
By the generalized chain rule, $f_{ij}$ equals the **sum of path products** over all directed paths $\Pi(i, j)$ from root $i$ to leaf $j$:
$$f_{ij} = \sum_{\pi \in \Pi(i, j)} \prod_{e \in \pi} \text{value}(e)$$

### 1.2 Factor Subgraphs: Dominance and Postdominance
To avoid exponential path enumeration, the $D^*$ algorithm identifies factorable subgraphs $[b, c]$ defined by dominance relationships:
- **Dominance ($b \text{ dom } c$)**: Node $b$ is on every path from node $c$ to every root. If $b$ has more than one child, $b$ is a **factor node**; if $c$ has more than one parent, $c$ is a **factor base**. The factor subgraph $[b, c]$ consists of $b$, $c$, and all intermediate nodes on any path between them.
- **Postdominance ($b \text{ pdom } c$)**: Node $b$ is on every path from node $c$ to every leaf. If $b$ has more than one parent, $b$ is a **factor node**; if $c$ has more than one child, $c$ is a **factor base**.

### 1.3 Factoring and Subgraph Edges
Factoring evaluates the sum of path products within $[b, c]$ and replaces the subgraph with a new **subgraph edge** $S$ connecting $b$ and $c$:
$$\text{value}(S) = \sum_{\pi \in \Pi(b, c)} \prod_{e \in \pi} \text{value}(e)$$
To preserve the overall sum of path products for paths that do not pass through both $b$ and $c$, edges $e \in [b, c]$ are selectively deleted or retained:
- **`dominatorTest` ($b \text{ dom } c$)**: An edge $e$ is deleted if and only if $c \text{ pdom } e.1$ (all paths from the bottom vertex of $e$ to leaves pass through $c$).
- **`postDominatorTest` ($b \text{ pdom } c$)**: An edge $e$ is deleted if and only if $c \text{ dom } e.2$ (all paths from the top vertex of $e$ to roots pass through $c$).

### 1.4 The Multi-Output Challenge ($f: \mathbb{R}^n \to \mathbb{R}^m$)
In `FastDifferentiation.jl`, edges in $\mathcal{G}$ are shared across constituent functions $f_{ij}$. To track reachability without duplicating edges, each `PathEdge` carries two bitmasks:
- `reachable_variables::BitVector` of length $n$: indicates which variables $\{x_1, \dots, x_n\}$ are reachable below the edge.
- `reachable_roots::BitVector` of length $m$: indicates which roots $\{f_1, \dots, f_m\}$ are reachable above the edge.

In a dominator subgraph, the *dominance mask* is `reachable_roots`, while the *non-dominance mask* is `reachable_variables`. In a postdominator subgraph, the roles are reversed.

---

## 2. Failure Modes Identified on `main`

Prior to this branch, factorization failed on expressions with complex shared structures (e.g. issues [#23](https://github.com/brianguenter/FastDifferentiation.jl/issues/23), [#65](https://github.com/brianguenter/FastDifferentiation.jl/issues/65), [#108](https://github.com/brianguenter/FastDifferentiation.jl/issues/108), and ADNLPModel `:helical`):

1. **Path Reachability Over-Constraining (Issue #23)**:
   In `old_edge_path`, candidate forward edges were filtered using `subset(reachable(subgraph), reachable_variables(x))`. Since `reachable(subgraph)` is the union of *all* variables reachable by the entire subgraph, single-variable paths (where `reachable_variables(x)` is a strict subset of the union) were rejected. This prematurely terminated path traversal, dropping valid derivative terms to zero (`0.0`).

2. **Reachability Collapsing & Mask Aliasing (Issues #65 & #108)**:
   In `evaluate_subgraph`, all paths through $[b, c]$ were summed into a single `Node` and assigned a single replacement `PathEdge` with the union mask. When paths through $[b, c]$ reached disjoint subsets of variables (or roots), collapsing them assigned derivative contributions to variables that had no topological connection to those paths. Furthermore, sharing mutable mask references caused in-place mask resets to corrupt sibling edges, resulting in:
   `AssertionError: Should only be one path from root R to variable V. Instead have 2 children from node N on the path`.

3. **Lossy Scalar Branching Fallback**:
   `is_branching` incorrectly treated paths for disjoint variables sharing a trunk edge as an internal branching diamond. Subgraphs flagged by `is_branching` were passed to `evaluate_branching_subgraph`, which discarded reachability bitmasks entirely and summed derivatives across all variables into an unpartitioned scalar.

4. **Fragile Internal Path Walking in Edge Reset**:
   `add_non_dom_edges!` and `reset_edge_masks!` relied on `edge_path` to walk paths from each forward edge to the dominator node. If internal branching occurred, `edge_path` encountered multiple candidate edges and threw `AssertionError: count ≤ 1`, crashing the factorization.

5. **Non-Idempotent Graph Reconstruction**:
   `DerivativeGraph` unconditionally wrapped each root in a `NoOp` node. When reconstructing a graph from existing roots (e.g. in `_symbolic_jacobian`), this created nested `NoOp(NoOp(...))` wrappers, altering postorder numbers and perturbing heap tie-breaking order in `compute_factorable_subgraphs`.

---

## 3. Function-by-Function Breakdown of Changes

The merged code differs in places from the first version of this branch.  Follow-up changes made factorization scale linearly in the number of inputs and outputs again, removed reachability bits that factoring had left behind, found factorizations that the first version missed, restored products shared between derivatives, and removed the code that the new factorization no longer calls.  This section describes the code as it is merged.  Functions that no longer exist are named only where they explain a design choice.

### 3.1 `src/DerivativeGraph.jl`

#### `DerivativeGraph(roots::AbstractVector, index_type::Type=Int64)`
```julia
# Old:
new_roots[i] = create_NoOp(root)

# New:
new_roots[i] = create_NoOp(is_NoOp(root) ? children(root)[1] : root)
```
- **Rationale**: Unwraps any existing `NoOp` node before creating a fresh root wrapper. This guarantees idempotency: reconstructing a `DerivativeGraph` from its roots produces identical postorder numbering and topological structure.

#### `initialize_edge_masks!`
- **New Function**:
  Computes exact reachability bitmasks directly from the graph topology in two linear-time passes:
  1. **Bottom-Up Pass**: Traverses nodes in postorder ($1 \dots \text{num\_nodes}$), propagating variable reachability upward from leaf variables.
  2. **Top-Down Pass**: Traverses nodes in reverse postorder ($\text{num\_nodes} \dots 1$), propagating root reachability downward from output roots.
  3. **Edge Mask Assignment**: Sets `reachable_variables(edge)` to the bottom vertex's variable mask and `reachable_roots(edge)` to the top vertex's root mask.
- **Invariant**: Every edge $e$ begins with accurate bitmasks matching its reachability in the full unfactored graph:
  $$\text{reachable\_variables}(e)_j = \text{true} \iff \exists \text{ path from } \text{bott\_vertex}(e) \text{ to } x_j$$
  $$\text{reachable\_roots}(e)_i = \text{true} \iff \exists \text{ path from } f_i \text{ to } \text{top\_vertex}(e)$$

---

### 3.2 `src/FactorableSubgraph.jl`

#### `add_non_dom_edges!` & `reset_edge_masks!`
- **Key Algorithmic Improvement (Boundary-Cut Factorization)**:
  - **Old Approach**: Attempted to walk the entire path from each forward edge to the dominator node using `edge_path(subgraph, s_edge)`. When subgraphs contained multiple internal branches, `edge_path` asserted on `count > 1`.
  - **New Approach**: Operates directly on boundary edges incident to `dominated_node(subgraph)`:
    ```julia
    for s_edge in forward_edges(subgraph, dominated_node(subgraph))
        if test_edge(subgraph, s_edge)
            edge_mask = reachable_dominance(subgraph, s_edge)
            diff = set_diff(edge_mask, reachable_dominance(subgraph))
            if any(diff)
                # Split non-dominance reachability on the boundary edge
                ...
                @. edge_mask &= !diff
            end
        end
    end
    ```
    Similarly, `reset_edge_masks!` clears dominance bits on `forward_edges(subgraph, dominated_node(subgraph))` without traversing internal paths:
    ```julia
    for fwd_edge in forward_edges(subgraph, dominated_node(subgraph))
        if test_edge(subgraph, fwd_edge)
            mask = non_dominance_mask(subgraph, fwd_edge)
            @. mask = mask & bypass_mask
            if can_delete(fwd_edge)
                push!(edges_to_delete, fwd_edge)
            end
        end
    end
    ```
- **Mathematical Justification**:
  In a factor subgraph $[b, c]$, all paths entering $[b, c]$ pass through $c$ (for dominator subgraphs, paths from leaves pass through $c$ before reaching $b$). Severing or splitting reachability at the boundary of $c$ removes the factored paths for $[b, c]$ while leaving internal edges intact for any bypass paths that enter the subgraph from outside $[b, c]$ at intermediate nodes. This eliminates the need for recursive internal path walking.

#### Path Walking (Removed)
- Earlier versions of this branch changed `next_valid_edge`, `isa_connected_path`, `edges_on_path`, and `PathIterator` so that a walk along a path keeps the reachability mask of its first edge, and they changed `subgraph_edges` to drop edges that cannot reach the dominating node.
- These functions served `subgraph_exists`, which decided whether a subgraph still needed factoring by walking its paths one at a time.  The merged `factor_subgraph!` decides this by counting paths instead (Section 3.3), so nothing calls them any more.  They have been removed, along with `subgraph_exists`, `check_edges`, `edge_path`, and `deconstruct_subgraph`.

---

### 3.3 `src/Factoring.jl`

#### `on_subgraph_path`
- **New Function**: `on_subgraph_path(subgraph, edge)` is true when the edge is live for at least one dominance bit of the subgraph (a root for a dominator subgraph, a variable for a postdominator subgraph) and for at least one of its non-dominance bits.
- **Rationale**: The older test, `test_edge`, also requires the edge's dominance mask to contain the *whole* dominance mask of the subgraph.  An earlier factorization can split an edge by root (or by variable) with `add_non_dom_edges!`, and the pieces can then fail that test although paths through them still belong to the subgraph.  Those paths were invisible, so some subgraphs that needed factoring were left unfactored.  `test_edge` is still used by `add_non_dom_edges!` and `reset_edge_masks!`, which operate on a partition subgraph whose dominance mask is that of a single replacement edge.

#### `subgraph_path_groups`
- **Replaces `_evaluate_subgraph_paths`**, which grouped the paths from the dominated node to the dominating node by their non-dominance mask alone.
- **Change**: Paths are grouped by the *pair* of masks for which every edge of the path is live.  For a path $\pi$ through $[b, c]$, let $D(\pi)$ be the intersection of the subgraph's dominance mask with the dominance masks of the edges of $\pi$, and let $N(\pi)$ be the same intersection of non-dominance masks.  Each group stores the sum of its path products and its number of paths:
  $$G_{D,N} = \sum_{\pi:\ D(\pi) = D,\ N(\pi) = N} \ \prod_{e \in \pi} \text{value}(e)$$
  Tracking both masks is what makes edges split by an earlier factorization visible.
- **Order of Multiplication**: Each product is formed from the top of the subgraph downward, which is the order in which `follow_path` forms products, so that products shared by several derivatives are built once.  The method for dominator subgraphs recurses upward from the dominated node with memoization.  The method for postdominator subgraphs pushes the groups downward through the nodes of the subgraph in topological order.
- **Hashing**: Dictionaries keyed by masks use the masks' chunk vectors (`mask_key`), because `hash` on a `BitVector` falls back to the generic element-by-element method, which is about ten times slower for 800 bits.

#### `index_classes` and `factored_edges`
- **Replace `evaluate_subgraph`**.
- **Partition of Root–Variable Pairs**: Each path group covers a rectangle of root–variable pairs, $D \times N$.  `factored_edges` partitions the dominance indices, and separately the non-dominance indices, according to which groups contain them (`index_classes`).  Every block of the resulting grid is covered by a single set of groups, so it becomes a single replacement edge, whose value is the sum of those groups' values and whose two masks are the block's two classes.  No pair receives a product from a path that does not reach it.
- **Path Count**: `factored_edges` also returns the largest number of paths that any single pair lies on, which `factor_subgraph!` uses to decide whether to factor the subgraph at all.
- **Scaling**: The first version built a signature vector for every variable (or root) of the whole graph, for every subgraph, which made `factor!` quadratic in the number of inputs or outputs.  `index_classes` refines its classes one mask at a time, so its cost depends on the numbers of groups and classes rather than on the number of indices.  The Hessian of the 2000-dimensional Rosenbrock function took 11.9 s with the first version, against 1.84 s on `main` and 1.88 s with the merged code (Section 5).

#### `factor_subgraph!`
- **When to Factor**: A subgraph is factored only if some root–variable pair lies on at least two of its paths.  Otherwise each pair already has at most one path through the subgraph, and factoring would gain nothing.  This test replaces `subgraph_exists`, which walked paths with `next_valid_edge` and was misled by the stale reachability bits described below.
- **Replacement**: For each replacement edge, `factor_subgraph!` builds a partition subgraph whose dominance and non-dominance masks are those of the edge, and it applies `add_non_dom_edges!` and `reset_edge_masks!` to that partition, so that the boundary edges at the dominated node lose the pairs that the replacement edge now carries.  It then deletes the boundary edges left without pairs and adds the replacement edge.  Finally it calls `prune_stale_masks!`.

#### `subgraph_region` and `prune_stale_masks!`
- **New Functions**.
- **Problem**: Factoring clears reachability only on the boundary edges at the dominated node, so the edges inside the factored subgraph kept the bits of the pairs whose paths had just been replaced.  These stale bits left dead-end edges in the graph.  A stale edge next to a live one also looked like a fork, which made `subgraph_exists` reject subgraphs that still needed factoring, so some root–variable pairs kept more than one path.
- **Fix**: `subgraph_region` records the nodes of the subgraph before it is factored.  Afterward, `prune_stale_masks!` visits those nodes, starting next to the dominated node, removes from each edge the bits that no continuing edge carries, and deletes edges left empty.  When the dominance bits of one edge continue along different edges, the edge is split into one edge per group, as `add_non_dom_edges!` does for boundary edges.  No live path is removed, because every bit that survives is continued by an edge that is live for the same dominance bit, whether that edge lies inside the subgraph or outside it.
- **Effect**: After `factor!` on 500 random expression graphs, the first version of this branch left more than one path for some root–variable pair in 21 graphs and left dead-end edges in 345.  The merged code leaves neither in any of them, and its factored graphs have 29% fewer edges.

#### `follow_path` and `sum_all_paths`
- **On `main`**, `follow_path` asserted that exactly one edge continues the path at each step, and it multiplied the edges of the path in decreasing order of the number of root–variable pairs that use them, so that products shared by many derivatives are formed once.
- **First Version of This Branch**: The function was replaced by a memoized sum over all paths, which cannot fail the assertion but multiplies in path order, so the shared products were lost.  On problems for which `main` was already correct, the results took 10–24% more operations.
- **Merged Code**: `follow_path` again walks the single remaining path and sorts its edges by use.  It falls back to the memoized sum, now called `sum_all_paths`, only if it meets a branch or a dead end, which the changes above should prevent.  On every problem in Section 5 for which `main` is correct, the operation counts now equal those of `main`.

#### Removed Code
- `old_edge_path`, `is_branching`, `evaluate_branching_subgraph`, and `compute_internal_idoms`, along with the helpers that only they used, were already not called by the first version of this branch.  They have been removed, as have `evaluate_subgraph`, `make_factored_edge`, `multiply_sequence`, `path_sort_order`, and `PathConstraint` with its methods.
- The failure of `old_edge_path` described in Section 2 explains wrong results on `main`.  The merged code groups paths instead of walking them, so that function has no successor.

---

### 3.4 `src/Jacobian.jl`

#### `hessian`
```julia
# Old:
tmp = DerivativeGraph(expression)
jac = _symbolic_jacobian!(tmp, variable_order)
tmp2 = DerivativeGraph(vec(jac))
return _symbolic_jacobian!(tmp2, variable_order)

# New:
gradient = jacobian([expression], variable_order)
return jacobian(vec(gradient), variable_order)
```
- **Rationale**: Constructs the gradient using `jacobian([expression], variable_order)` and computes the Hessian by differentiating the gradient vector as a single multi-output graph. This ensures maximum subexpression reuse across all elements of the Hessian matrix.

---

### 3.5 Other Unused Code

- Code elsewhere in the package that was already unused on `main` has also been removed.  Examples are `graph_statistics` (which needed the Statistics package, which the package does not load), `ConstrainedPathIterator`, the `@invariant` macro, and the files `FunctionDerivative.jl` and `UnspecifiedFunction.jl`, which the module never included.
- Before anything was removed, the ten registered packages that depend on FastDifferentiation, and Brian Guenter's other public repositories, were searched for every removed name.  The only one in use was `number_of_operations`, which `brianguenter/Benchmarks` calls, so it was kept.

---

### 3.6 `test/` Test Suite Updates

#### `test/runtests.jl`
- **`FD.compute_factorable_subgraphs test order`**:
  Fixed syntax typo `6_1` -> `_6_1`. Updated priority assertion for subgraphs with identical heap priority (`diff = 3`, `times_used = 1`) to accept either tie-break order.
- **`jacobian` Test Item**:
  Corrected analytical derivative reference:
  $$n_5 = n_2 n_4 = (x y)(x y^2) = x^2 y^3 \implies \frac{\partial n_5}{\partial y} = 3 x^2 y^2$$
  The old test compared against `4 * x^2 * y^2`, which was mathematically incorrect.
- **Tests of Removed Functions**: The test items for `isa_connected_path`, `edge_path`, `PathIterator`, `subgraph_edges`, `deconstruct_subgraph`, `path_sort_order`, and `multiply_sequence` were removed along with those functions, as were the empty `relation_edges` item and the `evaluate_subgraph` item, which made no assertions.
- **`subgraph_region` and `factored_edges`**: The `subgraph_edges` item became a test of `subgraph_region`, with node sets taken from the same reference edges, and the `make_factored_edge` item now checks the masks and path counts returned by `factored_edges`.

#### `test/test-complex-expressions.jl`
- Added rigorous test coverage for issues [#23](https://github.com/brianguenter/FastDifferentiation.jl/issues/23), [#65](https://github.com/brianguenter/FastDifferentiation.jl/issues/65), and [#108](https://github.com/brianguenter/FastDifferentiation.jl/issues/108), validating symbolic Jacobians and Hessians against `FiniteDifferences.central_fdm` down to tolerances of $10^{-6}$ and $10^{-12}$.

#### `test/test-differentiation-interface.jl`
- Added integration tests running `DifferentiationInterfaceTest.default_scenarios()` with `AutoFastDifferentiation()`. All 170 scenarios (7,278 assertions) pass 100%.

#### `test/test-optimization-problems.jl`
- Added stress tests across 52 non-linear optimization problems from `OptimizationProblems.ADNLPProblems`. Validates dense and sparse gradients, Hessians, and constraint Jacobians.
- **Comparative Benchmark**: On `main`, this test suite crashes on `:helical` with `AssertionError: Should only be one path from root 1 to variable 3. Instead have 2 children from node 43 on the path`. On this branch, all 52 problems pass 100%.

---

## 4. Correctness Invariants & Complexity Analysis

### 4.1 Chain Rule Invariant Preservation
Let $\mathcal{G}_0$ be the initial derivative graph and $\mathcal{G}_k$ be the graph after $k$ factorizations. For every root $i \in \{1, \dots, m\}$ and variable $j \in \{1, \dots, n\}$:
$$\sum_{\pi \in \Pi_{\mathcal{G}_k}(i, j)} \prod_{e \in \pi} \text{value}(e) = \sum_{\pi \in \Pi_{\mathcal{G}_0}(i, j)} \prod_{e \in \pi} \text{value}(e) = f_{ij}$$
- **Inside $[b, c]$**: `subgraph_path_groups` and `factored_edges` divide the root–variable pairs of the subgraph into blocks, such that all pairs in a block lie on the same set of paths.  The replacement edge $S_r$ for a block therefore carries, for each of its pairs, exactly the sum of path products that the subgraph contributed to that pair.
- **Boundary Cuts**: For each replacement edge, `add_non_dom_edges!` and `reset_edge_masks!` remove the edge's pairs from the boundary edges at the dominated node, so that the factored paths can be reached only through $S_r$.
- **Stale Bits**: `prune_stale_masks!` removes a bit from an interior edge only if no continuing edge carries it, so every path that is still live keeps all of its edges.
- **Bypass Paths**: Paths that enter intermediate nodes of $[b, c]$ from outside keep their bits, because the boundary cuts touch only the edges at the dominated node, and pruning keeps every bit that some continuing edge carries, whether inside the subgraph or outside it.

### 4.2 Computational Complexity
- **Mask Initialization**: $O((n + m) \cdot |\mathcal{V}| + |\mathcal{E}|)$ using linear two-pass topological DAG traversal.
- **Subgraph Region**: `subgraph_region` is a breadth-first search over the edges of the subgraph, with $O(|\mathcal{V}_{[b, c]}| + |\mathcal{E}_{[b, c]}|)$ edge tests.
- **Path Grouping**: `subgraph_path_groups` performs $O(|\mathcal{E}_{[b, c]}| \cdot g)$ mask operations, where $g$ is the number of distinct pairs of masks among the paths of the subgraph.
- **Partition**: `index_classes` performs $O(g \cdot k)$ mask operations, where $k$ is the number of classes it produces, and `factored_edges` then examines each block of the grid once.  No step loops over individual roots or variables, so apart from word-level operations on the masks themselves, the cost does not grow with $n$ or $m$.
- **Pruning**: `prune_stale_masks!` compares each edge of the subgraph with the edges that continue it.

---

## 5. Verification Summary

All results were measured on one laptop.  The full test runs used both Julia 1.12 and 1.13, and everything else used Julia 1.12.

### 5.1 Correctness

| Check | Merged code | `main` |
| :--- | :--- | :--- |
| **All 59 test items**, Julia 1.12 and 1.13 | **Pass** | — |
| **`test-complex-expressions.jl`** (#23, #65, #108) | **Pass** | Fails (zero derivatives and assertions) |
| **`test-optimization-problems.jl`** (52 ADNLPModels) | **Pass** | Fails on `:helical` (multi-path assertion) |
| **`test-differentiation-interface.jl`** (170 scenarios) | **Pass** | Pass |
| **Jacobians of 3,000 small random graphs**, compared with ForwardDiff | **All agree** | 210 wrong, 21 errors |
| **Jacobians of 900 larger random graphs**, compared with ForwardDiff | **All agree** | 374 wrong, 192 errors |
| **Spherical-harmonics finite-difference test** from the test suite, orders 8, 10, and 12 | **Pass** | Fails at all 125 points |
| **Spherical-harmonics Jacobian**, degrees up to 10, compared with ForwardDiff | Largest error $1.2 \times 10^{-16}$ of the largest entry | Largest error 5,500 times the largest entry |
| **Chain Hessian**, $n = 200$, compared with ForwardDiff | Largest error $2.2 \times 10^{-16}$ of the largest entry | Largest error $1.3 \times 10^{-2}$ of the largest entry |
| **Root–variable pairs with more than one path** after `factor!`, 500 random graphs | **None** | — |

The small random graphs have 1–4 inputs, 1–3 outputs, and 4–18 operations; the larger ones have 2–8 inputs, 2–8 outputs, and 20–70 operations.  Their operations are drawn from $a + b$, $a - b$, $a b$, $\sin a$, $\cos a$, $a / (1 + b^2)$, and $a^2$, with operands biased toward recent results so that subexpressions are shared.  The errors on `main` are `AssertionError`s and `KeyError`s.  The chain is $\sum_i \sin(x_i x_{i+1}) \exp(x_i)$, and the spherical harmonics are the unnormalized recurrences from the test suite.

### 5.2 Performance

Times are for the `jacobian` or `hessian` call, as the shortest of two to five runs, and operation counts are from `number_of_operations`.

| Problem | `main` | First version of this branch | Merged code |
| :--- | :--- | :--- | :--- |
| Spherical harmonics, degrees up to 20 | 0.155 s, 8,825 ops (wrong) | 0.175 s, 8,985 ops | 0.078 s, 6,188 ops |
| Rosenbrock Hessian, $n = 200$ | 0.051 s, 1,790 ops | 0.169 s, 1,790 ops | 0.063 s, 1,790 ops |
| Rosenbrock Hessian, $n = 2000$ | 1.84 s, 17,990 ops | 11.9 s, 17,990 ops | 1.88 s, 17,990 ops |
| Chain gradient, $n = 1000$ | 0.299 s, 9,989 ops | 0.751 s, 10,987 ops | 0.294 s, 9,989 ops |
| Chain Hessian, $n = 200$ | 0.059 s, 3,986 ops (wrong) | 0.314 s, 5,656 ops | 0.066 s, 4,975 ops |
| Dense map Jacobian, $n = 40$ | 0.075 s, 6,481 ops | 0.072 s, 8,041 ops | 0.069 s, 6,481 ops |

The dense map has components $x_i \sin\bigl(\sum_j x_j (i + j) / n\bigr)$.

The generated Rosenbrock Hessians have the same operation counts on `main`, on the first version of this branch, and in the merged code, and evaluating them in place takes the same time on all three to within measurement noise: about 0.25 μs for $n = 50$ and 2.5 μs for $n = 200$.  This measurement therefore does not reproduce the slower evaluation reported for `rosenbrock50` in `benchmark/BENCHMARK_RESULTS.md`.
