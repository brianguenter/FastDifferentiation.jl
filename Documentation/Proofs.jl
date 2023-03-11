"""
Given dom subgraph (top,bott) & bott' ∈ (top,bott) & top' = idom(bott'). Then dom subgraph (top',bott') ⊂ (top,bott).
This is used to find factorable subgraphs contained inside other factorable subgraphs. These interior subgraphs must be factored first since local methods to determine which edges to keep or prune depend on all interior subgraphs being factored.

Proof: assume that top' occurs higher in the graph than top, i.e., there is a path from bott'->top and then from top->top' (all paths from bott' must run through top because dom(top,bott')). If top' is somewhere on the upward path from bott' before top' then top' ∈ (top,bott). If top' is on the path from top to top' then top' cannot be the idom of bott, because dom(top,bott') and top' comes after top in the upward recursion toward the root. This violates the assumption top' = idom(bott') so top' must be in (top,bott).



Assume there is a node top' ∈ (top,bott) and a node bott' = pidom(top') with bott' ∈ (top,bott). Then pdom subgraph (bott',top') ⊂ (top,bott)


Definion of a factorable pdom subgraph: a pdom subgraph is a two nodes (bott,top) s.t. bott < top & pdom(bott,top) & num_children(top) > 1 & num_parents(bott) > 1 

"""
test = true

"""
Contained dom or pdom subgraphs must be factored before their containing subgraphs. The following properties determine when containment occurs. 

Definition of factorable dom subgraph (top,bott): a dom subgraph is two nodes top,bott s.t. top > bott & dom(top,bott) & num_children(top) > 1 & num_parents(bott) > 1.

Definition: a node n is inside dom subgraph (top,bott), n ∈ (top,bott), if n is on the upward path from bott to top. 

If n ∈ (top,bott) then dom(top,n). Proof: Trivial from definition of dom subgraph. All paths upward from bott must pass through top, so all paths from nodes on the path bott->top must also pass through top. Otherwise there would be a path from bott to root that did not go through top but this would violate the assumption dom(top,bott). Since all paths from n must pass through top then dom(top,n)

Given dom subgraph (top,bott) and node bott' inside (top,bott) (easily determined by recursing upward through the parents of bott). If there is a node top' s.t. top' = idom(bott') then top' is in (top,bott) and dom(top,top').

Given a factorable dom subgraph (top,bott) and nodes bott' with bott' ∈ (top,bott) s.t. top' = idom(bott'). If  and bott' = pidom(top'). Then (top',bott') ⊂ (top,bott) and (bott',top') ⊂ (top,bott). Call (top',bott')ᴰ a doubly factorable subgraph. 

Proof: since top' = idom(bott') and top ∈ (top,bott) (top',bott') ⊂ (top,bott).

Doubly factorable subgraphs, (top,bott)ᴰ have the property that they will not be destroyed by factoring any other subgraph. Destruction of subgraphs can only occur when there are frontier nodes. This will never occur for doubly factorable subgraphs because for all n ∈ (top,bott)ᴰ dom(top,n) & pdom(bott,n). So they can always be factored, and should be factored before their containing subgraph.

Singly factorable subgraphs should be factored before their enclosing subgraphs but they may be destroyed in the factorization process.
"""

"""
Frontier node test. Definition of a frontier node for the dom subgraph case (the pdom case is symmetric):
* If dom(a,b) & num_children(a) > 1 & num_parents(b) > 1 then subgraph (a,b) is a dom subgraph. 
* A node is contained in the subgraph if it is on the path from b to a. 
* A frontier node, nₚ, is node contained in (a,b) with !pdom(b,nₚ). Frontier nodes cause edges to be preserved after factorization.



"""

"""
The rules for factorization are based on a simple idea: the factored graph must include all possible paths between roots and variables. Edges can be deleted, (in this implementation have their root or variable bitmasks set to 0), only if doing so does not remove a necessary path. For example, in a factorable subgraph `(a,b)` where 
    `dom(a,b) = true & Rv(pdom(b,a)) = true` all edges contained in `(a,b)` can be removed from the graph and replaced by the sum of products evaluation of (a,b). This is because a path from any 
"""
struct EdgePreservationOverview end

"""
Rules for clearing root and variable bitmasks during factorization.


`R = roots(a)` and `V = variables(b)` are bitmasks with the property that `R[i] = true` for edge `e` then there is a path downward from root i to `top_vertex(e)`. Similarly if `V[j] = true` for edge `e` then there is a path upward from variable[j] to `bott_vertex(e)`.

Dom case, factorable subgraph `(a,b), a > b`, so that `dom(a,b) = true`.  The root path constraint `Rₚ` is the bitmask corresponding to the set of roots which have `(a,b)` as a factorable dom. `Rₚ(a,b)` is the set of edges reachable from `a` that are in the dom `(a,b)` and have all bits corresponding to `Rₚ` true. 

Evaluation of `(a,b)` and resetting bitmasks starts at node `b` and proceeds upward. For edge `e` encountered during this upward traversal, where `bott_vertex(e)≠b`, compute `tᵦ = bott_vertex(e)`. If `children(tᵦ)` constrained by `Rₚ` has only one element then `roots(e) = roots(e) & ! R`. In words, after factorization there will no longer be a path downward from the root set `R` through edge e. 

If `children(tᵦ)` constrained by `Rₚ` has more than one element then `roots(e)` and all edges above it until vertex `a` have their roots fields unchanged. 

Proof: 
Assume edge `e` has `bott_vertex(e)=b`. Then all paths downward from `e` must pass through `b`. No path from the root set `Rₚ` can pass through `e` without also passing through `b`. So `e` can be eliminated from the graph constrained by `Rₚ`: `reachable_roots(e)←reachable_roots(e) & !Rₚ`.

Assume now `e₁` is one edge upward from `b`, i.e. `bott_vertex(bott_vertex(e₁))=b`. If `length(children(Rₚ(e₁))) = 1` then there is only one path through `Rₚ(e₁)` to `b` and `reachable_roots(e₁)←reachable_roots(e₁) & !Rₚ`.

Proceed upward through the edges until `top_vertex(e)=a`. If any edge along the path fails the `length(children(Rₚ(e₁))) = 1` test then all edges above that point have `reachable_roots` unchanged because there is a path to some variables that does not pass through `b`.

A similar argument applies to downward edges for factorable pdom subgraph `(a,b), a<b`.

Starting from `a` and working downward if `length(parent_edges(top_vertex(e))) = 1` then there is no path downward to any variable except through `a`. Edge `e` will have its variable mask 




"""
struct BitmaskClearingRules end