#!/usr/bin/env julia
# Regenerates the fig-*.svg diagrams used by ALGORITHM_EXPLAINED.md to illustrate
# the derivative graph for the worked example in that document (@variables v1 v2
# v3; r1 = (p+s) + (n16*v1), see ALGORITHM_EXPLAINED.md section 2.1).
#
# Usage (from the repository root):
#     julia img/generate_figures.jl
# or, to write somewhere other than img/:
#     julia img/generate_figures.jl path/to/outdir
#
# All five diagrams are generated from ONE shared table of node (x, y)
# coordinates (POS, below), so a given node sits at the exact same pixel
# position in every picture it appears in - nothing needs Mermaid's
# auto-layout, which re-lays-out the whole graph on any change to the
# source and makes before/after pictures hard to compare by eye.
#
# The edge masks in STATE_A/STATE_B_OVERRIDES/STATE_C below are not
# hand-derived: they were read off FastDifferentiation's own DerivativeGraph
# for this example at three points - before any factoring, after the first
# factoring step ([b=2,c=7]), and after the complete `factor!` pass. If the
# worked example in ALGORITHM_EXPLAINED.md section 2.1 ever changes, re-derive
# these with a script such as:
#
#     import FastDifferentiation as FD
#     FD.@variables v1 v2 v3
#     d = v2 - v3; s = v1 * d; p = v1 + v2
#     n12 = v2 + s; n10 = p + s; n15 = n12 * p; n16 = v2 + n15; n18 = n16 * v1
#     r1 = n10 + n18
#     gr = FD.DerivativeGraph([r1])
#     # dump every edge's (reachable_roots, reachable_variables):
#     for node in 1:13
#         es = FD.node_edges(gr, node)
#         es === nothing && continue
#         for e in FD.children(es)
#             println(FD.top_vertex(e), "->", FD.bott_vertex(e), "  r:", FD.reachable_roots(e), "  v:", FD.reachable_variables(e))
#         end
#     end
#     # then pop one subgraph off FD.compute_factorable_subgraphs(gr) and call
#     # FD.factor_subgraph! on it (state B), or call FD.factor!(gr) directly
#     # (state C), and dump again.

const POS = Dict(
    13 => (430, 30),
    12 => (430, 120),
    11 => (150, 210),
    7  => (710, 210),
    10 => (150, 300),
    9  => (150, 390),
    8  => (150, 480),
    3  => (430, 480),
    6  => (430, 570),
    1  => (440, 660),
    5  => (620, 660),
    2  => (330, 750),
    4  => (610, 750),
)

# fill, stroke  (role coloring only; no masks in the box)
const STYLE = Dict(
    12 => ("#9fd3ff", "#00488f"),
    6  => ("#ffd27f", "#9a5b00"),
    1  => ("#c8f0c8", "#1f7a1f"),
    2  => ("#c8f0c8", "#1f7a1f"),
    4  => ("#c8f0c8", "#1f7a1f"),
)
default_style() = ("#f0f0f0", "#333333")
const LEAVES = (1, 2, 4)
const BOXW, BOXH = 70, 40

# ---- real reachable_roots/reachable_variables per edge, from the actual code ----
# STATE A: before any factoring
const STATE_A = Dict(
    (3,1)=>"v1", (3,2)=>"v2", (5,2)=>"v2", (5,4)=>"v3", (6,1)=>"v1", (6,5)=>"v2,v3",
    (7,3)=>"v1,v2", (7,6)=>"v1,v2,v3", (8,2)=>"v2", (8,6)=>"v1,v2,v3", (9,8)=>"v1,v2,v3",
    (9,3)=>"v1,v2", (10,2)=>"v2", (10,9)=>"v1,v2,v3", (11,10)=>"v1,v2,v3", (11,1)=>"v1",
    (12,7)=>"v1,v2,v3", (12,11)=>"v1,v2,v3", (13,12)=>"v1,v2,v3",
)
# STATE B: after step 1 ([2,7] factored) -- only entries that differ from A, plus the new edge
const STATE_B_OVERRIDES = Dict((7,3)=>"v1", (7,6)=>"v1,v3")
const STATE_B_NEW_EDGE = ((7,2), "v2")
# STATE C: after the entire factor! pass (verified against FD.factor!, not hand-derived)
const STATE_C = Dict(
    (6,5)=>"v3", (5,4)=>"v3", (13,12)=>"v1,v2,v3", (12,6)=>"v3", (12,2)=>"v2", (12,1)=>"v1",
)

const ALL_EDGES = [
    (13, 12), (12, 7), (12, 11), (7, 3), (7, 6), (11, 10), (11, 1),
    (10, 2), (10, 9), (9, 8), (9, 3), (8, 2), (8, 6), (6, 1), (6, 5),
    (3, 1), (3, 2), (5, 2), (5, 4),
]

r_label(varlabel) = "r:1 v:" * replace(varlabel, "v"=>"")

function node_svg(id)
    x, y = POS[id]
    fill, stroke = get(STYLE, id, default_style())
    label = string(id)
    if id in LEAVES
        shape = """<ellipse cx="$x" cy="$y" rx="$(BOXW÷2)" ry="$(BOXH÷2)" fill="$fill" stroke="$stroke" stroke-width="2"/>"""
    else
        shape = """<rect x="$(x-BOXW÷2)" y="$(y-BOXH÷2)" width="$BOXW" height="$BOXH" rx="8" ry="8" fill="$fill" stroke="$stroke" stroke-width="2"/>"""
    end
    text = """<text x="$x" y="$(y+5)" text-anchor="middle" font-family="Helvetica,Arial,sans-serif" font-size="16" font-weight="700" fill="#000000">$label</text>"""
    return shape * "\n    " * text
end

function edge_path(xa, ya, xb, yb; shorten_a=24, shorten_b=24)
    dx, dy = xb - xa, yb - ya
    len = sqrt(dx^2 + dy^2)
    ux, uy = dx / len, dy / len
    x1, y1 = xa + ux * shorten_a, ya + uy * shorten_a
    x2, y2 = xb - ux * shorten_b, yb - uy * shorten_b
    return x1, y1, x2, y2
end

"""
    edge_svg(a, b, posA, posB, label; color, width, label_frac, label_color, dashed, bow)

Returns `(line_svg, label_svg_or_nothing)` for one edge from node `a` at `posA`
to node `b` at `posB`. `label_frac` places the label at that fraction of the
way along the edge (0 = at `a`, 1 = at `b`); `bow` bows the edge into a
quadratic Bezier curve (perpendicular offset in pixels) to route it around a
node that would otherwise sit on the straight line between `a` and `b`.

Callers must draw ALL returned lines before ANY returned labels (see the
`diagram_*` functions below): if a label is drawn before every line, a later
edge's line can be drawn on top of an earlier edge's label wherever they
cross, obscuring it.
"""
function edge_svg(a, b, posA, posB, label; color="#333333", width=1.6, label_frac=0.5, label_color=nothing, dashed=false, bow=0.0)
    xa, ya = posA
    xb, yb = posB
    x1, y1, x2, y2 = edge_path(xa, ya, xb, yb)
    dash = dashed ? " stroke-dasharray=\"6,4\"" : ""
    local mx, my
    if bow == 0.0
        line = """<line x1="$x1" y1="$y1" x2="$x2" y2="$y2" stroke="$color" stroke-width="$width"$dash marker-end="url(#arrow)"/>"""
        mx, my = x1 + (x2 - x1) * label_frac, y1 + (y2 - y1) * label_frac
    else
        # quadratic bezier control point, offset perpendicular to the chord
        cxm, cym = (x1 + x2) / 2, (y1 + y2) / 2
        dx, dy = x2 - x1, y2 - y1
        len = sqrt(dx^2 + dy^2)
        px, py = -dy / len, dx / len   # unit perpendicular
        cx, cy = cxm + px * bow, cym + py * bow
        line = """<path d="M $x1,$y1 Q $cx,$cy $x2,$y2" fill="none" stroke="$color" stroke-width="$width"$dash marker-end="url(#arrow)"/>"""
        # point on the quadratic curve at parameter t=label_frac
        t = label_frac
        mx = (1-t)^2*x1 + 2*(1-t)*t*cx + t^2*x2
        my = (1-t)^2*y1 + 2*(1-t)*t*cy + t^2*y2
    end
    lbl = nothing
    if label !== nothing
        lc = label_color === nothing ? color : label_color
        w = length(label) * 6.2 + 8
        lbl = """<rect x="$(mx-w/2)" y="$(my-9)" width="$w" height="16" fill="#ffffff" opacity="0.92" stroke="$lc" stroke-width="0.6"/>
    <text x="$mx" y="$(my+3)" text-anchor="middle" font-family="Helvetica,Arial,sans-serif" font-size="10.5" fill="$lc">$label</text>"""
    end
    return line, lbl
end

const DEFS = """
  <defs>
    <marker id="arrow" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto" markerUnits="strokeWidth">
      <path d="M0,0 L6,3 L0,6 Z" fill="#333333"/>
    </marker>
  </defs>
"""

function wrap_svg(body; width=860, height=820)
    return """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 $width $height" width="$width" height="$height" font-family="Helvetica,Arial,sans-serif">
  <rect x="0" y="0" width="$width" height="$height" fill="#ffffff"/>
$DEFS
  <g>
$body
  </g>
</svg>
"""
end

allnodes() = sort(collect(keys(POS)))

default_label_fracs() = Dict(
    (11,1) => 0.3, (10,2) => 0.25, (9,3) => 0.65, (8,2) => 0.28,
    (12,7) => 0.68, (12,11) => 0.68, (7,3) => 0.22, (7,6) => 0.68,
)

"""
    diagram_plain(masks; extra_edge, highlight, focus, title)

The full 13-node graph: used (with `masks = STATE_A`) for "the graph before
factoring" and, with `focus`/`title` set, for the `[2,7]`-subgraph-highlighted
picture; used (with `masks` merging in `STATE_B_OVERRIDES` and `extra_edge =
STATE_B_NEW_EDGE`) for the picture after that step runs.
"""
function diagram_plain(masks; extra_edge=nothing, highlight=Set(), focus=nothing, title=nothing)
    bg, lines, labels, nodes = String[], String[], String[], String[]
    fracs = default_label_fracs()
    if focus !== nothing
        xs = [POS[n][1] for n in focus]
        ys = [POS[n][2] for n in focus]
        pad = 55; topPad = 40
        bx1, by1 = minimum(xs) - pad, minimum(ys) - topPad
        bx2, by2 = maximum(xs) + pad, maximum(ys) + pad
        push!(bg, """    <rect x="$bx1" y="$by1" width="$(bx2-bx1)" height="$(by2-by1)" rx="18" ry="18" fill="#7a1fa2" opacity="0.07" stroke="#7a1fa2" stroke-width="2" stroke-dasharray="10,6"/>""")
        if title !== nothing
            push!(bg, """    <text x="$bx1" y="$(by1-10)" font-family="Helvetica,Arial,sans-serif" font-size="13" font-weight="700" fill="#7a1fa2">$title</text>""")
        end
    end
    for (a, b) in ALL_EDGES
        lbl = r_label(masks[(a,b)])
        hl = (a,b) in highlight
        line, label = edge_svg(a, b, POS[a], POS[b], lbl;
            color = hl ? "#7a1fa2" : "#333333", width = hl ? 3.0 : 1.6,
            label_frac = get(fracs, (a,b), 0.5), label_color = hl ? "#7a1fa2" : "#555555")
        push!(lines, "    " * line)
        label !== nothing && push!(labels, "    " * label)
    end
    if extra_edge !== nothing
        (a, b), lbl = extra_edge
        line, label = edge_svg(a, b, POS[a], POS[b], r_label(lbl); color="#7a1fa2", width=3.0, label_frac=0.6, label_color="#7a1fa2")
        push!(lines, "    " * line)
        label !== nothing && push!(labels, "    " * label)
    end
    for id in allnodes()
        push!(nodes, "    " * node_svg(id))
    end
    return wrap_svg(join(vcat(bg, lines, labels, nodes), "\n"))
end

"""
    diagram_region()

The "region" between dominated node 6 and dominating node 12 that
`subgraph_region` finds for the `[12,6]` subgraph, in state B (after the
`[2,7]` step, before `[12,6]` runs): the short path through node 7 (blue),
the long path through 8-9-10-11 (red), and the dashed gray bypass edges to
the sibling operands (nodes 1, 2, 3) that sit outside the region.
"""
function diagram_region()
    # nodes in the [12,6] region plus the bypass targets shown dashed/outside
    spine = [(12,7),(7,6),(12,11),(11,10),(10,9),(9,8),(8,6)]
    bypass = [(7,3),(11,1),(10,2),(9,3),(8,2)]
    short_path = Set([(12,7),(7,6)])
    node_ids = (12,7,6,11,10,9,8,3,1,2)
    lines, labels, nodes = String[], String[], String[]
    fracs = Dict((12,7)=>0.68, (12,11)=>0.68, (7,6)=>0.72, (7,3)=>0.15,
                 (11,1)=>0.28, (10,2)=>0.24, (9,3)=>0.68, (8,2)=>0.26)
    masks = merge(STATE_A, STATE_B_OVERRIDES)
    for (a,b) in spine
        color = (a,b) in short_path ? "#0055aa" : "#cc2222"
        line, label = edge_svg(a, b, POS[a], POS[b], r_label(masks[(a,b)]);
            color=color, width=3.0, label_frac=get(fracs,(a,b),0.5), label_color=color)
        push!(lines, "    " * line)
        label !== nothing && push!(labels, "    " * label)
    end
    for (a,b) in bypass
        line, label = edge_svg(a, b, POS[a], POS[b], r_label(masks[(a,b)]);
            color="#999999", width=1.3, dashed=true, label_frac=get(fracs,(a,b),0.5), label_color="#777777")
        push!(lines, "    " * line)
        label !== nothing && push!(labels, "    " * label)
    end
    for id in node_ids
        push!(nodes, "    " * node_svg(id))
    end
    return wrap_svg(join(vcat(lines, labels, nodes), "\n"))
end

"""
    diagram_after_factor()

The graph after the complete `factor!` pass: only nodes 13, 12, 6, 5, 1, 2, 4
and 6 edges survive. `12->1` is drawn as a curve (`bow`) so it doesn't pass
directly through node 6's box, since both sit at nearly the same x coordinate.
"""
function diagram_after_factor()
    edges = [(13,12),(12,6),(12,2),(12,1),(6,5),(5,4)]
    node_ids = (13,12,6,5,1,2,4)
    bows = Dict((12,1) => -90.0)
    fracs = Dict((12,6) => 0.42, (12,1) => 0.55)
    lines, labels, nodes = String[], String[], String[]
    for (a,b) in edges
        line, label = edge_svg(a, b, POS[a], POS[b], r_label(STATE_C[(a,b)]);
            color="#333333", width=1.8, label_frac=get(fracs,(a,b),0.5), label_color="#555555",
            bow=get(bows,(a,b),0.0))
        push!(lines, "    " * line)
        label !== nothing && push!(labels, "    " * label)
    end
    for id in node_ids
        push!(nodes, "    " * node_svg(id))
    end
    return wrap_svg(join(vcat(lines, labels, nodes), "\n"))
end

outdir = isempty(ARGS) ? joinpath(@__DIR__) : ARGS[1]
mkpath(outdir)

write(joinpath(outdir, "fig-graph-before.svg"), diagram_plain(STATE_A))

focus = (7, 3, 6, 5, 2)
highlighted = Set([(7,3),(3,2),(7,6),(6,5),(5,2)])
write(joinpath(outdir, "fig-step1-highlight.svg"),
      diagram_plain(STATE_A; highlight=highlighted, focus=focus, title="subgraph about to be factored: [b=2, c=7]"))

masks_b = merge(STATE_A, STATE_B_OVERRIDES)
write(joinpath(outdir, "fig-step1-after.svg"),
      diagram_plain(masks_b; extra_edge=(STATE_B_NEW_EDGE[1], STATE_B_NEW_EDGE[2]), highlight=Set([(7,3),(7,6)])))

write(joinpath(outdir, "fig-region.svg"), diagram_region())
write(joinpath(outdir, "fig-after-factor.svg"), diagram_after_factor())

println("wrote 5 SVGs to ", outdir)
