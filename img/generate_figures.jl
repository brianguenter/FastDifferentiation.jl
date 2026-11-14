#!/usr/bin/env julia
# Regenerates the fig-*.svg diagrams in HOW_FASTDIFFERENTIATION_WORKS.md, which
# walk through differentiating the example
#
#     @variables x y
#     p = x + y; q = x + p; r = p + q
#     jacobian([q, r], [x, y])
#
# Usage (from the repository root):
#     julia --project=. img/generate_figures.jl            # writes into img/
#     julia --project=. img/generate_figures.jl outdir     # writes elsewhere
#
# Every diagram uses the same node positions (POS), so a node sits at the same
# place in every picture. The edge states below are written out by hand, but
# before writing anything the script rebuilds the graph with FastDifferentiation,
# runs the factoring steps one at a time, and checks that the drawn postorder
# numbers, edge states (initial, after step 1, after step 2, after step 5) and
# group tables match the code exactly. The two states of the earlier,
# one-edge-per-subgraph rule (MAIN_AFTER_2, MAIN_FINAL) were read off a run of
# FastDifferentiation 0.4.5 and cannot be checked against the current code.

import FastDifferentiation as FD

# ---------------------------------------------------------------- layout
const POS = Dict("r" => (200, 50), "q" => (400, 160), "p" => (200, 270), "y" => (80, 390), "x" => (400, 390))
const RAD = 23
const LABEL_POS = Dict(
    ("r", "p") => (140, 160), ("r", "q") => (318, 64), ("q", "p") => (322, 250), ("q", "x") => (462, 275),
    ("p", "x") => (262, 364), ("p", "y") => (96, 302), ("r", "x") => (548, 214))
const CAPTION_SIDE = Dict("r" => :left, "q" => :right, "p" => :left, "x" => :right, "y" => :left)

# light palette; the figures carry their own white background so they also read on dark pages
const C = Dict(
    :bg => "#ffffff", :fg => "#16202b", :muted => "#586473", :line => "#d4dbe3", :edge => "#3d4856",
    :faded => "#bcc4ce", :accent => "#0b7285", :q => "#1f4fd1", :r => "#b4520a", :var => "#3f7a12",
    :a => "#0b7285", :b => "#8b2fc9", :new => "#0b7285", :bad => "#c0262d")
const WIDTH = Dict(:normal => 1.8, :faded => 1.6, :a => 4.0, :b => 4.0, :new => 3.4, :bad => 3.4)
const NAME_COLOR = Dict('q' => C[:q], 'r' => C[:r], 'x' => C[:var], 'y' => C[:var])
const MONO = "'JetBrains Mono', Consolas, 'DejaVu Sans Mono', Menlo, monospace"
const SANS = "'Space Grotesk', 'Segoe UI', Helvetica, Arial, sans-serif"

# ---------------------------------------------------------------- edge states
# (upper, lower, value, root mask, variable mask, style; options)
struct Edge
    u::String; v::String; val::Int; R::String; V::String; style::Symbol
    label::Union{Nothing,Tuple{Int,Int}}   # label position override
    bend::Float64                          # perpendicular bow, pixels
    ctrl::Union{Nothing,Tuple{Int,Int}}    # explicit quadratic control point
end
E(u, v, val, R, V, style=:normal; label=nothing, bend=0.0, ctrl=nothing) = Edge(u, v, val, R, V, style, label, bend, ctrl)
restyle(e::Edge, s::Symbol) = Edge(e.u, e.v, e.val, e.R, e.V, s, e.label, e.bend, e.ctrl)

const INITIAL = [E("r", "p", 1, "r", "xy"), E("r", "q", 1, "r", "xy"), E("q", "x", 1, "qr", "x"),
                 E("q", "p", 1, "qr", "xy"), E("p", "x", 1, "qr", "x"), E("p", "y", 1, "qr", "y")]
const AFTER_1 = [E("r", "p", 2, "r", "xy", :new), E("r", "q", 1, "r", "x", :new), E("q", "p", 1, "q", "xy", :new),
                 E("q", "x", 1, "qr", "x"), E("p", "x", 1, "qr", "x"), E("p", "y", 1, "qr", "y")]
const AFTER_2 = [E("r", "p", 2, "r", "xy"), E("r", "q", 1, "r", "x"),
                 E("q", "x", 2, "q", "x", :new; label=(345, 300)), E("q", "x", 1, "r", "x", :new; bend=-62.0, label=(506, 272)),
                 E("q", "p", 1, "q", "y", :new), E("p", "x", 1, "r", "x", :new), E("p", "y", 1, "qr", "y")]
const FINAL = [E("r", "p", 2, "r", "y", :new), E("r", "x", 3, "r", "x", :new; ctrl=(620, 140)),
               E("q", "x", 2, "q", "x"; label=(345, 300)), E("q", "p", 1, "q", "y"), E("p", "y", 1, "qr", "y")]
const MAIN_AFTER_2 = [E("r", "p", 2, "r", "xy"), E("r", "q", 1, "r", "x"), E("q", "x", 2, "qr", "x", :bad),
                      E("q", "p", 1, "q", "y", :new), E("p", "x", 1, "r", "x", :new), E("p", "y", 1, "qr", "y")]
const MAIN_FINAL = [E("r", "p", 2, "r", "y"), E("r", "x", 4, "r", "x", :bad; ctrl=(620, 140)),
                    E("q", "x", 2, "q", "x"; label=(345, 300)), E("q", "p", 1, "q", "y"), E("p", "y", 1, "qr", "y")]

# ---------------------------------------------------------------- check against the code
function code_states()
    FD.@variables x y
    p = x + y; q = x + p; r = p + q
    gr = FD.DerivativeGraph([q, r])
    name = Dict(FD.postorder_number(gr, n) => s for (n, s) in ((x, "x"), (y, "y"), (p, "p"), (q, "q"), (r, "r")))
    # the postorder numbers drawn in fig-postorder.svg
    name == Dict(1 => "x", 2 => "y", 3 => "p", 4 => "q", 6 => "r") || error("postorder numbers differ from fig-postorder.svg: $name")
    (FD.is_NoOp(FD.node(gr, 5)) && FD.is_NoOp(FD.node(gr, 7))) || error("nodes 5 and 7 were expected to be the NoOp wrappers")
    function snapshot()
        out = Set{Tuple{String,String,Int,String,String}}()
        for i in keys(name)
            es = FD.node_edges(gr, i)
            es === nothing && continue
            for e in FD.children(es)
                haskey(name, FD.top_vertex(e)) || continue   # skip the NoOp wrapper edges above the roots
                R = join(["qr"[k] for k in 1:2 if FD.reachable_roots(e)[k]])
                V = join(["xy"[k] for k in 1:2 if FD.reachable_variables(e)[k]])
                push!(out, (name[FD.top_vertex(e)], name[FD.bott_vertex(e)], Int(FD.value(FD.value(e))), R, V))
            end
        end
        out
    end
    # The group tables G(u) of `subgraph_path_groups`, as (root mask, variable mask, V, n)
    # entries. For a dominator subgraph the memo holds every node's table; for a
    # postdominator subgraph only the final table, at the dominating node, is returned.
    function group_tables(sub)
        isdom = sub isa FD.FactorableSubgraph{Int64,FD.DominatorSubgraph}
        function entry(g)
            dm, nm, V, n = g
            roots, vars = isdom ? (dm, nm) : (nm, dm)
            (join(["qr"[k] for k in 1:2 if roots[k]]), join(["xy"[k] for k in 1:2 if vars[k]]), Int(FD.value(V)), n)
        end
        out = Dict{String,Set{Tuple{String,String,Int,Int}}}()
        if isdom
            memo = Dict{Int64,FD.PathGroups}()
            FD.subgraph_path_groups(sub, FD.dominated_node(sub), memo)
            for (u, groups) in memo
                out[name[u]] = Set(entry(g) for g in values(groups))
            end
        else
            groups = FD.subgraph_path_groups(sub, FD.subgraph_region(sub))
            out[name[FD.dominating_node(sub)]] = Set(entry(g) for g in values(groups))
        end
        return out
    end

    states = Dict(0 => snapshot())
    tables = Dict{Int,Dict{String,Set{Tuple{String,String,Int,Int}}}}()
    heap = FD.compute_factorable_subgraphs(gr)
    step = 0
    while !isempty(heap)
        sub = pop!(heap)
        step += 1
        step in (1, 2) && (tables[step] = group_tables(sub))
        FD.factor_subgraph!(sub)
        states[step] = snapshot()
    end
    return states, tables
end

as_set(edges) = Set((e.u, e.v, e.val, e.R, e.V) for e in edges)

# group tables drawn in fig-dp-step1.svg and fig-dp-step2.svg: node => [(root mask, variable mask, V, n)]
const DP_STEP1 = Dict("r" => [("r", "xy", 1, 1)], "q" => [("r", "xy", 1, 1)], "p" => [("r", "xy", 2, 2)])
const DP_STEP2 = Dict("q" => [("qr", "x", 1, 1)], "p" => [("q", "x", 1, 1)], "x" => [("qr", "x", 1, 1), ("q", "x", 1, 1)])

let (states, tables) = code_states()
    for (step, drawn) in ((0, INITIAL), (1, AFTER_1), (2, AFTER_2), (5, FINAL), (6, FINAL))
        states[step] == as_set(drawn) || error("figure state after step $step does not match the code:\n code:  $(sort(collect(states[step])))\n drawn: $(sort(collect(as_set(drawn))))")
    end
    states[3] == states[2] && states[4] == states[2] || error("steps 3 and 4 were expected to change nothing")
    # step 1 is a dominator subgraph, so every table can be checked; step 2 only at x
    tables[1] == Dict(k => Set(v) for (k, v) in DP_STEP1) || error("step 1 group tables do not match the code: $(tables[1])")
    tables[2] == Dict("x" => Set(DP_STEP2["x"])) || error("step 2 group table at x does not match the code: $(tables[2])")
end

# ---------------------------------------------------------------- drawing
f(x) = string(round(x, digits=1))

"Text spans for `s`, with the names q, r, x, y in their colours."
colored_spans(s) = join(c in keys(NAME_COLOR) ? """<tspan fill="$(NAME_COLOR[c])" font-weight="700">$c</tspan>""" :
                        """<tspan fill="$(C[:fg])">$c</tspan>""" for c in s)

"A box beside node `n` listing its group table: one line per entry, `roots | variables   V=…, n=…`."
function table_box(n, entries, side)
    (cx, cy) = POS[n]
    rows = ["$(join(collect(R), " ")) | $(join(collect(V), " "))   V=$val, n=$cnt" for (R, V, val, cnt) in entries]
    w = maximum(length, vcat(rows, ["G($n)"])) * 7.9 + 18
    h = 20 * (length(rows) + 1) + 8
    x0 = side == :left ? cx - RAD - 14 - w : cx + RAD + 14
    y0 = clamp(cy - h / 2, 4, 426 - h)
    parts = ["""<rect x="$(f(x0))" y="$(f(y0))" width="$(f(w))" height="$h" rx="6" fill="$(C[:bg])" stroke="$(C[:accent])" stroke-width="1.6"/>""",
             """<text x="$(f(x0 + 9))" y="$(f(y0 + 19))" font-family="$MONO" font-size="13" font-weight="700" fill="$(C[:accent])">G($n)</text>"""]
    for (i, row) in enumerate(rows)
        push!(parts, """<text x="$(f(x0 + 9))" y="$(f(y0 + 19 + 20i))" font-family="$MONO" font-size="13" xml:space="preserve">$(colored_spans(row))</text>""")
    end
    return join(parts, "\n  ")
end

function svg(edges::Vector{Edge}; rings=Tuple{String,String}[], values_only=false, tables=[], faded_labels=true)
    lines = String[]; labels = String[]; nodes = String[]
    for (n, role) in rings
        (cx, cy) = POS[n]
        dash = role == "dominating" ? " stroke-dasharray=\"5 4\"" : ""
        push!(nodes, """<circle cx="$cx" cy="$cy" r="$(RAD + 9)" fill="none" stroke="$(C[:accent])" stroke-width="2"$dash/>""")
        left = CAPTION_SIDE[n] == :left
        tx = left ? cx - (RAD + 16) : cx + (RAD + 16)
        push!(nodes, """<text x="$tx" y="$(cy + 4)" font-family="$MONO" font-size="12" letter-spacing="0.5" fill="$(C[:accent])" text-anchor="$(left ? "end" : "start")">$role</text>""")
    end
    for e in edges
        (ax, ay) = POS[e.u]; (bx, by) = POS[e.v]
        dx, dy = bx - ax, by - ay; len = hypot(dx, dy)
        ctrl = e.ctrl !== nothing ? e.ctrl :
               e.bend != 0 ? ((ax + bx) / 2 - dy / len * e.bend, (ay + by) / 2 + dx / len * e.bend) : nothing
        if ctrl === nothing
            s0 = (ax + dx / len * RAD, ay + dy / len * RAD)
            s1 = (bx - dx / len * (RAD + 2), by - dy / len * (RAD + 2))
            d = "M$(f(s0[1])),$(f(s0[2])) L$(f(s1[1])),$(f(s1[2]))"
        else
            u0 = (ctrl[1] - ax, ctrl[2] - ay); l0 = hypot(u0...)
            u1 = (ctrl[1] - bx, ctrl[2] - by); l1 = hypot(u1...)
            s0 = (ax + u0[1] / l0 * RAD, ay + u0[2] / l0 * RAD)
            s1 = (bx + u1[1] / l1 * (RAD + 2), by + u1[2] / l1 * (RAD + 2))
            d = "M$(f(s0[1])),$(f(s0[2])) Q$(f(ctrl[1])),$(f(ctrl[2])) $(f(s1[1])),$(f(s1[2]))"
        end
        col = e.style == :normal ? C[:edge] : C[e.style]
        push!(lines, """<path d="$d" fill="none" stroke="$col" stroke-width="$(WIDTH[e.style])" stroke-linecap="round" marker-end="url(#arrow-$(e.style))"/>""")
        (e.style == :faded && !faded_labels) && continue

        # label: value, then root mask | variable mask
        parts = Tuple{String,String,Bool}[(string(e.val), e.style == :bad ? C[:bad] : C[:fg], true)]
        if !values_only
            push!(parts, ("  ", C[:fg], false))
            for (i, ch) in enumerate(e.R)
                push!(parts, ((i > 1 ? " " : "") * ch, NAME_COLOR[ch], true))
            end
            push!(parts, (" | ", C[:muted], false))
            for (i, ch) in enumerate(e.V)
                push!(parts, ((i > 1 ? " " : "") * ch, NAME_COLOR[ch], true))
            end
        end
        nchar = sum(length(p[1]) for p in parts)
        w, h = nchar * 7.9 + 16, 22
        (lx, ly) = e.label !== nothing ? e.label : LABEL_POS[(e.u, e.v)]
        if values_only # a small label: put it just beside the edge's midpoint, on the same side as the full label
            mx, my = (ax + bx) / 2, (ay + by) / 2
            dl = hypot(lx - mx, ly - my)
            lx, ly = mx + (lx - mx) / dl * 18, round(Int, my + (ly - my) / dl * 18)
        end
        border = e.style in (:new, :bad, :a, :b) ? C[e.style] : C[:line]
        opacity = e.style == :faded ? " opacity=\"0.55\"" : ""
        spans = join(["""<tspan fill="$c"$(b ? " font-weight=\"700\"" : "")>$t</tspan>""" for (t, c, b) in parts])
        push!(labels, """<g$opacity><rect x="$(f(lx - w / 2))" y="$(ly - h / 2)" width="$(f(w))" height="$h" rx="5" fill="$(C[:bg])" stroke="$border" stroke-width="$(border == C[:line] ? 1 : 1.8)"/>""" *
                      """<text x="$(f(lx - w / 2 + 8))" y="$(ly + 5)" font-family="$MONO" font-size="13" xml:space="preserve">$spans</text></g>""")
    end
    for (n, (cx, cy)) in sort(collect(POS))
        isvar = n in ("x", "y"); isroot = n in ("q", "r")
        stroke = isroot ? C[Symbol(n)] : isvar ? C[:var] : C[:edge]
        if isvar
            push!(nodes, """<rect x="$(cx - RAD)" y="$(cy - RAD)" width="$(2RAD)" height="$(2RAD)" rx="9" fill="$(C[:bg])" stroke="$stroke" stroke-width="2"/>""")
        else
            push!(nodes, """<circle cx="$cx" cy="$cy" r="$RAD" fill="$(C[:bg])" stroke="$stroke" stroke-width="$(isroot ? 3 : 2)"/>""")
        end
        text_color = (isroot || isvar) ? stroke : C[:fg]
        push!(nodes, """<text x="$cx" y="$(cy + 7)" font-family="$SANS" font-size="19" font-weight="700" text-anchor="middle" fill="$text_color">$n</text>""")
    end
    markers = join(["""<marker id="arrow-$s" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="9" markerHeight="9" markerUnits="userSpaceOnUse" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="$(s == :normal ? C[:edge] : C[s])"/></marker>""" for s in keys(WIDTH)], "\n    ")
    return """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 640 430" width="640" height="430">
  <defs>
    $markers
  </defs>
  <rect x="0" y="0" width="640" height="430" rx="12" fill="$(C[:bg])"/>
  $(join(lines, "\n  "))
  $(join(nodes, "\n  "))
  $(join(labels, "\n  "))
  $(join([table_box(n, entries, side) for (n, entries, side) in tables], "\n  "))
</svg>
"""
end

"""
The k x k grid of scripts/scalable_examples.jl (g[i,j] = g[i-1,j] * g[i,j-1]),
drawn as a diamond with the root g[k,k] at the top and x = g[1,1] at the bottom.
One of its binomial(2k-2, k-1) paths is highlighted, and the two corner squares,
the first subgraphs factored (both have the smallest node difference), are
shaded. Used in Appendix A to show a graph with exponentially many paths that
factoring still handles a few at a time."""
function grid_svg(k=4)
    dx, dy, top = 72, 52, 46
    pos(i, j) = (300 + (j - i) * dx, top + ((k - i) + (k - j)) * dy)
    height = top + 2 * (k - 1) * dy + 44
    # highlighted path: alternate down-left and down-right from the root
    path = [(k, k)]
    while path[end] != (1, 1)
        (i, j) = path[end]
        step = (isodd(length(path)) && i > 1) || j == 1 ? (i - 1, j) : (i, j - 1)
        push!(path, step)
    end
    on_path = Set(zip(path[1:end-1], path[2:end]))
    parts = String[]
    # the first subgraphs factored: the 2 x 2 corners [g[k,k], g[k-1,k-1]] and [x, g[2,2]]
    for corner in ([pos(2, 2), pos(2, 1), pos(1, 1), pos(1, 2)], [pos(k, k), pos(k, k - 1), pos(k - 1, k - 1), pos(k - 1, k)])
        push!(parts, """<polygon points="$(join(["$(a),$(b)" for (a, b) in corner], " "))" fill="$(C[:accent])" fill-opacity="0.10" stroke="$(C[:accent])" stroke-width="1.2" stroke-dasharray="4 3"/>""")
    end
    for i in 1:k, j in 1:k
        for (ci, cj) in ((i - 1, j), (i, j - 1))
            (ci < 1 || cj < 1) && continue
            (ax, ay) = pos(i, j); (bx, by) = pos(ci, cj)
            ddx, ddy = bx - ax, by - ay; len = hypot(ddx, ddy)
            hl = ((i, j), (ci, cj)) in on_path
            col = hl ? C[:a] : C[:edge]
            push!(parts, """<path d="M$(f(ax + ddx / len * 15)),$(f(ay + ddy / len * 15)) L$(f(bx - ddx / len * 17)),$(f(by - ddy / len * 17))" stroke="$col" stroke-width="$(hl ? 3.5 : 1.6)" fill="none" stroke-linecap="round" marker-end="url(#arrow-$(hl ? "a" : "normal"))"/>""")
        end
    end
    for i in 1:k, j in 1:k
        (x, y) = pos(i, j)
        if (i, j) == (1, 1)
            push!(parts, """<rect x="$(x - 16)" y="$(y - 16)" width="32" height="32" rx="7" fill="$(C[:bg])" stroke="$(C[:var])" stroke-width="2"/>""")
            push!(parts, """<text x="$x" y="$(y + 6)" font-family="$SANS" font-size="16" font-weight="700" text-anchor="middle" fill="$(C[:var])">x</text>""")
        else
            isroot = (i, j) == (k, k)
            push!(parts, """<circle cx="$x" cy="$y" r="14" fill="$(C[:bg])" stroke="$(isroot ? C[:r] : C[:edge])" stroke-width="$(isroot ? 3 : 2)"/>""")
        end
    end
    (rx, ry) = pos(k, k)
    push!(parts, """<text x="$(rx + 24)" y="$(ry + 5)" font-family="$MONO" font-size="13" fill="$(C[:fg])">root g[$k,$k]</text>""")
    for (i, j) in ((2, 2), (k - 1, k - 1))
        (sx, sy) = pos(i, j)
        push!(parts, """<text x="$(sx + 24)" y="$(sy + 5)" font-family="$MONO" font-size="12" fill="$(C[:accent])">factored first</text>""")
    end
    push!(parts, """<text x="20" y="$(height - 14)" font-family="$MONO" font-size="13" fill="$(C[:a])">highlighted: one of $(binomial(2k - 2, k - 1)) paths from the root to x</text>""")
    markers = join(["""<marker id="arrow-$s" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="9" markerHeight="9" markerUnits="userSpaceOnUse" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="$(C[s == "a" ? :a : :edge])"/></marker>""" for s in ("a", "normal")], "")
    return """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 640 $height" width="640" height="$height">
  <defs>$markers</defs>
  <rect x="0" y="0" width="640" height="$height" rx="12" fill="$(C[:bg])"/>
  $(join(parts, "\n  "))
</svg>
"""
end
"""
The example graph with the NoOp wrapper nodes the code puts above each root, and
every node's postorder number in a badge. The numbers are checked against the
code before drawing (see the `let` block that calls `code_states`)."""
function postorder_svg()
    pos = Dict("w7" => (200, 45), "r" => (200, 150), "w5" => (480, 150), "q" => (400, 250),
               "p" => (200, 340), "y" => (80, 450), "x" => (400, 450))
    number = Dict("x" => 1, "y" => 2, "p" => 3, "q" => 4, "w5" => 5, "r" => 6, "w7" => 7)
    edges = [("w7", "r"), ("w5", "q"), ("r", "p"), ("r", "q"), ("q", "p"), ("q", "x"), ("p", "x"), ("p", "y")]
    parts = String[]
    for (u, v) in edges
        (ax, ay) = pos[u]; (bx, by) = pos[v]
        dx, dy = bx - ax, by - ay; len = hypot(dx, dy)
        r0 = startswith(u, "w") ? 17 : RAD; r1 = RAD + 2
        push!(parts, """<path d="M$(f(ax + dx / len * r0)),$(f(ay + dy / len * r0)) L$(f(bx - dx / len * r1)),$(f(by - dy / len * r1))" stroke="$(C[:edge])" stroke-width="1.8" fill="none" marker-end="url(#arrow-normal)"/>""")
    end
    for (n, (cx, cy)) in sort(collect(pos))
        if startswith(n, "w")
            push!(parts, """<rect x="$(cx - 34)" y="$(cy - 16)" width="68" height="32" rx="16" fill="$(C[:bg])" stroke="$(C[:muted])" stroke-width="1.6" stroke-dasharray="4 3"/>""")
            push!(parts, """<text x="$cx" y="$(cy + 5)" font-family="$MONO" font-size="13" text-anchor="middle" fill="$(C[:muted])">NoOp</text>""")
            bx, by = cx + 34, cy - 16
        else
            isvar = n in ("x", "y"); isroot = n in ("q", "r")
            stroke = isroot ? C[Symbol(n)] : isvar ? C[:var] : C[:edge]
            if isvar
                push!(parts, """<rect x="$(cx - RAD)" y="$(cy - RAD)" width="$(2RAD)" height="$(2RAD)" rx="9" fill="$(C[:bg])" stroke="$stroke" stroke-width="2"/>""")
            else
                push!(parts, """<circle cx="$cx" cy="$cy" r="$RAD" fill="$(C[:bg])" stroke="$stroke" stroke-width="$(isroot ? 3 : 2)"/>""")
            end
            text_color = (isroot || isvar) ? stroke : C[:fg]
            push!(parts, """<text x="$cx" y="$(cy + 7)" font-family="$SANS" font-size="19" font-weight="700" text-anchor="middle" fill="$text_color">$n</text>""")
            # beside the node, on a side no edge touches
            side = Dict("r" => -1, "p" => -1, "y" => -1, "q" => 1, "x" => 1)[n]
            bx, by = cx + side * (RAD + 13), cy
        end
        push!(parts, """<circle cx="$bx" cy="$by" r="11" fill="$(C[:accent])"/>""")
        push!(parts, """<text x="$bx" y="$(by + 4.5)" font-family="$MONO" font-size="12" font-weight="700" text-anchor="middle" fill="#ffffff">$(number[n])</text>""")
    end
    return """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 640 490" width="640" height="490">
  <defs><marker id="arrow-normal" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="9" markerHeight="9" markerUnits="userSpaceOnUse" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="$(C[:edge])"/></marker></defs>
  <rect x="0" y="0" width="640" height="490" rx="12" fill="$(C[:bg])"/>
  $(join(parts, "\n  "))
</svg>
"""
end

outdir = isempty(ARGS) ? (@__DIR__) : ARGS[1]
mkpath(outdir)
seg(edges, styles) = [restyle(e, get(styles, (e.u, e.v, e.val, e.R), :faded)) for e in edges]

figures = [
    "fig-graph.svg" => svg(INITIAL; values_only=true),
    "fig-masks.svg" => svg(INITIAL),
    "fig-step1.svg" => svg(seg(INITIAL, Dict(("r", "p", 1, "r") => :a, ("r", "q", 1, "r") => :b, ("q", "p", 1, "qr") => :b));
                           rings=[("r", "dominating"), ("p", "dominated")]),
    "fig-step1-after.svg" => svg(AFTER_1),
    "fig-step2.svg" => svg(seg([restyle(e, :normal) for e in AFTER_1], Dict(("q", "x", 1, "qr") => :a, ("q", "p", 1, "q") => :b, ("p", "x", 1, "qr") => :b));
                           rings=[("x", "dominating"), ("q", "dominated")]),
    "fig-step2-main.svg" => svg(MAIN_AFTER_2),
    "fig-step2-new.svg" => svg(AFTER_2),
    "fig-step5.svg" => svg(seg(AFTER_2, Dict(("r", "q", 1, "r") => :a, ("q", "x", 1, "r") => :a, ("r", "p", 2, "r") => :b, ("p", "x", 1, "r") => :b));
                           rings=[("r", "dominating"), ("x", "dominated")]),
    "fig-final.svg" => svg(FINAL),
    "fig-final-main.svg" => svg(MAIN_FINAL),
    "fig-grid.svg" => grid_svg(4),
    "fig-postorder.svg" => postorder_svg(),
    "fig-dp-step1.svg" => svg(seg(INITIAL, Dict(("r", "p", 1, "r") => :normal, ("r", "q", 1, "r") => :normal, ("q", "p", 1, "qr") => :normal));
                              faded_labels=false,
                              tables=[("r", DP_STEP1["r"], :left), ("q", DP_STEP1["q"], :right), ("p", DP_STEP1["p"], :left)]),
    "fig-dp-step2.svg" => svg(seg([restyle(e, :normal) for e in AFTER_1], Dict(("q", "x", 1, "qr") => :normal, ("q", "p", 1, "q") => :normal, ("p", "x", 1, "qr") => :normal));
                              faded_labels=false,
                              tables=[("q", DP_STEP2["q"], :right), ("p", DP_STEP2["p"], :left), ("x", DP_STEP2["x"], :right)]),
]
for (file, content) in figures
    write(joinpath(outdir, file), content)
end
println("checked the drawn states against FastDifferentiation and wrote $(length(figures)) SVGs to ", outdir)
