# One small SVG line chart of a measured quantity against the size of the
# derivative graph, shared by scripts/measure_dominators.jl and
# scripts/measure_segments.jl.
#
# Light palette, matching the other figures: the series uses categorical slot 1
# of the dataviz reference palette; text and axes use neutral text tokens.

const SURFACE = "#ffffff"
const SERIES = "#2a78d6"
const TEXT = "#16202b"
const MUTED = "#586473"
const GRID = "#e3e7ec"
const REF = "#9aa4b0"
const SANS = "'Space Grotesk', 'Segoe UI', Helvetica, Arial, sans-serif"
const MONO = "'JetBrains Mono', Consolas, 'DejaVu Sans Mono', Menlo, monospace"

"About five round tick values from 0 to at least `hi`; whole numbers only if `integer`."
function ticks(hi; integer=false)
    hi <= 0 && return [0]
    raw = hi / 4
    mag = 10.0^floor(log10(raw))
    step = first(s for s in (1, 2, 2.5, 5, 10) .* mag if s >= raw)
    integer && (step = max(1, round(Int, step)))
    return collect(0:step:(ceil(hi / step) * step))
end
function fmt(v)
    v < 10_000 && return replace(string(round(Int, v)), r"(?<=\d)(?=(\d{3})+$)" => ",")
    k = round(v / 1000, digits=v >= 100_000 ? 0 : 1)
    return (isinteger(k) ? string(Int(k)) : string(k)) * "K"
end
r1(v) = string(round(v, digits=1))

"""
    growth_chart(title, xs, ys; quantity, end_label, zero_message)

`xs` are graph sizes (nodes + edges), `ys` the measured quantity, smallest size
first. A dashed grey reference line grows in proportion to size and passes
through the largest example; it is left out when every value is 0, in which case
`zero_message` is shown instead (and no end label). `end_label` labels the largest
example. `integer` keeps the vertical ticks to whole numbers, for counts."""
function growth_chart(title, xs, ys; quantity, end_label, zero_message, integer=false)
    W, H = 480, 300
    L, R, T, B = 66, 22, 62, 46
    xt = ticks(maximum(xs))
    yt = ticks(1.15 * maximum(ys); integer)   # headroom so the end label clears the legend
    ytop = maximum(ys) > 0 ? last(yt) : 1
    X(v) = L + (W - L - R) * v / last(xt)
    Yc(v) = H - B - (H - T - B) * v / ytop
    p = String[]
    push!(p, """<rect width="$W" height="$H" rx="12" fill="$SURFACE"/>""")
    push!(p, """<text x="16" y="26" font-family="$SANS" font-size="15" font-weight="700" fill="$TEXT">$title</text>""")
    # legend: measured series, and the proportional reference
    push!(p, """<circle cx="22" cy="44" r="4" fill="$SERIES" stroke="$SURFACE" stroke-width="2"/>""")
    push!(p, """<text x="32" y="48" font-family="$SANS" font-size="12" fill="$MUTED">$quantity</text>""")
    legend_x = 32 + 7 * length(quantity) + 14
    if maximum(ys) > 0
        push!(p, """<line x1="$legend_x" y1="44" x2="$(legend_x + 22)" y2="44" stroke="$REF" stroke-width="1.5" stroke-dasharray="4 3"/>""")
        push!(p, """<text x="$(legend_x + 28)" y="48" font-family="$SANS" font-size="12" fill="$MUTED">proportional to size</text>""")
    end
    for v in yt
        push!(p, """<line x1="$L" y1="$(r1(Yc(v)))" x2="$(W - R)" y2="$(r1(Yc(v)))" stroke="$GRID" stroke-width="1"/>""")
        push!(p, """<text x="$(L - 8)" y="$(r1(Yc(v) + 4))" font-family="$MONO" font-size="11" text-anchor="end" fill="$MUTED">$(fmt(v))</text>""")
    end
    for v in xt
        push!(p, """<line x1="$(r1(X(v)))" y1="$(H - B)" x2="$(r1(X(v)))" y2="$(H - B + 4)" stroke="$MUTED" stroke-width="1"/>""")
        push!(p, """<text x="$(r1(X(v)))" y="$(H - B + 17)" font-family="$MONO" font-size="11" text-anchor="middle" fill="$MUTED">$(fmt(v))</text>""")
    end
    push!(p, """<line x1="$L" y1="$(H - B)" x2="$(W - R)" y2="$(H - B)" stroke="$MUTED" stroke-width="1"/>""")
    push!(p, """<text x="$(L + (W - L - R) / 2)" y="$(H - 8)" font-family="$SANS" font-size="12" text-anchor="middle" fill="$MUTED">size of the derivative graph (nodes + edges)</text>""")
    if maximum(ys) > 0
        push!(p, """<line x1="$(X(0))" y1="$(Yc(0))" x2="$(r1(X(last(xs))))" y2="$(r1(Yc(last(ys))))" stroke="$REF" stroke-width="1.5" stroke-dasharray="4 3"/>""")
    else
        push!(p, """<text x="$(L + (W - L - R) / 2)" y="$(r1(Yc(0.5)))" font-family="$SANS" font-size="12" text-anchor="middle" fill="$MUTED">$zero_message</text>""")
    end
    pts = join(["$(r1(X(x))),$(r1(Yc(y)))" for (x, y) in zip(xs, ys)], " ")
    push!(p, """<polyline points="$pts" fill="none" stroke="$SERIES" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>""")
    for (x, y) in zip(xs, ys)
        push!(p, """<circle cx="$(r1(X(x)))" cy="$(r1(Yc(y)))" r="4" fill="$SERIES" stroke="$SURFACE" stroke-width="2"/>""")
    end
    # direct label on the largest example only
    maximum(ys) > 0 && push!(p, """<text x="$(r1(X(last(xs)) - 8))" y="$(r1(Yc(last(ys)) - 10))" font-family="$MONO" font-size="11" text-anchor="end" fill="$TEXT">$end_label</text>""")
    return """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 $W $H" width="$W" height="$H" role="img" aria-label="$title: $quantity against graph size">
  $(join(p, "\n  "))
</svg>
"""
end
