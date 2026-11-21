# Turns the raw suite result lines in data/ into data/results.json and rebuilds ../analysis.html.
#   julia investigations/tracingjit/report/parse.jl
# Later files override earlier ones when a case appears twice (the rerun file corrects the first run of cases 78-87).
const HERE = @__DIR__
files = [joinpath(HERE, "data", f) for f in ("suiteA.txt", "suiteB_68.txt", "suiteB_87.txt", "suiteB_109.txt", "suiteB_127.txt", "suiteB_rerun_78_87.txt")]
rx = r"^(\S+)\s+(.*?)\s+(\S+)\s+err\s+(\S+)\s+ref\s+(\S+)\s+guards\s+(\d+)\s+trace\s+(\d+)\s+jac\s+(\d+)\s+\[(\w+)\s*\]\s+first\s+([\d.]+)s\s+hit\s+(\S+)us\s+FD\s+(\S+)us\s+\((\S+)\)\s+hits/misses\s+(\d+)/(\d+)\s+(\w+)\s+extra\[(.*?)\]\s*(.*)$"
rows = Dict{String,String}(); order = String[]
num(s) = (v = tryparse(Float64, s); v === nothing || !isfinite(v) ? "null" : string(v))
js(s) = "\"" * replace(replace(s, "\\" => "\\\\"), "\"" => "\\\"") * "\""
for f in files, line in eachline(f)
    m = match(rx, line)
    if m === nothing; println(stderr, "not parsed: ", first(line, 70)); continue; end
    g = m.captures
    name = String(g[2]); name in order || push!(order, name)
    rows[name] = "{\"status\":$(js(g[1])),\"name\":$(js(name)),\"pkg\":$(js(g[3])),\"err\":$(num(g[4])),\"ref\":$(js(g[5])),\"guards\":$(g[6]),\"trace\":$(g[7]),\"jac\":$(g[8]),\"method\":$(js(g[9])),\"first\":$(g[10]),\"hit\":$(num(g[11])),\"fd\":$(num(g[12])),\"speed\":$(num(replace(g[13], "x" => ""))),\"hits\":$(g[14]),\"misses\":$(g[15]),\"mode\":$(js(g[16])),\"extra\":$(js(g[17])),\"msg\":$(js(first(String(g[18]), 150)))}"
end
json = "[" * join((rows[n] for n in order), ",\n") * "]\n"
write(joinpath(HERE, "data", "results.json"), json)
println(length(order), " cases written to data/results.json")

# page = head + data + script (a "</" inside the data is escaped so it cannot close the script element)
page = read(joinpath(HERE, "head.html"), String) * "<script>const DATA=\n" * replace(json, "</" => "<" * string(Char(92)) * "/") * ";</script>\n" * read(joinpath(HERE, "tail.html"), String)
write(joinpath(HERE, "..", "analysis.html"), page)
println("wrote analysis.html (", sizeof(page), " bytes)")
