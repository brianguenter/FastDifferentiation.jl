using Documenter, FastDifferentiation

makedocs(
    modules=FastDifferentiation,
    sitename="FastDifferentiation.jl",
    pages=[
        "index.md",
        "Limitations" => "limitations.md",
        "How it works" => "howitworks.md",
        "Examples" => "examples.md",
        "How to use make_function" => "makefunction.md",
        "Benchmarks" => "benchmarks.md",
        "Symbolic processing" => "symbolicprocessing.md",
        "Future work" => "futurework.md",
        "API" => "api.md"
    ],
    warnonly=Documenter.except(:autodocs_block)
)

deploydocs(
    repo="github.com/brianguenter/FastDifferentiation.jl.git",
    devbranch="main"
)