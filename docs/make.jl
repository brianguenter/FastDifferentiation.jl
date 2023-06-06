using Documenter, FastDifferentiation, ElectronDisplay

makedocs(
    modules=FastDifferentiation,
    sitename="FastDifferentiation.jl",
    pages=[
        "Introduction" => "introduction.md",
        "Limitations" => "limitations.md",
        "How it works" => "howitworks.md",
        "Examples" => "examples.md",
        "Benchmarks" => "benchmarks.md",
        "Symbolic processing" => "symbolicprocessing.md"
    ]
)