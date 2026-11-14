using FastDifferentiation
using TestItemRunner
using TestItems

# Compares `jacobian` with exact forward-mode derivatives on random expression
# graphs (see test/random_graphs.jl). The seed is fixed so CI runs are
# reproducible; scripts/fuzz_factoring.jl runs the same checks for longer with
# any seed. Before the factoring fix, about 4% of the tiny graphs, 30% of the
# small ones and nearly all of the medium ones failed.
@testitem "random graphs" begin
    using Test
    using Random
    include(joinpath(@__DIR__, "random_graphs.jl"))
    using .RandomGraphs

    rng = MersenneTwister(20261114)
    for (label, n, nvars, nops, nroots) in [
        ("tiny", 2000, 2:3, 3:6, 1:3),
        ("small", 1000, 2:5, 8:15, 1:4),
        ("medium", 200, 3:8, 20:40, 2:6),
        ("large", 20, 4:12, 60:120, 3:10),
    ]
        failures, _ = RandomGraphs.fuzz(rng, n; nvars, nops, nroots)
        for (spec, status, msg) in first(failures, 3)
            @error "random graph $label: $status" code = spec_code(spec) msg
        end
        @test isempty(failures)
    end
end
