using FastDifferentiation
using TestItemRunner
using TestItems

@testitem "DifferentiationInterface default scenarios" begin
    using Test
    import FastDifferentiation as FD
    using DifferentiationInterface
    using DifferentiationInterfaceTest
    using ADTypes

    backend = AutoFastDifferentiation()
    scenarios = DifferentiationInterfaceTest.default_scenarios()
    @test length(scenarios) > 0

    test_differentiation(backend, scenarios; logging=false)
end
