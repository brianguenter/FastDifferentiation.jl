using FastDifferentiation
using TestItemRunner
using TestItems

@testitem "OptimizationProblems ADNLPProblems" begin
    using Test
    import FastDifferentiation as FD
    using OptimizationProblems
    using OptimizationProblems.ADNLPProblems
    using ADNLPModels
    using NLPModels
    using DifferentiationInterface
    using ADTypes
    using SparseArrays

    """
        test_adnlp_problem(nlp; atol=1e-4, rtol=1e-4, test_di=true, test_sparse=true)

    Test FastDifferentiation derivative computation against an ADNLPModel `nlp`.
    Compares objective gradient, dense/sparse Hessian, and dense/sparse constraint
    Jacobian evaluated at `nlp.meta.x0` against reference values from `NLPModels`.
    Also optionally verifies `DifferentiationInterface` with `AutoFastDifferentiation()`.
    """
    function test_adnlp_problem(nlp; atol=1e-4, rtol=1e-4, test_di=true, test_sparse=true)
        nvar = nlp.meta.nvar
        ncon = nlp.meta.ncon
        x0 = nlp.meta.x0

        vars = FD.make_variables(:x, nvar)
        f_node = nlp.f(vars)

        # 1. FastDifferentiation gradient
        grad_f = vec(FD.jacobian([f_node], vars))
        fg = FD.make_function(grad_f, vars)
        g_val = fg(x0)
        ad_g = grad(nlp, x0)

        if !any(isnan, ad_g) && !any(isinf, ad_g)
            @test isapprox(g_val, ad_g; atol=atol, rtol=rtol)
        end

        # 2. FastDifferentiation Hessian (dense and sparse)
        hess_f = FD.hessian(f_node, vars)
        fh = FD.make_function(hess_f, vars)
        h_val = fh(x0)
        ad_h = hess(nlp, x0)

        if !any(isnan, ad_h) && !any(isinf, ad_h)
            @test isapprox(h_val, ad_h; atol=atol, rtol=rtol)
            if test_sparse
                sp_h = FD.sparse_hessian(f_node, vars)
                f_sph = FD.make_function(sp_h, vars)
                sp_h_val = Matrix(f_sph(x0))
                @test isapprox(sp_h_val, ad_h; atol=atol, rtol=rtol)
            end
        end

        # 3. DifferentiationInterface with AutoFastDifferentiation
        if test_di
            g_di = DifferentiationInterface.gradient(nlp.f, AutoFastDifferentiation(), x0)
            if !any(isnan, ad_g) && !any(isinf, ad_g)
                @test isapprox(g_di, ad_g; atol=atol, rtol=rtol)
            end
            H_di = DifferentiationInterface.hessian(nlp.f, AutoFastDifferentiation(), x0)
            if !any(isnan, ad_h) && !any(isinf, ad_h)
                @test isapprox(H_di, ad_h; atol=atol, rtol=rtol)
            end
        end

        # 4. Constraints (if any)
        if ncon > 0
            cx = fill(FD.Node(0), ncon)
            NLPModels.cons!(nlp, vars, cx)
            jac_c = FD.jacobian(cx, vars)
            fj = FD.make_function(jac_c, vars)
            j_val = fj(x0)
            ad_j = Matrix(jac(nlp, x0))

            if !any(isnan, ad_j) && !any(isinf, ad_j)
                @test isapprox(j_val, ad_j; atol=atol, rtol=rtol)
                if test_sparse
                    sp_j = FD.sparse_jacobian(cx, vars)
                    f_spj = FD.make_function(sp_j, vars)
                    sp_j_val = Matrix(f_spj(x0))
                    @test isapprox(sp_j_val, ad_j; atol=atol, rtol=rtol)
                end
            end
        end

        return true
    end

    # Test the test_adnlp_problem function itself
    @test test_adnlp_problem(OptimizationProblems.ADNLPProblems.rosenbrock())
    @test test_adnlp_problem(OptimizationProblems.ADNLPProblems.hs10())

    # Diverse collection of unconstrained and constrained ADNLPModel problems
    problem_names = [
        # Classical unconstrained optimization problems
        :rosenbrock,
        :beale,
        :booth,
        :brownbs,
        :helical,
        :powellbs,
        :zangwil3,
        # Hock-Schittkowski benchmark problems (with bounds, equality, and inequality constraints)
        :hs1, :hs2, :hs3, :hs4, :hs5, :hs6, :hs7, :hs8, :hs9, :hs10,
        :hs11, :hs12, :hs13, :hs14, :hs15, :hs16, :hs17, :hs18, :hs19, :hs20,
        :hs21, :hs22, :hs23, :hs24, :hs26, :hs27, :hs28, :hs29, :hs30,
        :hs31, :hs32, :hs33, :hs34, :hs35, :hs36, :hs37, :hs38, :hs39, :hs40,
        :hs41, :hs42, :hs43, :hs44, :hs45,
        # Constrained test problem
        :bt1,
    ]

    for name in problem_names
        @testset "$name" begin
            nlp = getfield(OptimizationProblems.ADNLPProblems, name)()
            @test test_adnlp_problem(nlp)
        end
    end
end
