using FastDifferentiation
using TestItemRunner
using TestItems

@testitem "complex expressions" begin
    using Test
    import FiniteDifferences

    # issue 23
    @variables x02 x04 u1

    vars = [x02, x04, u1]

    c1 = -(cos((x02 + (x04 + ((-(cos(x02)) * (((-(x04) * sin(x02)) * x04) - u1)) + sin(x02))))))
    c2 = ((-((x04 + (-(cos((x02 + x04))) * (((-((x04 + (-(cos(x02)) * (((-(x04) * sin(x02)) * x04) - u1)))) * sin((x02 + x04))) * x04) - u1)))) * sin(x02)) * x04)

    j1 = jacobian([c1], vars)
    j2 = jacobian([c2], vars)
    j_sum = jacobian([c1+c2], vars)

    f_j1 = make_function(j1, vars)
    f_j2 = make_function(j2, vars)
    f_sum = make_function(j_sum, vars)

    eval_c1_c2(vals) = let (x02_v, x04_v, u1_v) = (vals[1], vals[2], vals[3])
        v_c1 = -(cos((x02_v + (x04_v + ((-(cos(x02_v)) * (((-(x04_v) * sin(x02_v)) * x04_v) - u1_v)) + sin(x02_v))))))
        v_c2 = ((-((x04_v + (-(cos((x02_v + x04_v))) * (((-((x04_v + (-(cos(x02_v)) * (((-(x04_v) * sin(x02_v)) * x04_v) - u1_v)))) * sin((x02_v + x04_v))) * x04_v) - u1_v)))) * sin(x02_v)) * x04_v)
        [v_c1 + v_c2]
    end

    test_pts_23 = [[0.3, 0.7, 0.2], [1.1, -0.5, 0.9], [-0.8, 0.4, -0.3]]
    for pt in test_pts_23
        @test isapprox(f_sum(pt), f_j1(pt) + f_j2(pt), atol=1e-12)
        fd_grad = FiniteDifferences.jacobian(FiniteDifferences.central_fdm(5, 1), eval_c1_c2, pt)[1]
        @test isapprox(f_sum(pt), fd_grad, atol=1e-6)
    end


    # issue 108
    𝐗 = make_variables(:𝐗, 2)
    𝐘 = make_variables(:𝐘, 3)
    eqs = FastDifferentiation.Node[ 𝐘[1] * (𝐘[2] + 𝐘[3]) * 𝐘[1] * (1 / 𝐗[1]) * 𝐘[3] * 𝐗[2] ^ 𝐘[1]
                                    𝐘[1] * (𝐘[2] + 𝐘[3]) * 𝐘[1] + (𝐘[2] + 𝐘[3])
                                    𝐘[1] * (𝐘[2] + 𝐘[3]) * 𝐘[1] * 𝐘[3] * 𝐗[2] ^ 𝐘[1] ]
    jac108 = jacobian(eqs, 𝐘)
    f_jac108 = make_function(jac108, [𝐗..., 𝐘...])

    eval_eqs108(X_vals, Y_vals) = [
        Y_vals[1] * (Y_vals[2] + Y_vals[3]) * Y_vals[1] * (1 / X_vals[1]) * Y_vals[3] * X_vals[2] ^ Y_vals[1],
        Y_vals[1] * (Y_vals[2] + Y_vals[3]) * Y_vals[1] + (Y_vals[2] + Y_vals[3]),
        Y_vals[1] * (Y_vals[2] + Y_vals[3]) * Y_vals[1] * Y_vals[3] * X_vals[2] ^ Y_vals[1]
    ]

    test_pts_108 = [([1.5, 2.0], [0.8, 1.2, 0.5]), ([2.5, 1.1], [1.1, 0.4, 0.9])]
    for (xv, yv) in test_pts_108
        eval_y(y) = eval_eqs108(xv, y)
        fd_jac = FiniteDifferences.jacobian(FiniteDifferences.central_fdm(5, 1), eval_y, yv)[1]
        @test isapprox(f_jac108([xv..., yv...]), fd_jac, atol=1e-6)
    end


    # my problem
    function logpdf_mvn2x2(x1, x2, s1, s2, ρ)
        det = s1^2*s2^2*(1-ρ^2)
        -log(2π) - 0.5*log(det) - 0.5*(x1^2*s2^2 + x2^2*s1^2 - 2*ρ*x1*x2*(s1*s2))/det
    end
    function logpdfn(x, s)
        -0.5*log(2π) -log(s) - 0.5*x^2/s^2
    end
    function _sharell(σϵ, σζ, ρ, ψ, σν, α, logBm, logBl,
                  sm, sl, m, ℓ)
        ϵ = (-sm + logBm + 0.5*σϵ*σϵ)
        ζ = sl #(sl - sm - logBl + logBm + 0.5*σζ*σζ)
        ν = ℓ #-(ℓ - m + α + logBm - logBl + ψ + σζ*(σζ/2))
        ll = logpdf_mvn2x2(ϵ, ζ, σϵ, σζ, ρ) #+ logpdfn(ν, σν)
    end

    @variables σϵ σζ ρ ψ σν α logBm logBl sm sl m ℓ
    vars_sharell = [σϵ, σζ, ρ, ψ, σν, α, logBm, logBl]
    f = _sharell(σϵ, σζ, ρ, ψ, σν, α, logBm, logBl,
                  sm, sl, m, ℓ)
    ∇²f = hessian(f, vars_sharell)
    f_hess = make_function(∇²f, [vars_sharell..., sm, sl, m, ℓ])

    eval_f(v, extra) = _sharell(v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], extra[1], extra[2], extra[3], extra[4])

    test_pts_sharell = [
        ([0.8, 0.7, 0.3, 0.1, 0.5, 0.2, 0.4, 0.3], [0.1, 0.2, 0.5, 0.3]),
        ([1.2, 0.9, -0.4, 0.3, 0.6, 0.1, 0.5, 0.2], [0.3, 0.1, 0.4, 0.2])
    ]
    for (v, extra) in test_pts_sharell
        fd_hess = FiniteDifferences.jacobian(FiniteDifferences.central_fdm(5, 1), x -> FiniteDifferences.grad(FiniteDifferences.central_fdm(5, 1), y -> eval_f(y, extra), x)[1], v)[1]
        @test isapprox(f_hess([v..., extra...]), fd_hess, atol=1e-5)
    end

    # another one
    @variables a1 a2 s1 s2 rho eps y

    mgf = 0.5 * s2^2 + log(1 + 2*a1*s2 + (a1^2 + 2*a2)*(s2^2 + 1) + 2*a1*a2*(s2^3 + 3*s2) + a2^2*(s2^4 + 6*s2^2 + 3)) - log(1 + a1^2 + 2*a2 + 3*a2^2)
    ζ = y + mgf
    z = (ζ / s2 - rho * eps / s1) / sqrt(1 - rho^2 + 1e-12)
    f = log((1 + a1 * z + a2 * z^2)^2)

    ∇²f = hessian(f, [a1, a2, s1, s2, rho, eps, y])
    f_hess = make_function(∇²f, [a1, a2, s1, s2, rho, eps, y])

    eval_last(v) = begin
        a1v, a2v, s1v, s2v, rhov, epsv, yv = v
        mgfv = 0.5 * s2v^2 +
                log(1 + 2 * a1v * s2v + (a1v^2 + 2 * a2v) * (s2v^2 + 1) +
                    2 * a1v * a2v * (s2v^3 + 3 * s2v) +
                    a2v^2 * (s2v^4 + 6 * s2v^2 + 3)) -
                log(1 + a1v^2 + 2 * a2v + 3 * a2v^2)
        ζv = yv + mgfv
        zv = (ζv / s2v - rhov * epsv / s1v) / sqrt(1 - rhov^2 + 1e-12)
        log((1 + a1v * zv + a2v * zv^2)^2)
    end

    test_pts_last = [
        [0.8, 0.2, 1.1, 0.9, 0.0, 0.3, 0.4],
        [1.1, -0.1, 0.7, 1.3, -0.3, 0.2, -0.5],
    ]
    for v in test_pts_last
        fd_hess = FiniteDifferences.jacobian(
            FiniteDifferences.central_fdm(5, 1),
            x -> FiniteDifferences.grad(FiniteDifferences.central_fdm(5, 1), eval_last, x)[1],
            v,
        )[1]
        @test isapprox(f_hess(v), fd_hess, atol=1e-5)
    end


end
