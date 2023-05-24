#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    @inline sigm(x) = tanh(x)

    @inline nlss(V, θ, σ) = 0.5 * (1 + sigm(0.5 * (V - θ) / σ))
    @inline nltc(V, θ, σ, τ0, τ1, τ2) = τ0 + τ1 * (1 - sigm(0.5 * (V - θ) / σ)^2) + τ2 * (1 + sigm(0.5 * (V - θ) / σ))

    function hhvf_4d_ctrl(u, oc, p, sc)
        V, m, h, n = u #voltage, sodium activation, sodium inactivation, potassium activation

        k = oc[1]

        Cm, Isa,
        EL, EK, ENa,
        gL, gK, gNa,
        θm, σm, τm0, τm1, τm2,
        θh, σh, τh0, τh1, τh2,
        θn, σn, τn0, τn1, τn2 = p

        Iapp, data = sc

        dV = (1 / Cm) * (gNa * m * m * m * h * (ENa - V) +
                         gK * n * n * n * n * (EK - V) +
                         gL * (EL - V) +
                         Iapp / Isa) +
             k * (data - V)
        dm = (nlss(V, θm, σm) - m) / nltc(V, θm, σm, τm0, τm1, τm2)
        dh = (nlss(V, θh, σh) - h) / nltc(V, θh, σh, τh0, τh1, τh2)
        dn = (nlss(V, θn, σn) - n) / nltc(V, θn, σn, τn0, τn1, τn2)

        [dV, dm, dh, dn]
    end

    U = make_variables(:u, 4)
    OC = make_variables(:oc, 1)
    P = make_variables(:p, 23)
    SC = make_variables(:sc, 2)

    vf_terms = hhvf_4d_ctrl(U, OC, P, SC)

    jacobian_Expr(vf_terms, vcat(U, OC, P))
end



