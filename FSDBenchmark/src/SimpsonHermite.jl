module SimpsonHermite
##### Set up residual functions
using DifferentialEquations
using Interpolations
using Plots
using FastSymbolicDifferentiation

struct VarAssimVectorField
    vf::Function #out-of-place vector field with type signature vf(u,oc,p,sc) where u is state, oc are optimized controls, p are parameters, sc are static controls (e.g. driving current, observed data)
    D::Int64 #state dim = dim(u)
    OC::Int64 #dim of optimized controls = dim(oc);
    P::Int64 #param count = dim(p)
    SC::Int64 #dim of static static controls = dim(sc);
end

function simpson_hermite_residual(vavf::VarAssimVectorField, x00, oc00, x01, oc01, x10, oc10, p, sc00, sc01, sc10, dt)
    f00 = vavf.vf(x00, oc00, p, sc00)
    f01 = vavf.vf(x01, oc01, p, sc01)
    f10 = vavf.vf(x10, oc10, p, sc10)

    simp = x10 .- x00 .- ((dt / 6) .* (f00 + 4 * f01 + f10))
    herm = x01 .- ((x00 .+ x10) ./ 2 .+ ((dt / 8) .* (f00 - f10)))

    vcat(simp, herm)
end

# Let XOC = [[x00 oc00] [x01 oc01] ... [xN0 ocN0]], P = p, SC = [sc00 sc01 sc10 ... scN0] in these orders.
# DTVEC is a length N vector where DTVEC[j] = the length t_j - t_{j-1} of the jth timestep

function C(vavf::VarAssimVectorField, XOC, P, SC, DTVEC)
    N = div(length(SC) - 1, 2)

    all_constraints = zeros(Node, 2 * vavf.D * N)

    for j = 1:N
        x00 = XOC[2j-1][1:vavf.D]
        oc00 = XOC[2j-1][vavf.D+1:vavf.D+vavf.OC]
        x01 = XOC[2j][1:vavf.D]
        oc01 = XOC[2j][vavf.D+1:vavf.D+vavf.OC]
        x10 = XOC[2j+1][1:vavf.D]
        oc10 = XOC[2j+1][vavf.D+1:vavf.D+vavf.OC]
        sc00 = SC[2j-1]
        sc01 = SC[2j]
        sc10 = SC[2j+1]

        all_constraints[(j-1)*2*vavf.D+1:j*2*vavf.D] .= simpson_hermite_residual(vavf, x00, oc00, x01, oc01, x10, oc10, P, sc00, sc01, sc10, DTVEC[j])
    end

    all_constraints
end


### Need derivatives of C with respect to XOC and P, but not with respect to SC or DTVEC


##### Test with 6D HH model

@inline nlss(V, θ, σ) = 0.5 * (1 + tanh(0.5 * (V - θ) / σ))
@inline nltc(V, θ, σ, τ0, τ1, τ2) = τ0 + τ1 * (1 - tanh(0.5 * (V - θ) / σ)^2) + τ2 * (1 + tanh(0.5 * (V - θ) / σ))

function hhvf_assim_controlled(x, oc, pvec, sc)
    V, m, h, n, z, p = x #voltage, sodium activation, sodium inactivation, potassium activation, M-current activation, persistent sodium activation

    k = oc[1]

    Cm, Isa,
    EL, EK, ENa,
    gL, gK, gNa, gM, gNaP,
    θm, σm, τm0, τm1, τm2,
    θh, σh, τh0, τh1, τh2,
    θn, σn, τn0, τn1, τn2,
    θz, σz, τz0, τz1, τz2,
    θp, σp, τp0, τp1, τp2 = pvec

    Iapp, data = sc

    dV = (1 / Cm) * (gNa * m * m * m * h * (ENa - V)
                     + gK * n * n * n * n * (EK - V)
                     + gNaP * p * (ENa - V)
                     + gM * z * (EK - V)
                     + gL * (EL - V)
                     + Iapp / Isa)
    +k * (data - V)
    dm = (nlss(V, θm, σm) - m) / nltc(V, θm, σm, τm0, τm1, τm2)
    dh = (nlss(V, θh, σh) - h) / nltc(V, θh, σh, τh0, τh1, τh2)
    dn = (nlss(V, θn, σn) - n) / nltc(V, θn, σn, τn0, τn1, τn2)
    dz = (nlss(V, θz, σz) - z) / nltc(V, θz, σz, τz0, τz1, τz2)
    dp = (nlss(V, θp, σp) - p) / nltc(V, θp, σp, τp0, τp1, τp2)

    [dV, dm, dh, dn, dz, dp]
end

vavf = VarAssimVectorField(hhvf_assim_controlled, 6, 1, 35, 2)

### Generate data for testing



## Lorenz applied current input

function lorenz(u, p, t)
    dx = 10.0 * (u[2] - u[1])
    dy = u[1] * (28.0 - u[3]) - u[2]
    dz = u[1] * u[2] - (8 / 3) * u[3]
    [dx, dy, dz]
end

function scale_timeseries(ts, minv, maxv)
    tsmax = maximum(ts)
    tsmin = minimum(ts)
    minv .+ (maxv - minv) * (ts .- tsmin) / (tsmax - tsmin)
end



gtp = (
    Cm=1.0,
    Isa=1.0,
    EL=-55.0,
    EK=-80.0,
    ENa=50.0,
    gL=3.0,
    gK=360.0,
    gNa=1000.0,
    gM=100.0,
    gNaP=30.0,
    θm=-40.0,
    σm=10.0,
    τm0=0.05,
    τm1=0.4,
    τm2=0.0,
    θh=-62.2,
    σh=-5.0,
    τh0=1.0,
    τh1=7.5,
    τh2=0.0,
    θn=-55.0,
    σn=20.0,
    τn0=1.0,
    τn1=5.0,
    τn2=0.0,
    θz=-39.0,
    σz=5.0,
    τz0=75.0,
    τz1=0.0,
    τz2=0.0,
    θp=-47.0,
    σp=3.0,
    τp0=0.05,
    τp1=0.4,
    τp2=0.0
)



##### Test C function
#Test C function by calling like this C(vavf, XOC, P, SC, DTVEC)

function SiH_test()
    ## Generate lorenz data

    # lor_tmax = 160.0
    # lor_dt = 0.5

    lor_tmax = 160.0
    lor_dt = 0.5

    lor_u0 = [-1.31; 0.8; 19.77]
    lor_tspan = (0.0, lor_tmax)
    lor_prob = ODEProblem(lorenz, lor_u0, lor_tspan)
    lor_sol = solve(lor_prob, RK4(), adaptive=true, dt=lor_dt)

    ## Linear interpolation of Lorenz x coordinate for use as HH stimulus
    lorfunc = linear_interpolation(scale_timeseries(lor_sol.t, 0.0, 1200.0),
        scale_timeseries(lor_sol[1, :], -200.0, 50.0),
        extrapolation_bc=0.0)


    function hhvf_for_datagen(x, pvec, t)
        hhvf_assim_controlled(x, [0.0], pvec, [lorfunc(t), 0.0])
    end

    # hh_tmax = 200.0
    # hh_saveat = 0.1

    hh_tmax = 200.0
    hh_saveat = 0.5

    hh_u0 = [-68.24221681836171
        0.056029230048653705
        0.7700232861002139
        0.3402655968929933
        0.0
        0.0]

    hh_prob = ODEProblem(hhvf_for_datagen, hh_u0, (0.0, hh_tmax), collect(gtp))
    @time hh_sol = solve(hh_prob, saveat=hh_saveat)

    ##### Set up inputs to C function

    SC_mat = vcat(lorfunc.(hh_sol.t)', hh_sol[1, :]')
    SC = [SC_mat[:, i] for i in 1:size(SC_mat, 2)]
    DTVEC = diff(hh_sol.t)
    P = collect(gtp)
    XOC_mat = vcat(Array(hh_sol), zeros(length(hh_sol.t))')
    XOC = [XOC_mat[:, i] for i in 1:size(XOC_mat, 2)]


    SC_var = Vector{Vector{Node}}(undef, length(SC))
    XOC_var = Vector{Vector{Node}}(undef, length(XOC))

    for i in eachindex(SC)
        SC_var[i] = make_variables(:SC_var, length(SC[i]))
    end

    for i in eachindex(XOC)
        XOC_var[i] = make_variables(:XOC_var, length(XOC[i]))
    end

    DTVEC_var = make_variables(:DTVEC_var, length(DTVEC))
    P_var = make_variables(:P_var, length(P))

    temp = C(vavf, XOC_var, P_var, SC_var, DTVEC_var)
    @info "Beginning derivative graph construction"
    gr = DerivativeGraph(temp)
    @info "roots = $(length(roots(gr))) variables = $(length(variables(gr))) nodes = $(length(nodes(gr)))"
    @info "Done with derivative graph"
    return gr
    # return nothing
end
export SiH_test

function time_test(gr)
    graph_statistics(gr)
    # Vis.draw_dot(gr, value_labels=false, reachability_labels=true, start_nodes=[54])

    # factor!(gr)
    # symbolic_jacobian!(gr)
    return gr
    return nothing
end
export time_test

end #module


const SiH = SimpsonHermite
export SiH

