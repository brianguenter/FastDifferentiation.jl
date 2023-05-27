#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 
using StaticArrays
using FiniteDifferences
using .FSDTests

using Interpolations


struct SimApp
    # variables
    nx::Int
    Δt::Float64
    Δx::Float64

    # coordinates
    x::AbstractRange
end
SimApp(nx, Δt) = SimApp(nx, Δt, 2π / nx, range(0.0, 2π * (1 - 1 / nx), nx))

function step(y, app::SimApp)
    # initialize interpolation
    itp = interpolate(y, BSpline(Linear(Periodic())))
    sitp = scale(itp, app.x)
    extp = extrapolate(sitp, Periodic())

    # interpolate solution at departure points
    x_back = app.x .- app.Δt * y
    y_intp = extp.(x_back)

    return y_intp
end

function test()
    @variables u[1:64]

    nu = [Node(u_i) for u_i ∈ u]

    app = SimApp(64, 0.02)

    fs = step(nu, app)


    fs = step(collect(u), app)
end



