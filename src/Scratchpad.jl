#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using Profile
using PProf
using .TestCases

function test()
    Symbolics.@variables x y

    A = [x^2+y 0 2x
        0 0 2y
        y^2+x 0 0]
    dag = expr_to_dag.(A)
    symbolics_answer = Symbolics.substitute.(A, Ref(Dict(x => 1.1, y => 2.3)))
    float_answer = similar(symbolics_answer, Float64)
    for index in eachindex(symbolics_answer)
        float_answer[index] = symbolics_answer[index].val
    end

    FSD_func = make_function.(dag)
    res = [FSD_func[1, 1](1.1, 2.3) FSD_func[1, 2]() FSD_func[1, 3](1.1)
        FSD_func[2, 1]() FSD_func[2, 2]() FSD_func[2, 3](2.3)
        FSD_func[3, 1](2.3, 1.1) FSD_func[3, 2]() FSD_func[3, 3]()
    ]
    @assert isapprox(res, float_answer)

end
export test

function profile()
    graph, qx, qy, qz = to_graph(25)
    Profile.clear()
    @profile symbolic_jacobian!(graph, [Node(qx), Node(qy), Node(qz)])
    pprof()
end
export profile

#changed
#change
export test
