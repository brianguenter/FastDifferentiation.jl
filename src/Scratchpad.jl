#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests

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

    FSD_func = make_function.(dag, Ref([x, y]))
    return FSD_func
    res = [FSD_func[1, 1](1.1, 2.3) FSD_func[1, 2](0.0, 0.0) FSD_func[1, 3](1.1, 0.0)
        FSD_func[2, 1](0.0, 0.0) FSD_func[2, 2](0.0, 0.0) FSD_func[2, 3](0, 2.3)
        FSD_func[3, 1](1.1, 2.3) FSD_func[3, 2](0, 0) FSD_func[3, 3](0, 0)
    ]
    @assert isapprox(res, float_answer)

    return FSD_func
end
export test

#changed
#change
export test
