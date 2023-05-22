@@ -1,70 +0,0 @@
function_array(a::Symbol, n::Int, var::T) where {T<:Node} = [function_of(Symbol(a, x), var) for x in 1:n]
export function_array

random_transform(joint_function::Node) = transformation(matrix(random_rotation(joint_function)), SVector{3}(Node.(rand(3))))
export random_transform

random_inertia() = transformation(Node.(collect(LinearAlgebra.Symmetric(rand(3, 3)))))
export random_inertia
struct Linkage
    Aᵢ::Vector{Matrix{Node}}
    joint_angles::Vector{Node}
    inertia_tensors::Vector{Matrix{Node}}
    time_var::Node
    mass::Vector{Float64}

    function Linkage(n::Int)
        Symbolics.@variables t
        tvar = Node(t)
        joint_angles = function_array(:q, n, tvar)
        inertia = [random_inertia() for _ in 1:n]
        new(random_transform.(joint_angles), joint_angles, inertia, tvar, rand(n))
    end
end
export Linkage

# variables(linkage::Linkage) = reduce(∪, variables.(linkage.Aᵢ))
# export variables
num_links(linkage::Linkage) = length(linkage.Aᵢ)
Jⱼ(linkage::Linkage, index::Int) = linkage.inertia_tensors[index]
export Jⱼ
mⱼ(linkage::Linkage, index::Int) = linkage.mass[index]
export mⱼ
rⱼ(linkage::Linkage, index::Int) = W(linkage, index)[1:3, 4]
export rⱼ

W(linkage::Linkage, index::Int) = Node.(reduce(*, linkage.Aᵢ[1:index], init=I()))
export W

function τᵢ(linkage::Linkage, index::Int)
    qᵢ = linkage.joint_angles[index]
    sum = Node(0)
    t = linkage.time_var
    g = rand(3) #random gravity vector


    for j in index:num_links(linkage)
        DWⱼ_qᵢ = derivative(W(linkage, j), qᵢ)
        DWⱼ_qᵢ[4, 4] = Node(1.0) #hack to make sure still homogeneous transformation
        Dtt = derivative(W(linkage, j), t, t)
        J = Jⱼ(linkage, index)
        grav = mⱼ(linkage, index) * (transpose(g) * rⱼ(linkage, index))
        trace = LinearAlgebra.tr(Node.(DWⱼ_qᵢ * J * Dtt)) #for some reason matrix multiplication can't determine the type of the output. Maybe I have to define an interface function.
        sum += trace - grav
    end

    return sum
end
export τᵢ

function lagrangian_dynamics(n::Integer)
    result = Node[]
    links = Linkage(n)
    for i in eachindex(links.Aᵢ)
        torque = τᵢ(links, i)

        push!(result, torque)
    end
    return result
end
export lagrangian_dynamics