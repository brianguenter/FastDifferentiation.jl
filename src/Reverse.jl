function reverse_AD(a::DerivativeGraph, variable_order::AbstractVector{<:Node})
    @assert length(roots(a)) == 1 #only works for Rⁿ->R¹ functions

    let visited = Dict{Int64,Tuple{Int64,Node}}()
        all_vars = Vector{Node}(undef, length(variables(a)))

        function _reverseAD(a::DerivativeGraph, curr_deriv::Node, curr_node::Int64, all_vars, visited)
            if (tmp = get(visited, curr_node, nothing)) === nothing
                visited[curr_node] = (1, curr_deriv)
            else
                visit_count, val = visited[curr_node]
                visited[curr_node] = (visit_count + 1, val + curr_deriv)
            end

            visit_count, val = visited[curr_node]

            if visit_count < length(parent_edges(a, curr_node))
                return
            else
                for c_edge in child_edges(a, curr_node)
                    _reverseAD(a, val * value(c_edge), bott_vertex(c_edge), all_vars, visited)
                end

                if is_variable(node(a, curr_node))
                    all_vars[variable_postorder_to_index(a, curr_node)] = val
                end
            end
        end


        _reverseAD(a, one(Node), root_index_to_postorder_number(a, 1), all_vars, visited)

        result = Vector{Node}(undef, length(variable_order))

        #now map variable values to variable_order
        for (i, node) in pairs(variable_order)
            result[i] = all_vars[variable_node_to_index(a, node)]
        end

        return result
    end
end
export reverse_AD

reverse_AD(a::AbstractVector{<:Node}, variable_order::AbstractVector{<:Node}) = reverse_AD(DerivativeGraph(a), variable_order)
