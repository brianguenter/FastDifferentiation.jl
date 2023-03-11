average_num_parents(graph) = mean(x -> length(parent_edges(graph, x)), nodes(graph))
export average_num_parents
