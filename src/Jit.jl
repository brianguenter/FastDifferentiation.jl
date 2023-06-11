

is_conditional(a::Node) = is_tree(a) && value(a) in BOOL_OPS

ifelse_nodes(nodes::Vector{Node}) = filter(x -> x === ifelse, nodes) #nodes(gr) is sorted by postorder number for results are as well.

bool_nodes(nodes::Vector{Node}) = filter(x -> x in BOOL_OPS, nodes) #nodes(gr) is sorted by postorder number for results are as well. 


