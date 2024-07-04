import MacroTools
import RuntimeGeneratedFunctions

"""
    extract_statements(ex::Expr)

Extract the statements from an expression.  The expression may represent a function, in
which case the function body alone is extracted.  It may also contain an `@inbounds` block,
in which case the contents of that block are extracted.  The returned quantity is an array
of `Expr`s.
"""
function extract_statements(ex::Expr)
    body = MacroTools.unblock(try
        MacroTools.splitdef(ex)[:body]
    catch
        ex
    end)
    Base.remove_linenums!(body)
    MacroTools.@capture(body, @inbounds begin statements__ end) && return statements
    return body.args
end

"""
    get_returned_expr(statements::Array{Expr})

Extract the expression that is returned by the last statement in the array of expressions.
"""
function get_returned_expr(statements::Array{Expr})
    MacroTools.@capture(statements[end], return returned_expr_) ||
        throw(ArgumentError("No return statement found in last line of statements:\n$statements"))
    return returned_expr
end
