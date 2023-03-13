function Chebyshev(n, x)
    if n == 0
        y = 1
    elseif n == 1
        y = x
    else
        for i = 1:xleng
            y = 2 * (x) * myChebyshevPoly1(n - 1, x) - myChebyshevPoly1(n - 2, x)
        end
    end
end
