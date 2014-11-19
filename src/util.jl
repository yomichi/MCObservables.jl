squared(x::Real) = x*x
squared(x::Vector) = x.*x
max0(x::Real) = max(x,zero(x))
