export Jackknife, JackknifeSet
export jackknife

type Jackknife <: MCObservable
  xs :: Vector{Float64}
end

Jackknife(jk::Jackknife, f::Function) = Jackknife(map(f,jk.xs))

function Jackknife(b::BinningObservable)
  if count(b) > 0
    xs = b.bins
    s = sum(xs)
    n = length(xs)
    jk = map(x->(s-x)/(n-1), xs)
    return Jackknife(jk)
  else
    return Jackknife(zeros(0))
  end
end

count(jk::Jackknife) = length(jk.xs)

function mean(jk::Jackknife)
  if isempty(jk) 
    return nan(Float64)
  else
    return mean(jk.xs)
  end
end
function stderror(jk::Jackknife)
  n = count(jk)
  if n == 0
    return nan(Float64)
  elseif n == 1
    return inf(Float64)
  else
    sums = mapreduce(x->[x,x*x], +, jk.xs)
    sums /= n
    sigma2 = sums[2] - sums[1]*sums[1]
    sigma2 *= n-1
    sigma2 = max(sigma2, 0.0)
    return sqrt(sigma2)
  end
end


unary_functions = (
  :-,
  :sin, :cos, :tan,
  :sind, :cosd, :tand,
  :sinpi, :cospi,
  :sinh, :cosh, :tanh,
  :asin, :acos, :atan,
  :asind, :acosd, :atand,
  :sec, :csc, :cot,
  :secd, :cscd, :cotd,
  :asec, :acsc, :acot,
  :asecd, :acscd, :acotd,
  :sech, :csch, :coth,
  :asinh, :acosh, :atanh,
  :asech, :acsch, :acoth,
  :sinc, :cosc,
  :log, :log2, :log10, :log1p,
  :exp, :exp2, :exp10, :expm1,
  :abs, :abs2,
  :sqrt, :cbrt,
  :erf, :erfc, :erfcx,
  :erfinv, :erfcinv,
  :gamma, :lgamma, :lfact,
  :digamma, :invdigamma, :trigamma,
  :airyai, :airyprime, :airyaiprime,
  :airybi, :airybiprime,
  :besselj0, :besselj1, 
  :bessely0, :bessely1,
  :eta, :zeta
)

for op in unary_functions
  eval( Expr(:import, :Base, op) )
  eval( Expr(:export, op) )
  @eval ($op)(jk::Jackknife) = Jackknife(jk, $op)
end

binary_functions = (
  :+, :-, :*, :/, :\
)

for op in binary_functions
  eval( Expr(:import, :Base, op) )
  eval( Expr(:export, op) )
  @eval ($op)(jk::Jackknife, rhs::Real) = Jackknife(jk, lhs->($op)(lhs,rhs))
  @eval ($op)(lhs::Real, jk::Jackknife) = Jackknife(jk, rhs->($op)(lhs,rhs))
  op_bw = symbol("."*string(op))
  @eval ($op)(lhs::Jackknife, rhs::Jackknife) = Jackknife( ($op_bw) (lhs.xs, rhs.xs))
end

typealias JackknifeSet MCObservableSet{Jackknife}

jackknife(bin::BinningObservable) = Jackknife(bin)
function jackknife(obsset :: BinningObservableSet)
  JK = JackknifeSet()
  for (k,v) in obsset
    JK[k] = Jackknife(v)
  end
  return JK
end

