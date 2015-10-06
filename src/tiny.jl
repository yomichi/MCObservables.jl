import Base: zero, zeros, deepcopy, mean, show

export TinyObservable, stddev

type TinyObservable <: ScalarObservable
  num :: Int64
  sum :: Float64
  sum2 :: Float64
end

TinyObservable() = TinyObservable(0, 0.0, 0.0)
zero(::Type{TinyObservable}) = TinyObservable()
zero(o::TinyObservable) = TinyObservable()
zeros(::Type{TinyObservable}, dims...) = reshape([zero(TinyObservable) for i in 1:prod(dims)],dims)

function reset!(obs :: TinyObservable)
  obs.num = 0
  obs.sum = 0.0
  obs.sum2 = 0.0
  return obs
end

count(obs::TinyObservable) = obs.num

function push!(obs :: TinyObservable, value) 
  obs.num += 1
  obs.sum += value
  obs.sum2 += value^2
  return obs
end

function mean(obs::TinyObservable)
  if obs.num > 0
    return obs.sum / obs.num
  else
    return NaN
  end
end

function var(obs::TinyObservable)
  if obs.num  > 1
    v = (obs.sum2 - obs.sum*obs.sum/obs.num)/(obs.num-1)
    return maxzero(v)
  elseif obs.num == 1
    return Inf
  else
    return NaN
  end
end
stddev(obs::TinyObservable) = sqrt(var(obs))
stderror(obs::TinyObservable) = sqrt(var(obs)/count(obs))
function confidence_interval(obs::TinyObservable, confidence_rate :: Real)
  if count(obs) < 2
    return Inf
  end
  q = 0.5 + 0.5confidence_rate
  correction = quantile( TDist(obs.num - 1), q)
  serr = stderror(obs)
  return correction * serr
end

function confidence_interval(obs::TinyObservable, confidence_rate_symbol::Symbol = :sigma1)
  n = parsesigma(confidence_rate_symbol)
  return confidence_interval(obs, erf(0.5n*sqrt(2.0)))
end

function merge!(obs::TinyObservable, other::TinyObservable)
  obs.num += other.num
  obs.sum += other.sum
  obs.sum2 += other.sum2
  return obs
end
merge(lhs::TinyObservable, rhs::TinyObservable) = merge!(deepcopy(lhs), rhs)

export TinyObservableSet
typealias TinyObservableSet MCObservableSet{TinyObservable}

function merge!(obs::TinyObservableSet, other::TinyObservableSet)
  obs_names = Set(keys(obs))
  union!(obs_names, Set(keys(other)))
  for name in obs_names
    if !haskey(obs, name)
      # in other only
      obs[name] = deepcopy(other[name])
    elseif haskey(other, name)
      # in both
      merge!(obs[name], other[name])

    # else
    # in obs only
    # NOTHING to do
    end
  end
  return obs
end
merge(lhs::TinyObservableSet, rhs::TinyObservableSet) = merge!(deepcopy(lhs), rhs)
