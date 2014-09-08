import Base: zero, zeros, deepcopy, mean, show

export SimpleObservable, stddev

type SimpleObservable <: MCObservable
  num :: Int64
  sum :: Float64
  sum2 :: Float64
end

SimpleObservable() = SimpleObservable(0, 0.0, 0.0)
zero(::Type{SimpleObservable}) = SimpleObservable()
zero(o::SimpleObservable) = SimpleObservable()
zeros(::Type{SimpleObservable}, dims...) = reshape([zero(SimpleObservable) for i in 1:prod(dims)],dims)

function reset!(obs :: SimpleObservable)
  obs.num = 0
  obs.sum = 0.0
  obs.sum2 = 0.0
  return obs
end

count(obs::SimpleObservable) = obs.num

function push!(obs :: SimpleObservable, value) 
  obs.num += 1
  obs.sum += value
  obs.sum2 += value^2
  return obs
end

function mean(obs::SimpleObservable)
  if obs.num > 0
    return obs.sum / obs.num
  else
    return NaN
  end
end

function var(obs::SimpleObservable)
  if obs.num  > 1
    v = (obs.sum2 - obs.sum*obs.sum/obs.num)/(obs.num-1)
    return max(v, 0.0)
  elseif obs.num == 1
    return 1.0/0.0
  else
    return NaN
  end
end
stddev(obs::SimpleObservable) = sqrt(var(obs))
stderror(obs::SimpleObservable) = sqrt(var(obs)/count(obs))
function confidence_interval(obs::SimpleObservable, confidence_rate :: Real)
  q = 0.5 + 0.5confidence_rate
  correction = quantile( TDist(obs.num - 1), q)
  serr = stderror(obs)
  return correction * serr
end

function confidence_interval(obs::SimpleObservable, confidence_rate_symbol::Symbol = :sigma1)
  n = parsesigma(confidence_rate_symbol)
  return confidence_interval(obs, erf(0.5n*sqrt(2.0)))
end


export SimpleObservableSet
typealias SimpleObservableSet MCObservableSet{SimpleObservable}
