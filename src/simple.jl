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

function add!(obs :: SimpleObservable, value) 
  obs.num += 1
  obs.sum += value
  obs.sum2 += value^2
  return obs
end

function mean(obs::SimpleObservable)
  if obs.num > 0
    return obs.sum / obs.num
  else
    throw(DomainError())
  end
end

function var(obs::SimpleObservable)
  if obs.num  > 1
    return obs.sum2/(obs.num-1) - obs.sum*obs.sum/(obs.num*(obs.num-1))
  elseif obs.num == 1
    return 1.0/0.0
  else
    throw(DomainError())
  end
end
stddev(obs::SimpleObservable) = sqrt(var(obs))
stderr(obs::SimpleObservable) = sqrt(var(obs)/obs.num)


export SimpleObservableSet
typealias SimpleObservableSet MCObservableSet{SimpleObservable}
