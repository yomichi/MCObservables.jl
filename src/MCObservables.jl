module MCObservables

import Base: show, <<, push!, mean, var, count, isempty
export MCObservable, mean, var, stderror, show, dump_plot

abstract MCObservable

isempty(obs::MCObservable) = count(obs) == 0

function show(io::IO, obs::MCObservable)
  if !isempty(obs)
    print(io, mean(obs), " +/- ", stderror(obs))
  else
    print(io, "No Entries")
  end
end

function dump_plot(io::IO, obs::MCObservable)
  print(io, mean(obs), " ", stderror(obs))
end

include ("observableset.jl")
include ("simple.jl")
include ("binning.jl")
include ("jackknife.jl")
include ("second_jackknife.jl")

## these three definitions are needed to resolve ambiguousness with Base.<<(Any,Integer)
<<(obs::MCObservable, x::Int32) = push!(obs,x)
<<(obs::MCObservable, x::Int64) = push!(obs,x)
<<(obs::MCObservable, x::Integer) = push!(obs,x)

<<(obs::MCObservable, x) = push!(obs,x)

end # module

