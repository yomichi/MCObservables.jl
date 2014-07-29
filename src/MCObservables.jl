module MCObservables

import Base: show, <<, count
export MCObservable, mean, var, stderr, show, dump_plot

abstract MCObservable

function show(io::IO, obs::MCObservable)
  try
    print(io, mean(obs), " +/- ", stderr(obs))
  catch
    print(io, "No Entries")
  end
end

function dump_plot(io::IO, obs::MCObservable)
  print(io, mean(obs), " ", stderr(obs))
end

include ("observableset.jl")
include ("simple.jl")
include ("binning.jl")
include ("jackknife.jl")
include ("second_jackknife.jl")

## these three definitions are needed to resolve ambiguousness with Base.<<(Any,Integer)
<<(obs::MCObservable, x::Int32) = add!(obs,x)
<<(obs::MCObservable, x::Int64) = add!(obs,x)
<<(obs::MCObservable, x::Integer) = add!(obs,x)

<<(obs::MCObservable, x) = add!(obs,x)

end # module

