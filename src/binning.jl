import Base: zeros, mean

export BinningObservable, push!, tau, reset!
export BinningObservableSet

type BinningObservable <: MCObservable
  ## index of these vectors denotes the level of bins (each bin stores the mean of 2^(i-1) values)
  bins :: Vector{Vector{Float64}}  ## Time series of bins
  sum :: Vector{Float64}           ## Summation of stored values
  sum2 :: Vector{Float64}          ## Summation of square of bin's mean
  entries :: Vector{Int}           ## Number of bins
  incompletes :: Vector{Float64}   ## Incomplete bin
  nincomplete :: Vector{Float64}   ## Number of values in incomplete bin
end

BinningObservable() = BinningObservable( (Vector{Float64})[(Float64)[]], zeros(1), zeros(1), zeros(Int,1), zeros(1), zeros(Int,1))

count(b::BinningObservable) = length(b.entries[1])

function push!(b::BinningObservable, x::Real)
  push!(b.bins[1], x)
  b.sum[1] += x
  b.sum2[1] += x*x
  c = b.entries[1]
  b.entries[1] += 1
  level = 2
  binlength = 2
  while c > 0
    if level > length(b.sum)
      push!(b.bins, zeros(0))
      push!(b.sum, 0.0)
      push!(b.sum2, 0.0)
      push!(b.entries, 0)
      push!(b.incompletes, b.sum[1])
      push!(b.nincomplete, b.entries[1])
    else
      b.incompletes[level] += x
      b.nincomplete[level] += 1
    end
    if b.nincomplete[level] == binlength
      ## last bin just filled
      filled = b.incompletes[level] / binlength
      push!(b.bins[level], filled)
      b.sum[level] += filled
      b.sum2[level] += filled*filled
      b.entries[level] += 1
      b.incompletes[level] = 0.0
      b.nincomplete[level] = 0
    end
    c >>= 1
    level += 1
    binlength <<= 1
  end
end

#=
function push!(b::BinningObservable, x)
  push!(b.bins[1], x)
  b.sum[1] += x
  b.sum2[1] += x*x

  i = b.count
  b.count += 1
  b.entries[1] += 1
  blen = 1
  bin = 1

  # push new filled bins
  while i&1 == 1
    blen *= 2
    bin += 1
    if bin > length(b.sum)
      push!(b.bins, zeros(0))
      push!(b.sum, 0.0)
      push!(b.sum2, 0.0)
      push!(b.entries, 0)
    end

    x1 =  b.sum[1] - b.sum[bin]  ## summations of not stored values (sum of new bin)
    x1 /= blen  ## mean value of new bin
    push!(b.bins[bin], x1)
    y1 = x1*x1

    b.sum2[bin] += y1
    b.sum[bin] = b.sum[1]  ## == b.sum[bin] += x1
    b.entries[bin] += 1

    i >>= 1
  end
end
=#

maxlevel(b::BinningObservable) = length(b.sum) > 7 ? length(b.sum)-7 : 1

function binvariance_detail(b::BinningObservable, bindex::Int)
  binlength = 1<<(bindex-1)
  entries = b.entries[bindex]
  s = b.sum[bindex]
  s2 = b.sum2[bindex]
  nin = b.nincomplete[bindex]
  if nin > 0
    incomplete = b.incompletes[bindex] / nin
    s += incomplete
    s2 += incomplete*incomplete
    entries += 1
  end
  return (s2 - s*s/entries)/(entries-1), entries
end

binvariance(b::BinningObservable, bindex::Int = maxlevel(b)) = binvariance_detail(b,bindex)[1]

function mean(b::BinningObservable)
  if b.entries[1] == 0
    throw(DomainError())
  else
    return b.sum[1]/b.entries[1]
  end
end

function var(b::BinningObservable)
  if b.entries[1] > 1
    return (b.sum2[1] - b.sum[1]*b.sum[1]/b.entries[1])/(b.entries[1]-1)
  elseif b.count == 1
    return 1.0/0.0
  else
    throw(DomainError())
  end
end

function stderr(b::BinningObservable, bindex::Int = maxlevel(b))
  bvar,entries = binvariance_detail(b,bindex)
  return sqrt(bvar/entries)
end

function tau(b::BinningObservable, bindex::Int = maxlevel(b))
  bvar,bentries = binvariance_detail(b,bindex)
  return 0.5*( (bvar*b.entries[1])/(var(b)*bentries) - 1.0)
end

function reset!(b::BinningObservable)
  b.bins = (Vector{Float64})[]
  b.sum = zeros(1)
  b.sum2 = zeros(1)
  b.entries = zeros(Int, 1)
  b.count = 0
  return b
end

function show(io::IO, obs::BinningObservable)
  try
    print(io, mean(obs), " +/- ", stderr(obs), "; tau = ", tau(obs))
  catch
    print(io, "No entries")
  end
end


typealias BinningObservableSet MCObservableSet{BinningObservable}

