import Base: zeros, mean

export BinningObservable, add!, tau, reset!
export BinningObservableSet

type BinningObservable <: MCObservable
  ## index of these vectors denotes the level of bins (each bin stores 2^(i-1) values)
  bins :: Vector{Vector{Float64}}  ## Bins
  sum :: Vector{Float64}           ## Summations of stored values
  sum2 :: Vector{Float64}          ## Summations of square of bin's mean
  entries :: Vector{Int}           ## Numbers of bin 
  count :: Int                     ## Number of stored values
end

BinningObservable() = BinningObservable( (Vector{Float64})[(Float64)[]], zeros(1), zeros(1), zeros(Int,1), 0)

function add!(b::BinningObservable, x)
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

maxlevel(b::BinningObservable) = length(b.sum) > 7 ? length(b.sum)-7 : 1

function binmean(b::BinningObservable, bindex::Int = maxlevel(b))
  return b.sum[bindex] / (b.entries[bindex] * 1<<(bindex-1))
end

function binvariance(b::BinningObservable, bindex::Int = maxlevel(b))
  m = b.sum[bindex] / (1<<(bindex-1))  ## summation of bin's mean
  return (b.sum2[bindex] - m*m/b.entries[bindex])/(b.entries[bindex]-1)
end

function mean(b::BinningObservable)
  if b.count == 0
    throw(DomainError())
  else
    return b.sum[1]/b.count
  end
end

function var(b::BinningObservable)
  if b.count > 1
    return (b.sum2[1] - b.sum[1]*b.sum[1]/b.count)/(b.count-1)
  elseif b.count == 1
    return 1.0/0.0
  else
    throw(DomainError())
  end
end

function stderr(b::BinningObservable, bindex::Int = maxlevel(b))
  return sqrt(binvariance(b,bindex)/b.entries[bindex])
end

function tau(b::BinningObservable, bindex::Int = maxlevel(b))
  return 0.5*((binvariance(b,bindex)*b.entries[1])/(binvariance(b,1)*b.entries[bindex]) - 1.0)
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

