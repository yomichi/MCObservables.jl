function random_samples(rho::Float64, N::Int)
  if !(0.0 <= rho < 1.0)
    throw(DomainError())
  end
  rhobar = sqrt(1.0-rho*rho)
  ret = zeros(N)
  ret[1] = randn()
  for i in 2:N
    ret[i] = ret[i-1]*rho + rhobar * randn()
  end
  return ret
end

function testBinningObservable(numtest, MCS, rho = 0.5)
  inform("Start testing BinningObservable")

  nsigma = 3
  counts = zeros(Int, nsigma)

  progress = 0.05
  for i in 1:numtest
    sobs = SimpleObservable()
    bobs = BinningObservable()
    for x in random_samples(rho, MCS)
      sobs << x
      bobs << x
    end
    m = mean(bobs)
    @test_approx_eq_eps mean(sobs) m 1e-10
    m = abs(m)
    er = stderror(bobs)
    for isigma in 1:nsigma
      if m < isigma * er
        counts[isigma] += 1
      end
    end

    if i/numtest >= progress
      inform(round(Int, progress * 100), "% done.")
      progress += 0.05
    end
  end

  sigma = counts ./ numtest
  
  @show sigma

  ## 1sigma, 2sigma, 3sigma
  ## 0.6827, 0.9545, 0.9973
  @test 0.67 < sigma[1] < 0.69
  @test 0.95 < sigma[2] < 0.96
  @test 0.99 < sigma[3]
  inform("Succeed testing BinningObservable")
end

testBinningObservable(10_000, 16384)
