using MCObservables

function test_simpleObservable(numtest, MCS)
  inform("Start testing SimpleObservable")

  nsigma = 3

  counts = zeros(Int, nsigma)

  for i in 1:numtest
    obs = SimpleObservable()
    for j in 1:MCS
      obs << randn()
    end
    m = abs(mean(obs))
    er = stderror(obs)
    for isigma in 1:nsigma
      if m < isigma*er
        counts[isigma] += 1
      end
    end
  end
  sigma = counts ./ numtest
  @show sigma

  @test 0.68 < sigma[1] < 0.69
  @test 0.95 < sigma[2] < 0.96
  @test 0.99 < sigma[3]
  inform("Success testing SimpleObservable")
end


result_simpleObservable = test_simpleObservable(100_000, 8192)
