function testJackknife(numtest, MCS)
  inform("Start testing Jackknife")

  nsigma = 3
  counts = zeros(Int, nsigma)

  for itest in 1:numtest
    bobs = BinningObservable()
    for mcs in 1:MCS
      bobs << randn()
    end
    jk = jackknife(bobs)
    
    @test_approx_eq mean(bobs) mean(jk)
    @test_approx_eq stderror(bobs) stderror(jk)
  end
  inform("Success testing Jackknife")
end

testJackknife(1_000, 8192)
