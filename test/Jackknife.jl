function testJackknife(numtest, MCS, obs)
  inform("Start testing Jackknife / $obs")

  nsigma = 3
  counts = zeros(Int, nsigma)

  for itest in 1:numtest
    sobs = obs()
    for mcs in 1:MCS
      sobs << randn()
    end
    jk = jackknife(sobs)
    
    @test_approx_eq_eps mean(sobs) mean(jk) 1.0e-10
    @test_approx_eq_eps stderror(sobs) stderror(jk) 1.0e-10
  end
  inform("Success testing Jackknife / $obs")
end

testJackknife(1_000, 8192, SimpleObservable)
testJackknife(1_000, 8192, BinningObservable)
