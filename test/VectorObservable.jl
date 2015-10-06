using MCObservables

function test_VectorObservable(MCS, SObs, VObs)
  inform("Start testing $VObs")

  n = 2

  s = zeros(SObs, n)
  v = VObs()

  for i in 1:MCS
    x = zeros(n)
    for j in 1:n
      x[j] = randn()
      s[j] << x[j]
    end
    v << x
  end

  vjk = jackknife(v)

  for i in 1:n
    @test mean(s[i]) ==  mean(v)[i] 
    @test stderror(s[i]) == stderror(v)[i] 
    sjk = jackknife(s[i])
    @test_approx_eq_eps mean(sjk) mean(vjk)[i] 1e-10
    @test_approx_eq_eps stderror(sjk) stderror(vjk)[i] 1e-10
  end

  inform("Success testing $VObs")
end


result_simpleObservable = test_VectorObservable(8192, SimpleObservable, SimpleVectorObservable)
result_simpleObservable = test_VectorObservable(8192, BinningObservable, BinningVectorObservable)
