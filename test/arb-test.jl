@testset "contains_pi" begin
    RR = RealField(64)
    # Fully determined cases
    @test HighestCuspedWave.contains_pi(RR(1), RR(3)) == (false, false)
    @test HighestCuspedWave.contains_pi(RR(1), RR(4)) == (false, true)
    @test HighestCuspedWave.contains_pi(RR(6), RR(7)) == (true, false)
    @test HighestCuspedWave.contains_pi(RR(-1), RR(1)) == (true, false)
    @test HighestCuspedWave.contains_pi(RR(-4), RR(-3)) == (false, true)
    @test HighestCuspedWave.contains_pi(RR(-7), RR(-6)) == (true, false)
    @test HighestCuspedWave.contains_pi(RR(-1), RR(7)) == (true, true)
    # Indeterminate cases
    @test HighestCuspedWave.contains_pi(RR(π), RR(π)) == (false, true)
    @test HighestCuspedWave.contains_pi(RR(0), RR(0)) == (true, false)
    @test HighestCuspedWave.contains_pi(2RR(π), 2RR(π)) == (true, false)
    @test HighestCuspedWave.contains_pi(RR(0), RR(π)) == (true, true)
    @test HighestCuspedWave.contains_pi(-RR(π), RR(0)) == (true, true)
    @test HighestCuspedWave.contains_pi(RR(0), 2RR(π)) == (true, true)
    @test HighestCuspedWave.contains_pi(RR(π), 2RR(π)) == (true, true)
end
