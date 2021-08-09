@testset "contains_pi" begin
    # Fully determined cases
    @test HighestCuspedWave.contains_pi(Arb(1), Arb(3)) == (false, false)
    @test HighestCuspedWave.contains_pi(Arb(1), Arb(4)) == (false, true)
    @test HighestCuspedWave.contains_pi(Arb(6), Arb(7)) == (true, false)
    @test HighestCuspedWave.contains_pi(Arb(-1), Arb(1)) == (true, false)
    @test HighestCuspedWave.contains_pi(Arb(-4), Arb(-3)) == (false, true)
    @test HighestCuspedWave.contains_pi(Arb(-7), Arb(-6)) == (true, false)
    @test HighestCuspedWave.contains_pi(Arb(-1), Arb(7)) == (true, true)
    # Indeterminate cases
    @test HighestCuspedWave.contains_pi(Arb(π), Arb(π)) == (false, true)
    @test HighestCuspedWave.contains_pi(Arb(0), Arb(0)) == (true, false)
    @test HighestCuspedWave.contains_pi(2Arb(π), 2Arb(π)) == (true, false)
    @test HighestCuspedWave.contains_pi(Arb(0), Arb(π)) == (true, true)
    @test HighestCuspedWave.contains_pi(-Arb(π), Arb(0)) == (true, true)
    @test HighestCuspedWave.contains_pi(Arb(0), 2Arb(π)) == (true, true)
    @test HighestCuspedWave.contains_pi(Arb(π), 2Arb(π)) == (true, true)
end
