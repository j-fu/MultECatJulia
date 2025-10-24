using Test
using MultECatJulia
using MySubPackage

MultECatJulia.greet()

@testset "MultECatJulia" begin
    @test true
end
