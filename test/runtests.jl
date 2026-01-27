using QWignerSymbols
using Test, TestExtras
using HalfIntegers

@testset "QWignerSymbols.jl" verbose = true begin
    @testset "q_clebsch_gordan" begin
        include("q_clebsch_gordan.jl")
    end

    @testset "q_sector" begin
        include("q_sector.jl")
    end
end
