include("../src/SU2k_sector.jl")
using Random
using Base.Iterators: take

Random.seed!(123456789)

smallset(::Type{I}) where {I<:Sector} = take(values(I), length(values(I)))

@testset "Unitarity SU2_$level F symbols" for level=2:7 begin
    I = SU2kIrrep{level}
        for a in smallset(I), b in smallset(I), c in smallset(I)
            for d in ⊗(a, b, c)
                es = collect(intersect(⊗(a, b), map(dual, ⊗(c, dual(d)))))
                fs = collect(intersect(⊗(b, c), map(dual, ⊗(dual(d), a))))
                @test length(es) == length(fs)
                F = [Fsymbol(a, b, c, d, e, f) for e in es, f in fs]
                @test isapprox(F' * F, one(F); atol=1e-12, rtol=1e-12)
            end
        end
    end
end

@testset "SU2_$level pentagons" for level=2:7 begin
        I = SU2kIrrep{level}
        for a in smallset(I), b in smallset(I), c in smallset(I), d in smallset(I)
            @test TensorKitSectors.pentagon_equation(a, b, c, d; atol=1e-12, rtol=1e-12)
        end
    end
end

@testset "SU2_$level Frobenius Schur" for level=2:12 begin
        I = SU2kIrrep{level}
        for (s, a) in enumerate(smallset(I))
            @test TensorKitSectors.frobenius_schur_phase(a) == (-1)^twice(a.j)
        end
    end
end