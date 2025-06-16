using HalfIntegers
using TensorKit
using QWignerSymbols
using Base: HasLength

struct SU2kIrrep{k} <: TensorKit.Sector
    j::HalfInt
    function SU2kIrrep{k}(j) where {k}
        j >= zero(j) || error("j=$j is not a positive half-integer")
        j <= half(k) || error("j=$j can be at most k/2=$(half(k))")
        return new{k}(j)
    end
end

function  SU2kIrrep(j, k::Int)
    k >= 1 || error("Level k must be positve")
    return SU2kIrrep{k}(j)
end
k(::SU2kIrrep{K}) where {K} = K
q(s::SU2kIrrep{K}) where {K} = Phase(1/(k(s) + 2)) # exp(2*π*im/(k(s)+2))
Base.convert(T::Type{<:SU2kIrrep}, j::Real) = T(j)

Base.one(::Type{T}) where {T<:SU2kIrrep} = T(zero(HalfInt))
Base.conj(s::SU2kIrrep) = s
Base.isless(s1::T, s2::T) where {T<:SU2kIrrep} = isless(s1.j, s2.j)

function TensorKit.:⊗(s1::T, s2::T) where {T<:SU2kIrrep}
    return TensorKit.SectorSet{T}(abs(s1.j - s2.j):min(s1.j + s2.j, k(s1) - s1.j - s2.j))
end
Base.IteratorSize(::Type{<:TensorKit.SectorValues{<:SU2kIrrep}}) = HasLength()
Base.length(::TensorKit.SectorValues{SU2kIrrep{k}}) where {k} = k + 1
function Base.iterate(::TensorKit.SectorValues{SU2kIrrep{k}}, i=0) where {k}
    return (SU2kIrrep{k}(half(i)), i + 1)
end
function Base.getindex(::TensorKit.SectorValues{SU2kIrrep{k}}, i::Int) where {k}
    if i <= k + 1
        return SU2kIrrep(half(i - 1), k)
    else
        Throw(BoundsError(k,i))
    end
end
TensorKit.findindex(::TensorKit.SectorValues{SU2kIrrep{k}}, s::SU2kIrrep{k}) where {k} = twice(s.j) + 1

TensorKit.dim(s::SU2kIrrep) = q_number(twice(s.j) + 1, q(s))
TensorKit.FusionStyle(::Type{<:SU2kIrrep}) = SimpleFusion()
TensorKit.BraidingStyle(::Type{<:SU2kIrrep}) = TensorKit.Anyonic()

δ(j₁, j₂, j₃, k) = (j₃ <= j₁ + j₂) && (j₁ <= j₂ + j₃) && (j₂ <= j₃ + j₁) && (j₁ + j₂ + j₃ <= k) && isinteger(j₁+j₂+j₃)

function TensorKit.Nsymbol(sa::T, sb::T, sc::T) where {T<:SU2kIrrep}
    return δ(sa.j, sb.j, sc.j, k(sa))
end

function TensorKit.Fsymbol(s1::T, s2::T, s3::T,
    s4::T, s5::T, s6::T) where {T<:SU2kIrrep}
    Nsymbol(s1,s2,s5) && Nsymbol(s5,s3,s4) && Nsymbol(s2,s3,s6) && Nsymbol(s1,s6,s4) || return 0.0
return sqrt(dim(s5) * dim(s6)) * q_racahW(s1.j, s2.j, s4.j, s3.j, s5.j, s6.j, ComplexF64(q(s1)))
end

function TensorKit.Rsymbol(a::T, b::T, c::T) where {T<:SU2kIrrep}
    Nsymbol(a, b, c) || return 0.0
    factor = q(a)^((c.j * (c.j + 1) - a.j * (a.j + 1) - b.j * (b.j + 1)) / 2)
    return isodd(convert(Int, a.j + b.j - c.j)) ? -factor : factor
end