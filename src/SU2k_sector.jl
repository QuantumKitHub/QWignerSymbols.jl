"""SU2kIrrep{k}(j)

Represents SU(2)_k irreps.

Type parameters
- `k` : integer level (k >= 1).

Fields
- `j::HalfInt` : spin label (half-integer) satisfying 0 <= j <= k/2.

Notes
- The allowed irreps for a given level `k` are j = 0, 1/2, 1, 3/2, ..., k/2.
- This type implements the `TensorKitSectors.Sector` API so it can be used as a sector in TensorKit tensors.
- This is a special case of q-deformed SU(2) irreps with q^2 = exp(2πi/(k+2)) a root of unity.
"""
struct SU2kIrrep{k} <: TensorKitSectors.Sector
    j::HalfInt
    function SU2kIrrep{k}(j) where {k}
        j >= zero(j) || error("j=$j is not a positive half-integer")
        j <= half(k) || error("j=$j can be at most k/2=$(half(k))")
        return new{k}(j)
    end
end

function SU2kIrrep(j, k::Int)
    k >= 1 || error("Level k must be positve")
    return SU2kIrrep{k}(j)
end

k(::SU2kIrrep{K}) where {K} = K # TODO: give better name?
q(s::SU2kIrrep{K}) where {K} = exp(2*π*im/(k(s)+2)) # TODO: give better name? this is actually q^2

Base.hash(s::SU2kIrrep, h::UInt) = hash(s.j, h)
Base.isless(s1::T, s2::T) where {T <: SU2kIrrep} = isless(s1.j, s2.j)
Base.convert(T::Type{<:SU2kIrrep}, j::Real) = T(j)

TensorKitSectors.unit(::Type{T}) where {T <: SU2kIrrep} = T(zero(HalfInt))
TensorKitSectors.dual(s::SU2kIrrep) = s

function TensorKitSectors.:⊗(s1::T, s2::T) where {T<:SU2kIrrep}
    return TensorKitSectors.SectorSet{T}(abs(s1.j - s2.j):min(s1.j + s2.j, k(s1) - s1.j - s2.j))
end

Base.IteratorSize(::Type{<:TensorKitSectors.SectorValues{<:SU2kIrrep}}) = Base.HasLength()
Base.length(::TensorKitSectors.SectorValues{SU2kIrrep{k}}) where {k} = k + 1
function Base.iterate(::TensorKitSectors.SectorValues{SU2kIrrep{k}}, i=0) where {k}
    if i >= k + 1
        return nothing
    end
    return (SU2kIrrep{k}(half(i)), i + 1)
end

function Base.getindex(::TensorKitSectors.SectorValues{SU2kIrrep{k}}, i::Int) where {k}
    if i <= k + 1
        return SU2kIrrep(half(i - 1), k)
    else
        throw(BoundsError(k,i))
    end
end
findindex(::TensorKitSectors.SectorValues{SU2kIrrep{k}}, s::SU2kIrrep{k}) where {k} = twice(s.j) + 1

TensorKitSectors.FusionStyle(::Type{<:SU2kIrrep}) = SimpleFusion()
TensorKitSectors.sectorscalartype(::Type{<:SU2kIrrep}) = ComplexF64
TensorKitSectors.BraidingStyle(::Type{<:SU2kIrrep}) = TensorKitSectors.Anyonic()

function WignerSymbols.δ(j₁, j₂, j₃, k)
    return (j₃ <= j₁ + j₂) && (j₁ <= j₂ + j₃) && (j₂ <= j₃ + j₁) && (j₁ + j₂ + j₃ <= k) && isinteger(j₁+j₂+j₃)
end

TensorKitSectors.dim(s::SU2kIrrep) = q_number(twice(s.j) + 1, q(s))

function TensorKitSectors.Nsymbol(sa::T, sb::T, sc::T) where {T<:SU2kIrrep}
    return δ(sa.j, sb.j, sc.j, k(sa))
end

function TensorKitSectors.Fsymbol(s1::T, s2::T, s3::T,
    s4::T, s5::T, s6::T) where {T<:SU2kIrrep}
    Nsymbol(s1,s2,s5) && Nsymbol(s5,s3,s4) && Nsymbol(s2,s3,s6) && Nsymbol(s1,s6,s4) || return 0.0
    return sqrt(dim(s5) * dim(s6)) * q_racahW(s1.j, s2.j, s4.j, s3.j, s5.j, s6.j, q(s1))
end

function TensorKitSectors.Rsymbol(a::T, b::T, c::T) where {T<:SU2kIrrep}
    Nsymbol(a, b, c) || return zero(TensorKitSectors.sectorscalartype(T))
    factor = q(a)^((c.j * (c.j + 1) - a.j * (a.j + 1) - b.j * (b.j + 1)) / 2)
    return isodd(convert(Int, a.j + b.j - c.j)) ? -factor : factor
end

# does this still hold for SU2k?
function TensorKitSectors.fusiontensor(a::T, b::T, c::T) where {T<:SU2kIrrep}
    da = twice(a.j) + 1
    db = twice(b.j) + 1
    dc = twice(c.j) + 1
    C = Array{Float64}(undef, da, db, dc, 1)
    ja, jb, jc = a.j, b.j, c.j

    for kc in 1:dc, kb in 1:db, ka in 1:da
        C[ka, kb, kc, 1] = q_clebschgordan(ja, ja + 1 - ka, jb, jb + 1 - kb, jc,
                                           jc + 1 - kc, q(T))
    end

    return C
end
