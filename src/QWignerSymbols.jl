module QWignerSymbols

# imports
# --------
using WignerSymbols
using WignerSymbols: δ, reorder6j
using HalfIntegers
using TensorKitSectors
import TensorKitSectors:
    Nsymbol, Fsymbol, Rsymbol, fusiontensor,
    unit, dual, dim,
    FusionStyle, BraidingStyle,
    findindex

# exports
# -------
export q_number, q_factorial, q_binomial
export q_wigner3j, q_clebschgordan, q_wigner6j, q_racahW
export SU2qIrrep, SU2kIrrep

# includes
# --------
include("SU2k_sector.jl")

# Q-numbers
# ---------
function q_number(n::Integer, q::Number)
    if isa(q, Real)
        return Float64(isone(q) ? n : sum(i -> q^((n + 1) / 2 - i), 1:n))
    elseif isa(q, ComplexF64)
        norm(q) ≈ 1.0 || throw(DomainError(norm(q), "q must be either real or a U₁ phase"))
        return real(isone(q) ? n : sum(i -> q^((n + 1) / 2 - i), 1:n))
    end
end
q_number(n::Number, q::Number) = q_number(Int(n), q)

q_factorial(n::Integer, q::Number) = prod(n -> q_number(n, q), 1:n; init = 1.0)
q_factorial(n::Number, q::Number) = q_factorial(Int(n), q)

function q_binomial(n::Integer, m::Integer, q::Number)
    return q_factorial(n, q) / (q_factorial(n, q) * q_factorial(n - m, q))
end

# Wigner symbols
# --------------
function q_wigner3j(j₁, j₂, j₃, m₁, m₂, m₃, q::Number)
    if !δ(j₁, j₂, j₃) || !iszero(m₁ + m₂ + m₃)
        return 0.0
    end
    factor = q^(((j₁ + j₂ - j₃) * (j₁ + j₂ + j₃ + 1) + 2 * (j₁ * m₂ - j₂ * m₁)) / 4) *
        Δ(j₁, j₂, j₃, q) *
        *(sqrt.((q_factorial.((j₁ - m₁, j₁ + m₁, j₂ - m₂, j₂ + m₂, j₃ - m₃, j₃ + m₃), q)))...)
    iszero(factor) && return factor

    term = zero(factor)
    nrange = ceil(max(0, -(j₃ - j₂ + m₁), -(j₃ - j₁ - m₂))):floor(min(j₁ + j₂ - j₃, j₁ - m₁, j₂ + m₂))
    for n in nrange
        term += (-1)^n * q^(-n * (j₁ + j₂ + j₃ + 1) / 2) /
            *(q_factorial.((n, j₁ - m₁ - n, j₂ + m₂ - n, j₁ + j₂ - j₃ - n, j₃ - j₂ + m₁ + n, j₃ - j₁ - m₂ + n), q)...)
    end
    result = factor * term
    return isodd(Int(j₁ - j₂ - m₃)) ? -result : result
end

function q_clebschgordan(j₁, m₁, j₂, m₂, j₃, m₃, q::Number)
    s = q_wigner3j(j₁, j₂, j₃, m₁, m₂, -m₃, q)
    iszero(s) && return s
    s *= sqrt(q_number(2j₃ + one(j₃), q))
    return isodd(Int(j₁ - j₂ + m₃)) ? -s : s
end

function q_wigner6j(
        j₁, j₂, j₃,
        j₄, j₅, j₆, q
    )
    α̂₁ = (j₁, j₂, j₃)
    α̂₂ = (j₁, j₆, j₅)
    α̂₃ = (j₂, j₄, j₆)
    α̂₄ = (j₃, j₄, j₅)

    # check triangle conditions
    (δ(α̂₁...) && δ(α̂₂...) && δ(α̂₃...) && δ(α̂₄...)) || return 0.0

    # reduce
    α₁ = convert(UInt, +(α̂₁...))
    α₂ = convert(UInt, +(α̂₂...))
    α₃ = convert(UInt, +(α̂₃...))
    α₄ = convert(UInt, +(α̂₄...))
    β₁ = convert(UInt, j₁ + j₂ + j₄ + j₅)
    β₂ = convert(UInt, j₁ + j₃ + j₄ + j₆)
    β₃ = convert(UInt, j₂ + j₃ + j₅ + j₆)

    (β₁, β₂, β₃, α₁, α₂, α₃, α₄) = reorder6j(β₁, β₂, β₃, α₁, α₂, α₃, α₄)

    s = Δ(α̂₁..., q) * Δ(α̂₂..., q) * Δ(α̂₃..., q) * Δ(α̂₄..., q)

    return s *= compute6jseries(β₁, β₂, β₃, α₁, α₂, α₃, α₄, q)
end

function q_racahW(j₁, j₂, J, j₃, J₁₂, J₂₃, q::Number)
    s = q_wigner6j(j₁, j₂, J₁₂, j₃, J, J₂₃, q)
    if !iszero(s) && isodd(convert(Int, j₁ + j₂ + j₃ + J))
        return -s
    else
        return s
    end
end

function Δ(j₁, j₂, j₃, q)
    δ(j₁, j₂, j₃) || return 0.0

    return *(sqrt.(q_factorial.((j₁ + j₂ - j₃, j₁ - j₂ + j₃, -j₁ + j₂ + j₃), q))...) /
        sqrt(q_factorial(j₁ + j₂ + j₃ + 1, q))
end

function compute6jseries(β₁, β₂, β₃, α₁, α₂, α₃, α₄, q)
    s = 0.0
    krange = max(α₁, α₂, α₃, α₄):min(β₁, β₂, β₃)
    for n in krange
        num = iseven(n) ? q_factorial(n + 1, q) : -q_factorial(n + 1, q)
        den = *(q_factorial.((n - α₁, n - α₂, n - α₃, n - α₄, β₁ - n, β₂ - n, β₃ - n), q)...)
        s += num / den
    end
    return s
end

# TensorKitSectors extension
# Rep(su(2)_q): q-deformed SU(2) irreps
# -------------------
struct SU2qIrrep{Q} <: Sector
    j::HalfInt
    function SU2qIrrep{Q}(j) where {Q}
        j >= zero(j) || throw(DomainError(j, "Not a valid SU₂ irrep"))
        return new{Q}(j)
    end
end

SU2qIrrep(j, q::Number) = SU2qIrrep{q}(j)
q(::Type{SU2qIrrep{Q}}) where {Q} = Q

Base.convert(T::Type{<:SU2qIrrep}, j::Real) = T(j)
Base.hash(s::SU2qIrrep, h::UInt) = hash(s.j, h)
Base.isless(s1::T, s2::T) where {T <: SU2qIrrep} = isless(s1.j, s2.j)

TensorKitSectors.sectorscalartype(::Type{<:SU2qIrrep}) = Float64

# sector values
Base.IteratorSize(::Type{SectorValues{I}}) where {I <: SU2qIrrep} = Base.IsInfinite()
Base.iterate(::SectorValues{I}, i::Int = 0) where {I <: SU2qIrrep} = (I(half(i)), i + 1)
function Base.getindex(::SectorValues{I}, i::Int) where {I <: SU2qIrrep}
    return 1 <= i ? I(half(i - 1)) : throw(BoundsError(values(I), i))
end
findindex(::SectorValues{I}, s::I) where {I <: SU2qIrrep} = twice(s.j) + 1

# sector fusion
FusionStyle(::Type{<:SU2qIrrep}) = SimpleFusion()

Nsymbol(sa::I, sb::I, sc::I) where {I <: SU2qIrrep} = δ(sa.j, sb.j, sc.j)

const SU2qIrrepProdIterator{Q} = TensorKitSectors.SectorProductIterator{SU2qIrrep{Q}}
Base.IteratorSize(::Type{<:SU2qIrrepProdIterator}) = Base.HasLength()
Base.length(it::SU2qIrrepProdIterator) = length(abs(it.a.j - it.b.j):(it.a.j + it.b.j))
function Base.iterate(it::SU2qIrrepProdIterator, state = abs(it.a.j - it.b.j))
    return state > (it.a.j + it.b.j) ? nothing : (eltype(it)(state), state + 1)
end

unit(::Type{T}) where {T <: SU2qIrrep} = T(zero(HalfInt))
dual(s::SU2qIrrep) = s
dim(s::SU2qIrrep) = q_number(twice(s.j) + 1, q(typeof(s)))

Fsymbol(s1::I, s2::I, s3::I, s4::I, s5::I, s6::I) where {I <: SU2qIrrep} =
    sqrt(dim(s5) * dim(s6)) * q_racahW(s1.j, s2.j, s4.j, s3.j, s5.j, s6.j, q(I))

# sector braiding
BraidingStyle(::Type{<:SU2qIrrep}) = Anyonic()

function Rsymbol(a::I, b::I, c::I) where {I <: SU2qIrrep}
    Nsymbol(a, b, c) || return 0.0
    factor = q(I)^((c.j * (c.j + 1) - a.j * (a.j + 1) - b.j * (b.j + 1)) / 2)
    return isodd(convert(Int, a.j + b.j - c.j)) ? -factor : factor
end

# TODO: this seems to not be compatible with testsuite
# function fusiontensor(a::I, b::I, c::I) where {I <: SU2qIrrep}
#     da = twice(a.j) + 1
#     db = twice(b.j) + 1
#     dc = twice(c.j) + 1
#     C = Array{sectorscalartype(I)}(undef, da, db, dc, 1)
#     ja, jb, jc = a.j, b.j, c.j
#
#     for kc in 1:dc, kb in 1:db, ka in 1:da
#         C[ka, kb, kc, 1] = q_clebschgordan(
#             ja, ja + 1 - ka, jb, jb + 1 - kb, jc, jc + 1 - kc, q(I)
#         )
#     end
#
#     return C
# end

# additional specialisations because dim does not return Int
# function Base.axes(V::GradedSpace{I}, c::I) where {I <: SU2qIrrep}
#     offset = 0
#     for c′ in sectors(V)
#         c′ == c && break
#         offset += (twice(c′.j) + 1) * dim(V, c′)
#     end
#     return (offset + 1):(offset + (twice(c.j) + 1) * dim(V, c))
# end

end
