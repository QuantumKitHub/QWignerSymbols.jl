# TensorKitSectors extension
# Rep(su(2)_q): q-deformed SU(2) irreps
# -------------------

"""
    abstract type SU2_{k} <: Group

Abstract type for the ``q``-deformed version of ``SU(2)`` at level ``k``.
"""
abstract type SU2_{k} <: TensorKitSectors.Group end

const SU₂_₂ = SU2_{2}
const SU₂_₃ = SU2_{3}
const SU₂_₄ = SU2_{4}
const SU₂_₅ = SU2_{5}
type_repr(::Type{SU₂_₂}) = "SU₂_₂"
type_repr(::Type{SU₂_₃}) = "SU₂_₃"
type_repr(::Type{SU₂_₄}) = "SU₂_₄"
type_repr(::Type{SU₂_₅}) = "SU₂_₅"

"""
    struct SU2qIrrep{Q} <: Sector
    SU2qIrrep{Q}(j::Real)
    SU2qIrrep(j::Real, Q::Number)

Represents q-deformed irreps of the group ``SU₂``, i.e. the irreps of the ``Rep(su(2)_q)`` quantum group.
The irrep is labelled by a half integer `j` which can be entered as an abitrary `Real`,
but is stored as a `HalfInt` from the HalfIntegers.jl package.

`Q::Number` is the deformation parameter `q`, which can be either a real number or a root of unity (i.e. ``|q| = 1``).
In the latter case, there is an equivalence between these irreps and the ones of the ``SU(2)_k``, see also [`SU2kIrrep`](@ref).
The convention used here is that ``q = exp(2πi / (k + 2))``, with ``k`` the level of ``SU(2)_k``. 

## Fields
- `j::HalfInt`: the label of the irrep, which can be any non-negative half integer.
"""
struct SU2qIrrep{Q} <: Sector
    j::HalfInt
    function SU2qIrrep{Q}(j) where {Q}
        _isvalid_irrep(j, Q)
        return new{Q}(j)
    end
end

@noinline _isvalid_irrep(j::Real, ::Number) =
    (zero(j) ≤ j || throw(DomainError(j, "Not a valid SU₂ irrep")); nothing)
@noinline function _isvalid_irrep(j::Real, q::RootOfUnity)
    zero(j) ≤ j || throw(DomainError(j, "Not a valid SU₂ irrep"))
    k = q.level
    1 ≤ k || throw(DomainError(k, "Level must be positive"))
    j ≤ half(k) || throw(DomainError(j, lazy"j can be at most k/2 = $(half(k))"))
    return nothing
end

SU2qIrrep(j, q::Number) = SU2qIrrep{q}(j)

"""
    SU2kIrrep(j::Real, k::Integer)

Construct the `SU2qIrrep` where `q` is the primitive root of unity corresponding to level ``k``, 
i.e. ``q = exp(2πi / (k + 2))``. See also [`RootOfUnity`](@ref).
"""
SU2kIrrep(j, k::Integer) = SU2qIrrep(j, RootOfUnity(_root(k)))

q(a::SU2qIrrep) = q(typeof(a))
q(::Type{SU2qIrrep{Q}}) where {Q} = Q

Base.convert(T::Type{<:SU2qIrrep}, j::Real) = T(j)
Base.hash(s::SU2qIrrep, h::UInt) = hash(s.j, h)
Base.isless(s1::T, s2::T) where {T <: SU2qIrrep} = isless(s1.j, s2.j)

FusionStyle(::Type{<:SU2qIrrep}) = SimpleFusion()
BraidingStyle(::Type{<:SU2qIrrep}) = Anyonic()

TensorKitSectors.fusionscalartype(::Type{I}) where {I <: SU2qIrrep} = float(typeof(q(I)))
TensorKitSectors.sectorscalartype(::Type{I}) where {I <: SU2qIrrep} = float(typeof(q(I)))
TensorKitSectors.braidingscalartype(::Type{I}) where {I <: SU2qIrrep} = float(typeof(q(I)))

# ------------------------------------------------------------------------------------
#
unit(::Type{T}) where {T <: SU2qIrrep} = T(zero(HalfInt))
dual(s::SU2qIrrep) = s
dim(s::SU2qIrrep) = _dim(s.j, q(s))
_dim(j, q::Number) = q_number(twice(j) + 1, q)
_dim(j, q::RootOfUnity) = sinpi((2j + 1) / q.root) / sinpi(1 / q.root)

# ------------------------------------------------------------------------------------

Base.IteratorSize(::Type{SectorValues{I}}) where {I <: SU2qIrrep} =
    q(I) isa RootOfUnity ? Base.HasLength() : Base.IsInfinite()

Base.iterate(::SectorValues{I}, i::Int = 0) where {I <: SU2qIrrep} = _iterate(I, q(I), i)
_iterate(::Type{I}, q::Number, i::Int) where {I} = I(half(i)), i + 1
_iterate(::Type{I}, q::RootOfUnity, i::Int) where {I} = i >= length(values(I)) ? nothing : (I(half(i)), i + 1)

function Base.length(::SectorValues{I}) where {I <: SU2qIrrep}
    q(I) isa RootOfUnity || throw(ArgumentError("length is infinite"))
    return q(I).level + 1
end

Base.checkbounds(::Type{Bool}, vals::SectorValues{I}, i::Int) where {I <: SU2qIrrep} =
    Base.IteratorSize(typeof(vals)) === Base.IsInfinite() ? 1 ≤ i : 1 ≤ i ≤ length(vals)
Base.checkbounds(vals::SectorValues{I}, i::Int) where {I <: SU2qIrrep} =
    checkbounds(Bool, vals, i) || throw(BoundsError(vals, i))

@inline function Base.getindex(vals::SectorValues{I}, i::Int) where {I <: SU2qIrrep}
    @boundscheck Base.checkbounds(vals, i)
    return I(half(i - 1))
end
findindex(::SectorValues{I}, s::I) where {I <: SU2qIrrep} = twice(s.j) + 1

# ------------------------------------------------------------------------------------

const SU2qIrrepProdIterator{Q} = TensorKitSectors.SectorProductIterator{SU2qIrrep{Q}}
Base.IteratorSize(::Type{<:SU2qIrrepProdIterator}) = Base.HasLength()

_stop(a, b, q) = a.j + b.j
_stop(a, b, q::RootOfUnity) = min(a.j + b.j, q.level - a.j - b.j)

Base.length(it::SU2qIrrepProdIterator{Q}) where {Q} = length(abs(it.a.j - it.b.j):_stop(it.a, it.b, Q))
function Base.iterate(it::SU2qIrrepProdIterator{Q}, state = abs(it.a.j - it.b.j)) where {Q}
    return state > _stop(it.a, it.b, Q) ? nothing : (eltype(it)(state), state + 1)
end

# ------------------------------------------------------------------------------------

Nsymbol(sa::I, sb::I, sc::I) where {I <: SU2qIrrep} = q_δ(sa.j, sb.j, sc.j, q(I))

function Fsymbol(s1::I, s2::I, s3::I, s4::I, s5::I, s6::I) where {I <: SU2qIrrep}
    T = fusionscalartype(I)
    return T(sqrt(dim(s5) * dim(s6)) * q_racahW(s1.j, s2.j, s4.j, s3.j, s5.j, s6.j, q(I)))
end

function Rsymbol(a::I, b::I, c::I) where {I <: SU2qIrrep}
    T = braidingscalartype(I)
    Nsymbol(a, b, c) || return zero(T)
    factor = T(q(I)^((c.j * (c.j + 1) - a.j * (a.j + 1) - b.j * (b.j + 1)) / 2))
    return isodd(convert(Int, a.j + b.j - c.j)) ? -factor : factor
end

# ------------------------------------------------------------------------------------

Base.getindex(::TensorKitSectors.IrrepTable, G::Type{<:SU2_{k}}) where {k} =
    SU2qIrrep{RootOfUnity(_root(k))}

type_repr(::Type{SU2qIrrep{Q}}) where {Q} =
    Q isa RootOfUnity ? "Irrep[$(type_repr(SU2_{Q.level}))]" : "SU2qIrrep{$Q}"

function Base.show(io::IO, s::SU2qIrrep)
    I = typeof(s)
    if get(io, :typeinfo, nothing) !== I
        print(io, type_repr(I), "(", s.j, ")")
    else
        print(io, s.j)
    end
    return nothing
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
