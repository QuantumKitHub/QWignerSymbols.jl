# SU(2)_k irreps
# -------------------
const SU2kIrrep{k} = SU2qIrrep{RootOfUnity{k}}
# const SU2kIrrep{k} = SU2qIrrep{RootOfUnity(k)} # won't work

Base.convert(::Type{ComplexF64}, ::RootOfUnity{k}) where {k} = cispi(2 // (k + 2)) # this is q^2
# ComplexF64(q::RootOfUnity)= cispi(2 // (q.k + 2))

q(::Type{SU2kIrrep{k}}) where {k} = RootOfUnity{k}() # might need to switch this, but influences q_number function
# q(a::SU2kIrrep{k}) where {k} = convert(ComplexF64, q(typeof(a)))
# q(::Type{SU2kIrrep{k}}) where {k} = ComplexF64(RootOfUnity(k))
SU2kIrrep(j, k::Int) = SU2kIrrep{k}(j)
level(::SU2kIrrep{k}) where {k} = k

# TensorKitSectors.sectorscalartype(::Type{<:SU2kIrrep}) = ComplexF64

# sector values
Base.IteratorSize(::Type{SectorValues{SU2kIrrep{k}}}) where {k} = Base.HasLength()
Base.length(::SectorValues{SU2kIrrep{k}}) where {k} = k + 1
function Base.iterate(::SectorValues{SU2kIrrep{k}}, i::Int = 0) where {k}
    i >= k + 1 && return nothing
    return (SU2kIrrep{k}(half(i)), i + 1)
end
function Base.getindex(::SectorValues{SU2kIrrep{k}}, i::Int) where {k}
    if i <= k + 1
        return SU2kIrrep{k}(half(i - 1))
    else
        throw(BoundsError(k, i))
    end
end

# sector fusion
function _δ(j₁, j₂, j₃, k)
    return (j₃ <= j₁ + j₂) && (j₁ <= j₂ + j₃) && (j₂ <= j₃ + j₁) && (j₁ + j₂ + j₃ <= k) && isinteger(j₁ + j₂ + j₃)
end

function Nsymbol(sa::T, sb::T, sc::T) where {T <: SU2kIrrep}
    return _δ(sa.j, sb.j, sc.j, level(sa))
end

const SU2kIrrepProdIterator{k} = TensorKitSectors.SectorProductIterator{SU2kIrrep{k}}
Base.IteratorSize(::Type{<:SU2kIrrepProdIterator}) = Base.HasLength()
Base.length(it::SU2kIrrepProdIterator) = length(abs(it.a.j - it.b.j):min(it.a.j + it.b.j, level(it.a) - it.a.j - it.b.j))
function Base.iterate(it::SU2kIrrepProdIterator, state = abs(it.a.j - it.b.j))
    return state > min(it.a.j + it.b.j, level(it.a) - it.a.j - it.b.j) ? nothing : (eltype(it)(state), state + 1)
end

# dim(s::SU2kIrrep) = Float64(q_number(twice(s.j) + 1, q(typeof(s)))) # q^2 needed

# function Fsymbol(
#         s1::T, s2::T, s3::T, s4::T, s5::T, s6::T
#     ) where {T <: SU2kIrrep}
#     Nsymbol(s1, s2, s5) && Nsymbol(s5, s3, s4) && Nsymbol(s2, s3, s6) && Nsymbol(s1, s6, s4) || return zero(TensorKitSectors._Fscalartype(T))
#     return sqrt(dim(s5) * dim(s6)) * q_racahW(s1.j, s2.j, s4.j, s3.j, s5.j, s6.j, q(typeof(s1))^2) # q^2 needed
# end

# function Rsymbol(a::I, b::I, c::I) where {I <: SU2kIrrep}
#     Nsymbol(a, b, c) || return zero(sectorscalartype(I))
#     factor = q(I)^((c.j * (c.j + 1) - a.j * (a.j + 1) - b.j * (b.j + 1))) # q^2 required
#     return isodd(convert(Int, a.j + b.j - c.j)) ? -factor : factor
# end
