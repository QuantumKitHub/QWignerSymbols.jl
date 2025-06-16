module QWignerSymbols

using WignerSymbols: δ, reorder6j
using HalfIntegers
using TensorKit

export Phase
export q_number, q_factorial, q_binomial
export q_wigner3j, q_clebschgordan, q_wigner6j, q_racahW
export SU2qIrrep

struct Phase
    p::Number
    function Phase(p::Number)
        return new(p % (2*π))
    end
end

Base.:*(p::Phase, q::Phase) = Phase(p.p + q.p)
Base.:+(p::Phase, q::Phase) = cis(p.p) + cis(q.p)
Base.:+(p::Phase, x::Number) = cis(p.p) + x
Base.:+(x::Number, p::Phase) = cis(p.p) + x
Base.:-(p::Phase, x::Number) = cis(p.p) - x
Base.:-(x::Number, p::Phase) = cis(p.p) - x
Base.:^(p::Phase, n::Number) = Phase(n * p.p)
Base.isone(p::Phase) = iszero(p.p) || iszero(p.p % (2*π))
Base.real(p::Phase) = cos(p.p)
Base.ComplexF64(p::Phase) = cis(p.p)


# Q-numbers
# ---------
q_number(n::Integer, q::Number) = Float64(isone(q) ? n : sum(i -> q^((n + 1) / 2 - i), 1:n))
q_number(n::Integer, q::Phase) = real(sum(i -> q^((n + 1) / 2 - i), 1:n))
q_number(n::Number, q::Number) = q_number(Int(n), q)
q_number(n::Number, q::Phase) = q_number(Int(n), q)

q_factorial(n::Integer, q::Number) = prod(n -> q_number(n, q), 1:n; init=1.0)
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
             *(sqrt.((q_factorial.((j₁ - m₁, j₁ + m₁, j₂ - m₂, j₂ + m₂, j₃ - m₃, j₃ + m₃),
                                   q)))...)
    iszero(factor) && return factor
    term = zero(factor)
    for n in
        ceil(max(0, -(j₃ - j₂ + m₁), -(j₃ - j₁ - m₂))):floor(min(j₁ + j₂ - j₃, j₁ - m₁,
                                                                 j₂ + m₂))
        term += (-1)^n * q^(-n * (j₁ + j₂ + j₃ + 1) / 2) /
                *(q_factorial.((n, j₁ - m₁ - n, j₂ + m₂ - n, j₁ + j₂ - j₃ - n,
                                j₃ - j₂ + m₁ + n, j₃ - j₁ - m₂ + n), q)...)
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

function q_wigner6j(j₁, j₂, j₃,
                    j₄, j₅, j₆, q)
    α̂₁ = (j₁, j₂, j₃)
    α̂₂ = (j₁, j₆, j₅)
    α̂₃ = (j₂, j₄, j₆)
    α̂₄ = (j₃, j₄, j₅)

    # check triangle conditions
    if !(δ(α̂₁...) && δ(α̂₂...) && δ(α̂₃...) && δ(α̂₄...))
        return 0.0
    end
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
    if !δ(j₁, j₂, j₃)
        return 0.0
    end

    return *(sqrt.(q_factorial.((j₁ + j₂ - j₃, j₁ - j₂ + j₃, -j₁ + j₂ + j₃), q))...) /
           sqrt(q_factorial(j₁ + j₂ + j₃ + 1, q))
end

function compute6jseries(β₁, β₂, β₃, α₁, α₂, α₃, α₄, q)
    krange = max(α₁, α₂, α₃, α₄):min(β₁, β₂, β₃)

    s = 0.0
    for n in max(α₁, α₂, α₃, α₄):min(β₁, β₂, β₃)
        num = iseven(n) ? q_factorial(n + 1, q) : -q_factorial(n + 1, q)
        den = *(q_factorial.((n - α₁, n - α₂, n - α₃, n - α₄, β₁ - n, β₂ - n, β₃ - n),
                             q)...)
        s += num / den
    end
    return s
end

# TensorKit extension
# -------------------
struct SU2qIrrep{Q} <: TensorKit.Sector
    j::HalfInt
    function SU2qIrrep{Q}(j) where {Q}
        j >= zero(j) || error("Not a valid SU₂ irrep")
        return new{Q}(j)
    end
end

SU2qIrrep(j, q::Number) = SU2qIrrep{q}(j)
q(::Type{SU2qIrrep{Q}}) where {Q} = Q
Base.convert(T::Type{<:SU2qIrrep}, j::Real) = T(j)

Base.one(::Type{T}) where {T<:SU2qIrrep} = T(zero(HalfInt))
Base.conj(s::SU2qIrrep) = s
function TensorKit.:⊗(s1::T, s2::T) where {T<:SU2qIrrep}
    return TensorKit.SectorSet{T}(abs(s1.j - s2.j):(s1.j + s2.j))
end
Base.IteratorSize(::Type{<:TensorKit.SectorValues{<:SU2qIrrep}}) = Base.IsInfinite()
function Base.iterate(::TensorKit.SectorValues{SU2qIrrep{Q}}, i=0) where {Q}
    return (SU2qIrrep{Q}(half(i)), i + 1)
end

TensorKit.FusionStyle(::Type{<:SU2qIrrep}) = SimpleFusion()
TensorKit.BraidingStyle(::Type{<:SU2qIrrep}) = TensorKit.Anyonic()
Base.isreal(::Type{<:SU2qIrrep{Q}}) where {Q} = isreal(Q)
function TensorKit.Nsymbol(sa::T, sb::T, sc::T) where {T<:SU2qIrrep}
    return δ(sa.j, sb.j, sc.j)
end

function TensorKit.Fsymbol(s1::T, s2::T, s3::T,
                           s4::T, s5::T, s6::T) where {T<:SU2qIrrep}
    return sqrt(dim(s5) * dim(s6)) * q_racahW(s1.j, s2.j, s4.j, s3.j, s5.j, s6.j, q(T))
end

function TensorKit.Rsymbol(a::T, b::T, c::T) where {T<:SU2qIrrep}
    Nsymbol(a, b, c) || return 0.0
    factor = q(T)^((c.j * (c.j + 1) - a.j * (a.j + 1) - b.j * (b.j + 1)) / 2)
    return isodd(convert(Int, a.j + b.j - c.j)) ? -factor : factor
end

TensorKit.dim(s::SU2qIrrep) = q_number(twice(s.j) + 1, q(typeof(s)))

function TensorKit.fusiontensor(a::T, b::T, c::T) where {T<:SU2qIrrep}
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

Base.hash(s::SU2qIrrep, h::UInt) = hash(s.j, h)
Base.isless(s1::T, s2::T) where {T<:SU2qIrrep} = isless(s1.j, s2.j)

# additional specialisations because dim does not return Int
function Base.axes(V::GradedSpace{I}, c::I) where {I<:SU2qIrrep}
    offset = 0
    for c′ in sectors(V)
        c′ == c && break
        offset += (twice(c′.j) + 1) * dim(V, c′)
    end
    return (offset + 1):(offset + (twice(c.j) + 1) * dim(V, c))
end

end
