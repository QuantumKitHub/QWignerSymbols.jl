# Q-numbers
# ---------
q_number(n::Integer, q::Number) = sum(i -> q^((n + 1) / 2 - i), 1:n)
q_number(n::Number, q::Number) = q_number(Int(n), q)

q_factorial(n::Integer, q::Number) = prod(Base.Fix2(q_number, q), 1:n; init = one(float(q)))
q_factorial(n::Number, q::Number) = q_factorial(Int(n), q)

function q_binomial(n::Integer, m::Integer, q::Number)
    return q_factorial(n, q) / (q_factorial(n, q) * q_factorial(n - m, q))
end

# Wigner symbols
# --------------
q_δ(j₁, j₂, j₃, q::Number) = δ(j₁, j₂, j₃)
q_δ(j₁, j₂, j₃, q::RootOfUnity) = δ(j₁, j₂, j₃) && j₁ + j₂ + j₃ ≤ level(q)

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
    (q_δ(α̂₁..., q) && q_δ(α̂₂..., q) && q_δ(α̂₃..., q) && q_δ(α̂₄..., q)) || return 0.0

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
    q_δ(j₁, j₂, j₃, q) || return 0.0

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
