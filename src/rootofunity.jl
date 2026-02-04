"""
    RootOfUnity(k::Integer)

Custom number type to represent the ``k``-th root of unity, i.e. `RootOfUnity(k)^k ≈ 1`
"""
struct RootOfUnity <: Number
    k::Int
end

level(q::RootOfUnity) = q.k

# TODO: level and root of unity convention are clashing here?
Base.float(q::RootOfUnity) = cispi(2 // (level(q) + 2))
Base.float(::Type{RootOfUnity}) = ComplexF64

Base.convert(::Type{T}, q::RootOfUnity) where {T <: Number} = T === RootOfUnity ? q : T(float(q))

Base.promote_rule(::Type{RootOfUnity}, ::Type{T}) where {T <: Number} = complex(float(T))
