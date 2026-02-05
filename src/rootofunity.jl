"""
    RootOfUnity(k::Integer)

Custom number type to represent the ``n``-th root of unity, i.e. `RootOfUnity(n)^n ≈ 1`
"""
struct RootOfUnity <: Number
    root::Int
end

Base.getproperty(q::RootOfUnity, field::Symbol) =
    field === :level ? _level(q.root) : getfield(q, field)

Base.propertynames(q::RootOfUnity) = (:root, :level)

# utility to map between "nth" root and "kth" level
_level(root::Integer) = root - 2
_root(level::Integer) = level + 2

Base.float(q::RootOfUnity) = cispi(2 // q.root)
Base.float(::Type{RootOfUnity}) = ComplexF64

Base.convert(::Type{T}, q::RootOfUnity) where {T <: Number} = T === RootOfUnity ? q : T(float(q))

Base.promote_rule(::Type{RootOfUnity}, ::Type{T}) where {T <: Number} = complex(float(T))

Base.:(*)(x₁::RootOfUnity, x₂::RootOfUnity) = float(x₁) * float(x₂)
Base.:(+)(x₁::RootOfUnity, x₂::RootOfUnity) = float(x₁) + float(x₂)
