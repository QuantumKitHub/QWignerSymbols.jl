"""
    RootOfUnity(n::Integer)

Custom number type to represent the ``n``-th root of unity, i.e. `RootOfUnity(n)^n ≈ 1`.
For `SU2qIrrep` with the deformation parameter `q` being a root of unity, the convention is
that ``q = exp(2πi / (k + 2))``, with ``k`` the level of ``SU(2)_k``. 
Therefore, the ``n``-th root of unity corresponds to level ``n - 2``.
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
