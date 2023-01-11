###bargman state ###

function bargmanstate(::Type{T}, b::FockBasis, alpha::Number) where T
    result = Ket(T, b)
    bargmanstate!(result, b, alpha)
    return result
end
bargmanstate(b, alpha) = bargmanstate(ComplexF64, b, alpha)


function bargmanstate!(ket::Ket, b::FockBasis, alpha::Number)
    C = eltype(ket)
    T = real(C)
    alpha = C(alpha)
    data = ket.data
    data[1] = 1 #this is the only change with respect to coherentstate!

    # Compute coefficient up to offset
    offset = b.offset
    @inbounds for n=1:offset
        data[1] *= alpha/sqrt(T(n))
    end

    # Write coefficients to state
    @inbounds for n=1:b.N-offset
        data[n+1] = data[n]*alpha/sqrt(T(n+offset))
    end

    return ket
end


###conversions###
function barg_basis(b::QuantumOptics.FockBasis, α::Number, ord::Int)
    adag = create(b)
    barg₊ = bargmanstate(b, α)
    barg₋ = bargmanstate(b, -α)
    basis_even = [adag^i*(barg₊  + (-1)^i *barg₋) for i in 0:ord]
    basis_odd = [adag^i*(barg₊ - (-1)^i *barg₋) for i in 0:ord]
    b_sub = SubspaceBasis(b, [basis_even; basis_odd])
    return b_sub
end

function barg_basis(b::QuantumOptics.CompositeBasis, α::Vector{<:Number}, ord::Vector{Int})
    nmodes = length(α)
    subbases = [barg_basis(b.bases[n], α[n], ord[n]) for n ∈ 1:nmodes]
    return tensor(subbases...)
end


function to_barg_basis(state::QuantumOptics.Operator{T, T, Matrix{F}} where {T <: FockBasis, F <: Number}, 
                        α::Number, ord::Int)
    b_Fock = state.basis_l
    b_sub = barg_basis(b_Fock, α, ord)
    Sdata = tobasisaugmented(1, [α]; ord=[ord])
    S⁻¹ = Operator(b_sub, b_sub, inv(Matrix(Sdata)))
    proj = projector(b_sub, b_Fock)
    return S⁻¹*proj*state*proj'*S⁻¹
end


function to_barg_basis(state::QuantumOptics.Operator{T, T, Matrix{F}} where {T <: CompositeBasis, F <: Number},
                        α::Vector{<:Number}, ord::Vector{Int})
    composite_basis = state.basis_l
    nmodes = length(α)
    b_sub = barg_basis(composite_basis, α, ord)
    Sdata = reshape_basis(tobasisaugmented(1, α; ord=ord), nmodes)
    S⁻¹ = Operator(b_sub, b_sub, inv(Matrix(Sdata)))
    proj = tensor([projector(sub, b) for (sub, b) ∈ zip(b_sub.bases, composite_basis.bases)]...)
    return S⁻¹*proj*state*proj'*S⁻¹
end


function to_cat_basis(state::QuantumOptics.Operator{T, T, Matrix{F}} where {T <: CompositeBasis, F <: Number},
    α::Vector{<:Number}, ord::Vector{Int})
    composite_basis = state.basis_l
    nmodes = length(α)
    b_sub = barg_basis(composite_basis, α, ord)
    proj = tensor([projector(QuantumOpticsBase.orthonormalize(sub), b) for (sub, b) ∈ zip(b_sub.bases, composite_basis.bases)]...)
    return proj*state*proj'
end


# function to_barg_basis(u::AbstractVector, b::QuantumOptics.FockBasis)
#     if isinteger(sqrt(length(u)-1))
#         dim = Int(sqrt(length(u)-1))
#         ord = Int(dim/2 - 1)
#         ρ = reshape(u[2:end], (dim, dim))
#     else
#         #@error("Density matrix is not square. Aborting.")
#         idx = SVector{length(u)-1}(2:length(u))
#         ρ = convert(Matrix{eltype(u)}, SHermitianCompact(u[idx]))
#         ord = Int(size(ρ)[1]/2 -1)
#     end
#     #@error("Density matrix is not square. Aborting.")
#     α = u[1]
#     barg_sub = barg_basis(b, α, ord)
#     return Operator(barg_sub, barg_sub, ρ)
# end

function barg_to_fock(op::QuantumOptics.Operator{T, T, Matrix{F}} where {T <: SubspaceBasis, F <: Number})
    b_sub = op.basis_l
    b_fock = b_sub.superbasis
    return projector(b_fock, b_sub)*op*projector(b_sub, b_fock)
end


function barg_to_fock(op::QuantumOptics.Operator{T, T, Matrix{F}} where {T <: CompositeBasis, F <: Number})
    bases = op.basis_l.bases
    proj = tensor([projector(b.superbasis, b) for b ∈ bases]...)
    return proj*op*proj'
end