struct TDVPSystem
    hamiltonian::Union{QuantumCumulants.QNumber, QuantumCumulants.SNuN}
    lossoperators::Vector{T} where T <:Union{QuantumCumulants.QNumber, QuantumCumulants.SNuN}
    rates::Vector{K} where K <: Number
    ord::Union{Vector{Int}, Int}
    nmodes::Int
end

function 
    TDVPSystem(hamiltonian, lossoperators, ord::Int; rates::Vector{<:Number}=ones(length(lossoperators)))
    return TDVPSystem(hamiltonian, lossoperators, rates, ord, 1)
end

function TDVPSystem(hamiltonian, lossoperators, ord::Vector{Int}; rates::Vector{<:Number}=ones(length(lossoperators)))
     return TDVPSystem(hamiltonian, lossoperators, rates, ord, length(ord))
end

function Base.show(io::IO, sys::TDVPSystem)
    printstyled("\tTDVPSystem:"; color = :green, bold=true)
    printstyled("\n Hamiltonian:\t\t"; color = :blue, bold=true)
    Base.show(io, sys.hamiltonian)
    printstyled("\nLoss operators:\t\t"; color = :blue, bold=true)
    Base.show(io, sys.lossoperators)
    printstyled("\nLoss operator rates:\t"; color = :blue, bold=true)
    Base.show(io, sys.rates)
    printstyled("\nOrder of the basis:\t"; color = :blue, bold=true)
    Base.show(io, sys.ord)
    printstyled("\nNumber of modes:\t"; color = :blue, bold=true)
    Base.show(io, sys.nmodes)
end


"""
Stores all Matrix primitives that only depend on α.
"""
mutable struct TDVPMatrices
    const system::TDVPSystem

#    α::Array{ComplexF64}
#    ρ::Array{ComplexF64, 2}

    S::SpArray{ComplexF64, 2}
    H::SpArray{ComplexF64, 2}
    κᵣ::SpArray{ComplexF64, 3}
    κₗᵣ::SpArray{ComplexF64, 4}
    J::SpArray{ComplexF64, 3}
#    ∂Liou::SpArray{ComplexF64, 3}
#    Liou::SpArray{ComplexF64, 3}
#    QGT::SpArray{ComplexF64, 4}
    Jdag::SpArray{ComplexF64, 3}
    JdagJ::SpArray{ComplexF64, 3}
    a::SpArray{ComplexF64, 3}
    aH::SpArray{ComplexF64, 3}
    aJ::SpArray{ComplexF64, 4}
    aJdagJ::SpArray{ComplexF64, 4}
end


function Base.show(io::IO, state::TDVPMatrices)
    Base.show(io, state.system)
    printstyled("\n\tTDVPMatrices:\n"; color = :green, bold=true)
    for field in [:S, :H, :κᵣ, :κₗᵣ, :J, :Jdag, :JdagJ, :a, :aH, :aJ, :aJdagJ]
        printstyled(string(field)*"\t"; color = :blue, bold=true)
    end
    #Base.show(io, state.α)
    #printstyled("\n ρ:\t\t"; color = :blue, bold=true)
    #Base.show(io, state.ρ)
end


struct TDVPFunction
    system::TDVPSystem

    S::RGF
    H::RGF
    κᵣ::RGF
    κₗᵣ::RGF
    J::RGF
    Jdag::RGF
    JdagJ::RGF
    a::RGF
    aH::RGF
    aJ::RGF
    aJdagJ::RGF
end

function Base.show(io::IO, state::TDVPFunction)
    printstyled("\tTDVPFunction of runtime-generated functions:\n"; color = :green, bold=true)
    Base.show(io, state.system)
end

function shape(sys::TDVPSystem)
    r = 2 .* (sys.ord .+ 1)
    return (r..., r...)
end