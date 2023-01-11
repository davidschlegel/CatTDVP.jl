
<img src="https://user-images.githubusercontent.com/24460573/211853363-d4b685f6-46ae-485a-a76b-658fabd11cf5.svg" alt="drawing" width="200"/>

# CatTDVP.jl

CatTDVP.jl is a julia framework to efficiently simulate the quantum dynamics of bosonic cat qubits.
Using a time-dependent variational ansatz in the basis of cat states, interacting cat-like systems can by simulated. Symbolic bosonic opperators are specified using the QuantumCumulants.jl package.

To obtain a numerical solution, equations derived with CatTDVP.jl can be solved with DifferentialEquations.jl.

This package is in an experimental state, please submit an issue if you encounter unexpected behavior.

## Installation

To install Symbolics.jl, use the Julia package manager:

```julia
julia> using Pkg
julia> Pkg.add("CatTDVP.jl")
```


## Example

To briefly illustrate how CatTDVP.jl works, here's how you can implement a coupled two-mode system.

```julia
using QuantumCumulants, QuantumOptics, CatTDVP

hf = FockSpace(:bosonic_mode)
a₁ = Destroy(hf, Symbol("a",1), 1)
a₂ = Destroy(hf, Symbol("a",2), 2)

G₁, G₂, η₁, η₂, J = [4.0, 4.0, 1.0, 1.0, 1.0]
# Hamiltonian: Two-photon drive in each mode plus a coherent hopping
H = G₁*(a₁^2 + (a₁')^2) + G₂*(a₂^2 + (a₂')^2) + J*(a₁*a₂' + a₁'*a₂)
# Dissipators: Two-photon loss
J = [a₁^2, a₂^2]
rates = [η₁, η₂]

# truncation order of the basis, i.e. the number of creation operators applied on the basis of each mode
order = [1, 1]; dim = prod(2 .* (order .+1))
sys = TDVPSystem(H, J, order; rates=rates)  #create TDVPSystem

#Creation of the initial state from FockBasis unsing QuantumOptics.jl
Nfock = 24
b_fock = FockBasis(Nfock)
α  = [sqrt(-1im*2*G₁/η₁), sqrt(-1im*2*G₂/η₂)] #steady-state α
# Prepare a cat in each individual mode
ψa = normalize(coherentstate(b_fock, α[1]) + coherentstate(b_fock, -α[1]))
ψb = normalize(coherentstate(b_fock, α[2]) + coherentstate(b_fock, -α[2]))
ρfock = dm(tensor(ψa, ψb))  # tensor product state

ρ_barg = to_barg_basis(ρfock, α, sys.ord)  # conversion to TDVP basis

tspan = (0.0, 2.0)  # define time span
t_list = range(0.0, 2.0, 100)  # define times to save the output
# create actual ODEProblem
using DifferentialEquations
prob = TDVPProblem(sys, α, ρ_barg.data, tspan; saveat=t_list)

#solve the system
sol = solve(prob, DFBDF(autodiff=false); initializealg=BrownFullBasicInit(), abstol=1e-6, reltol=1e-6);
```

Once the solution is obtained, we can analyze the solution either in the original basis or convert back to the Fock basis. 

