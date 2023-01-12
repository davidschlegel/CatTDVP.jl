function liouvillian!(L::AbstractArray{M}, H::SparseArray{K, 2}, S::SparseArray{K, 2},
                    J::SparseArray{K, 3}, Jdag::SparseArray{K, 3},
                    JdagJ::SparseArray{K, 3}, γ::AbstractVector{<:Number},
                    ρ::AbstractArray{M}) where {K<:Number, M<:Number}
    #This function makes problems!!!
    comm =  -im*(H*ρ*S - S*ρ*H)
    dissipator = sum([γ[μ]*(J[:, :, μ]*ρ*Jdag[:, :, μ] -1/2*JdagJ[:,:, μ]*ρ*S -1/2*S*ρ*JdagJ[:, :, μ]) for μ ∈ eachindex(γ)])
    L .= comm + dissipator
    return nothing
end


function liouvillian!(L::AbstractArray{M},  m::TDVPMatrices, ρ::AbstractArray{M}) where M<:Number
    return liouvillian!(L, m.H, m.S, m.J, m.Jdag, m.JdagJ, m.system.rates, ρ)
end


function liouvillianderivative!(∂L::AbstractArray{M}, H::SparseArray{K, 2}, S::SparseArray{K, 2},
                    a::SparseArray{K, 3}, aH::SparseArray{K, 3},
                    aJ::SparseArray{K, 4}, Jdag::SparseArray{K, 3},
                    JdagJ::SparseArray{K, 3}, aJdagJ::SparseArray{K, 4},
                    γ::AbstractVector{<:Number},
                    ρ::AbstractArray{M}, shape::Tuple) where {M, K}

    @inbounds for k ∈ 1:Int(length(shape)/2)
        tmp = -im*(aH[:,:,k]*ρ*S - a[:,:,k]*ρ*H) 
        @inbounds for μ ∈ eachindex(γ)
            tmp +=    γ[μ]*(aJ[:, :, μ, k]*ρ*Jdag[:, :, μ] +
                        - 1/2*aJdagJ[:, :, μ, k]*ρ*S +
                        - 1/2*a[:, :, k]*ρ*JdagJ[:, :, μ])
        end
        ∂L[:,:,k] .= flip_parity_resized(tmp, k, shape)
    end
    return nothing
end




function flip_parity_resized(mat, k, shape)
    mat_desh = reshape(mat, shape)
    #nmodes = Int(length(shape)/2)
    kidx = falses(length(shape))
    kidx[k] = true #flip on the right
    arr = similar(mat_desh)
    for key in keys(mat_desh)
        newkey = flip_parity(key, Tuple(kidx), size(mat_desh)) #get the flipped index
        arr[newkey.I...] = mat_desh[key] #assign
    end
    return reshape_basis(arr, Int(length(shape)/2)) 
end


function liouvillianderivative!(∂L::AbstractArray{M},  m::TDVPMatrices, ρ::AbstractArray{M}) where M<:Number
    return liouvillianderivative!(∂L, m.H, m.S, m.a, m.aH, m.aJ, m.Jdag, m.JdagJ, m.aJdagJ, m.system.rates, ρ, shape(m.system))
end



function update!(matrices::TDVPMatrices, func::TDVPFunction, α::AbstractVector{<:Number})
    @inbounds for field ∈ [:S, :H, :κᵣ, :κₗᵣ, :J, :Jdag, :JdagJ, :a, :aH, :aJ, :aJdagJ]
        f! = getfield(func, field)
        f!(getfield(matrices, field), α)
    end
end



function make_ODE_problem(sys::TDVPSystem)
    sys_f = TDVPFunction(sys)
    mat = TDVPMatrices(sys)

    dims = 2 .* (sys.ord .+1)
    prdims = prod(dims)
    n = sys.nmodes
    ∂L = zeros(ComplexF64, prdims, prdims, n)#dense
    L = zeros(ComplexF64, prdims, prdims)#dense
    τ = zeros(ComplexF64, prdims, prdims)#dense
    C = zeros(ComplexF64, n, n)#dense
    Y = zeros(ComplexF64, n)#dense

    QGT = SparseArray{ComplexF64}(undef, prdims, prdims, n, n)#sparse

    function f(res::AbstractVector{ComplexF64}, du::AbstractVector{ComplexF64}, u::AbstractVector{ComplexF64}, p, t; verbose=false)
        #println("t:\t", t)
        α = u[1:n]#view(u, 1:n)
        dα = du[1:n]#view(du, 1:n)
        ρ = reshape(u[(n+1):(prdims^2+n)], prdims, prdims)
        #reshape(view(u, (n+1):(prdims^2+n)), prdims, prdims) #and reshape

        #update all necessary matrices: order is important here!!!
        update!(mat, sys_f, α)
        liouvillian!(L, mat, ρ)
        verbose && println("L:\t", L)
        liouvillianderivative!(∂L, mat, ρ)
        verbose && println("∂L:\t", ∂L)
        #extract them
        #S, H, κᵣ, κₗᵣ, L, Ldag, LdagL, a, aH, aL, aLdagL = [getfield(matrices, field) for field ∈ fields]
        S, κᵣ, κₗᵣ = [mat.S, mat.κᵣ, mat.κₗᵣ]
        S⁻¹ = SparseArray(pinv(Matrix(S))) #this is not the most efficient, but pinv is only defined on dense matrices

        τ = sum([κᵣ[:,:,k]*dα[k] for k ∈ 1:n])
        verbose && println("τ\t", τ)

        @tensor QGT[i,j, k, l] = κₗᵣ[i,j,k,l] - conj(κᵣ[β,i,k])*S⁻¹[β, γ]*κᵣ[γ,j,l]
        verbose && println("QQT:\t", QGT)

        #right-hand-side for α
        @inbounds for β ∈ 1:n, γ ∈ 1:n
            C[β, γ] = tr(QGT[:,:,β,γ]*ρ*S*ρ)
        end
        verbose && println("C:\t", C)
        κₗ = conj.(permutedims(κᵣ, (2,1,3)))
        verbose && println("κₗ:\t", κₗ)
        @inbounds for β ∈ 1:n
            Y[β] = tr((∂L[:,:, β] - κₗ[:,:,β]*S⁻¹*L)*ρ)
        end
        verbose && println("Y:\t", Y)
        #println(C)
        #println(Y)
        #rhs_α = abs(det(C)) < 1e-6 ? zeros(ComplexF64, n) : inv(C)*Y
        rhs_α = pinv(C)*Y
        verbose && println("α:\t", α)
        verbose && println("ρ:\t", ρ)
        verbose && println("rhs_α:\t", rhs_α)
        #right-hand-side for ρ
        rhs_ρ = S⁻¹*L*S⁻¹ - S⁻¹*τ*ρ - ρ*τ'*S⁻¹
        verbose && println("rhs_ρ:\t", rhs_ρ)
        verbose && println(rhs_ρ)

        res[1:n] = du[1:n] - rhs_α
        res[n+1:end] = du[n+1:end] - reshape(rhs_ρ, prdims^2)
        #verbose && println("residue:\t", res)
        return nothing
    end
    return f
end



function TDVPProblem(sys::TDVPSystem, u0::AbstractVector{ComplexF64}, tspan::Tuple{K, K} where K <:Real; kwargs...)
    f_dae = make_ODE_problem(sys)  # create ODEProblem
    return DAEProblem{true}(f_dae, zeros(eltype(u0), length(u0)), u0, tspan, differential_vars=trues(length(u0)), kwargs...)
end


function TDVPProblem(sys::TDVPSystem, α0::AbstractVector{ComplexF64}, ρ0::AbstractMatrix{ComplexF64}, tspan::Tuple{K, K} where K <:Real; kwargs...)\
    ishermitian(ρ0)
    return TDVPProblem(sys, [α0; ρ0...], tspan; kwargs...)
end