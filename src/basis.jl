Base.zero(t::Expr) = 0 #this should be already covered by SymbolicUtils
Base.zero(t::Type{Expr}) = 0 #this should be already covered by SymbolicUtils
Base.zero(t::Symbol) = 0
Base.one(x::Type{Any}) = 1
Base.zero(::Type{Any}) = 0
Base.zero(t::Type{T} where T <: Union{Symbolics.Term, SymbolicUtils.Mul, SymbolicUtils.Add}) = 0
Base.zero(::Type{DataType}) = 0
Base.iszero(x::Union{Symbolics.Term, SymbolicUtils.Mul, SymbolicUtils.Add, SymbolicUtils.Sym}) = x === 0


function overlapmatrix(α::T where T)
    return [4*cosh(abs2(α)) 0; 0 4*sinh(abs2(α))]
end


function overlapmatrix(α::Vector{T}) where T
    l = length(α)
    dims = [2 for i in 1:(2*l)]
    arr = SparseArray{T}(undef, dims...)
    return overlapmatrix!(arr, α)
end

function overlapmatrix(α::Vector{T}) where T <: SymbolicUtils.Sym
    l = length(α)
    dims = [2 for i in 1:(2*l)]
    arr = SparseArray{Any}(undef, dims...)
    return overlapmatrix!(arr, α)
end

function overlapmatrix!(arr::SparseArray{K}, α::Vector{T}) where {T,K}
    l = length(α)
    t = [[1,2] for i in 1:l]
    for idx_l ∈ Base.product(t...)
        #only non-zero contribution if parity on the left = parity on the right
        arr[idx_l..., idx_l...] = overlapmatrixelement(idx_l, α)
    end
    return arr
end

function overlapmatrixelement(idx::Tuple, α::Vector{T} where T)
    return prod([id == 1 ? 4*cosh(abs2(β)) : 4*sinh(abs2(β))  for (id, β) in zip(idx, α)])
end



function tocatbasis(expr::Union{QuantumCumulants.QNumber, QuantumCumulants.SNuN}, α, S)
    #S = catoverlapmatrix(α)
    f = SymbolicUtils.operation(expr)
    if f===(+) || f===(-) # linearity
        args = map(x -> tocatbasis(x, α, S), SymbolicUtils.arguments(expr))
        return f(args...)
    elseif f === (*)
        c = expr.arg_c #prefactor
        nc = expr.args_nc #operators
        αpoly, perm = to_alpha(nc, α)
        nzkeys = nonzero_keys(S)
        newkeys = [index_perm(x,perm) for x ∈ nzkeys]
        Snew = similar(S)
        for (newk, oldk) ∈ zip(newkeys, nzkeys)
            Snew[newk] = S[oldk]
        end
        return c*αpoly*Snew
    else
        error("Unknown function $f")
    end
end


tocatbasis(expr::T where T <: QuantumCumulants.SNuN, α, S) = expr*S


function tocatbasis(op::Union{QuantumCumulants.Destroy, QuantumCumulants.Create}, α, S)
    αpoly, perm = to_alpha([op], α)
    nzkeys = nonzero_keys(S)
    newkeys = [index_perm(x,perm) for x ∈ nzkeys]
    Snew = similar(S)
    for (newk, oldk) ∈ zip(newkeys, nzkeys)
        Snew[newk] = S[oldk]
    end
    return αpoly*Snew
end

function tocatbasis(op::Union{QuantumCumulants.Create, QuantumCumulants.Destroy, QuantumCumulants.SNuN, QuantumCumulants.QNumber}, α)
    return tocatbasis(op, α, catoverlapmatrix(α))
end


function to_alpha(expr::Vector{K} where K, α::Vector{T} where T)
    out = 1;#Array{SymbolicUtils.Term}(undef, length(α))
    idxsa = zeros(Int, length(α))
    idxsadag = zeros(Int, length(α))
    for op ∈ expr
        op.aon > length(α) ? @error("Expression "*string(expr)*" contains operator that is incompatible with Dim "*string(length(α))*" of α") : nothing
        if isa(op, QuantumCumulants.Destroy)
        out *= α[op.aon]
        idxsa[op.aon] += 1
        elseif isa(op, QuantumCumulants.Create)
            out *= conj(α[op.aon])
            idxsadag[op.aon] += 1
        else
            @error("Operator "*string(op)*" is not a creation or annihilation operator!")
        end
    end
    return out, CartesianIndex([idxsadag; idxsa]...)
end


function index_perm(index1::CartesianIndex, index2::CartesianIndex; mod = 2)
    # Index addition modulo mod
    return CartesianIndex((index1.I .+ index2.I .- 1) .% mod .+1)
end


function flip_parity(index::CartesianIndex, k::NTuple{K, Bool}, shape::NTuple{K, Int}) where K
    ord = Int.(shape ./ 2) .-1 #get order from shape
    flip_idx = k .* (ord .+1) # multiply k-Tuple with order
    par_idx = CartesianIndex(flip_idx...)
    return index_perm(index, par_idx; mod = shape)
end



###up to here everything works###

function tobasisaugmented!(arr::Union{SubArray{K}, SparseArray{K}}, expr::Union{QuantumCumulants.QNumber, QuantumCumulants.SNuN}, 
                            α::Vector{T} where T, S::AbstractArray{K};
                            ord::Vector{Int}=zeros(length(α)), verbose=false) where K
    nmodes = length(α)
    a_ops = [QuantumCumulants.Destroy(FockSpace(:bosonic_mode), Symbol("a",aon), aon) for aon in 1:nmodes];
    rangesparity = [1:2 for i in 1:nmodes]
    rangesord = [0:i for i in ord]
    # four for-loops for left-parity, left-order, right-parity, right-order indices, respectively
    @inbounds for idxpar_l in Base.product(rangesparity...)
        @inbounds for idxord_l in Base.product(rangesord...)
            @inbounds for idxpar_r in Base.product(rangesparity...)
                @inbounds for idxord_r in Base.product(rangesord...)
                    idxl = idxord_l .+ (idxpar_l .-1) .* (ord .+ 1) .+ 1
                    idxr = idxord_r .+ (idxpar_r .-1) .* (ord .+ 1) .+ 1
                    arr_idx = CartesianIndex(idxl..., idxr...)
                    verbose && println("array index: ", arr_idx)    # print
                    newexpr = simplify(prod(a_ops .^ idxord_l)*expr*prod(adjoint.(a_ops) .^ idxord_r))
                    verbose && println("expression: ", newexpr)     # print
                    ord_index =  CartesianIndex(idxord_l..., idxord_r...)
                    par_index = CartesianIndex(idxpar_l..., idxpar_r...)
                    verbose && println("order index: ", ord_index)  # print
                    verbose && println("parity index: ", par_index) # print
                    catbasis_idx = index_perm(par_index, ord_index)
                    verbose && println("cat index: ", catbasis_idx) # print
                    arr[arr_idx] = tocatbasis(newexpr, α, S)[catbasis_idx]
                    verbose && println("val: ", arr[arr_idx])       # print
                end
            end
        end
    end
end


function tobasisaugmented(expr::Union{QuantumCumulants.QNumber, QuantumCumulants.SNuN}, 
                            α::Vector{T} where T, S::AbstractArray{K};
                            ord::Vector{Int}=zeros(length(α))) where K
    dims = 2 .* (ord .+ 1)
    arr = SparseArray{K}(undef, dims...,  dims... )
    tobasisaugmented!(arr, expr, α, S; ord)
    return arr
end



function tobasisaugmented(exprarr::AbstractArray{M} where M <: Union{QuantumCumulants.QNumber, QuantumCumulants.SNuN}, 
                        α::Vector{T} where T, S::AbstractArray{K};
                        ord::Vector{Int}=zeros(length(α))) where K
    dims = 2 .* (ord .+ 1)
    dimsexpr = size(exprarr)
    arr = SparseArray{K}(undef, dims...,  dims..., dimsexpr...)
    @inbounds for key in keys(exprarr)
        subarr = @view arr[CartesianIndices((dims..., dims...)), key] # @view is very important here
        tobasisaugmented!(subarr, exprarr[key], α, S; ord)
    end
    return arr
end

#tobasisaugmented(expr::Union{QuantumCumulants.QNumber, QuantumCumulants.SNuN}, α::T where T, S::AbstractArray{K} where K;
#                    ord::Int=0) = tobasisaugmented(expr, [α], S; ord = [ord])

#tobasisaugmented(expr::Union{QuantumCumulants.QNumber, QuantumCumulants.SNuN}, α::T where T;
#                    ord::Int=0) = tobasisaugmented(expr, [α], overlapmatrix([α]); ord=[ord])

tobasisaugmented(expr::Union{K, AbstractArray{K}} where K <: Union{QuantumCumulants.QNumber, QuantumCumulants.SNuN}, 
                    α::Vector{T} where T;
                    ord::Vector{Int}=zeros(length(α)))  = tobasisaugmented(expr, α, overlapmatrix(α); ord)




function right_derivative_matrix!(arr::SparseArray{K}, α::Vector{T}) where {K, T}
    arr_shape = size(arr)
    #println(arr_shape)
    nmodes = length(α)
    ord = collect(Int.(arr_shape[1:(end-1-nmodes)] ./ 2) .-1) # get order from arr
    adag_ops = [QuantumCumulants.Create(FockSpace(:bosonic_mode), Symbol("a",aon), aon) for aon in 1:nmodes];
    S = overlapmatrix(α)
    for k ∈ 1:nmodes
        #somehow, in-place does not work here with subarrays...
        newarr = tobasisaugmented(adag_ops[k], α, S; ord = ord)
        #loop to flip parity:
        for key ∈ nonzero_keys(newarr)
            #     #generate new key that flips the parity of mode k:
                 kidx = falses(nmodes*2)
                 kidx[nmodes + k] = true #flip on the right
                 newkey = flip_parity(key, Tuple(kidx), size(newarr)) #get the flipped index
                 arr[newkey.I..., k] = newarr[key] #assign
        end
    end
end


function right_derivative_matrix(α::Vector{T};  ord::Vector{Int}=zeros(length(α))) where T <: SymbolicUtils.Sym
    shape = 2 .* (ord .+ 1)
    arr = SparseArray{Any}(undef, shape... , shape..., length(α))
    right_derivative_matrix!(arr, α)
    return arr
end


# This is actually not needed in any computation
function left_derivative_matrix!(arr::SparseArray{K}, α::Vector{T}) where {K, T}
    arr_shape = size(arr)
    #println(arr_shape)
    nmodes = length(α)
    ord = collect(Int.(arr_shape[1:(end-1-nmodes)] ./ 2) .-1) # get order from arr
    a_ops = [QuantumCumulants.Destroy(FockSpace(:bosonic_mode), Symbol("a",aon), aon) for aon in 1:nmodes];
    S = overlapmatrix(α)
    for k ∈ 1:nmodes
        #somehow, in-place does not work here with subarrays...
        newarr = tobasisaugmented(a_ops[k], α, S; ord = ord)
        #loop to flip parity:
        for key ∈ nonzero_keys(newarr)
            #     #generate new key that flips the parity of mode k:
                 kidx = falses(nmodes*2)
                 kidx[k] = true #flip on the left
                 newkey = flip_parity(key, Tuple(kidx), size(newarr)) #get the flipped index
                 arr[newkey.I..., k] = newarr[key] #assign
        end
    end
end

# This is actually not needed in any computation
function left_derivative_matrix(α::Vector{T};  ord::Vector{Int}=zeros(length(α))) where T <: SymbolicUtils.Sym
    shape = 2 .* (ord .+ 1)
    arr = SparseArray{Any}(undef, shape... , shape..., length(α))
    left_derivative_matrix!(arr, α)
    return arr
end


function leftright_derivative_matrix!(arr::SparseArray{K}, α::Vector{T}) where {K, T}
    arr_shape = size(arr)
    #println(arr_shape)
    nmodes = length(α)
    ord = collect(Int.(arr_shape[1:(end-2-nmodes)] ./ 2) .-1) # get order from arr
    a_ops = [QuantumCumulants.Destroy(FockSpace(:bosonic_mode), Symbol("a",aon), aon) for aon in 1:nmodes];
    S = overlapmatrix(α)
    for kl ∈ 1:nmodes
        for kr ∈ 1:nmodes
            #somehow, in-place does not work here with subarrays...
            newarr = tobasisaugmented(a_ops[kl]*a_ops[kr]', α, S; ord = ord)
            #loop to flip parity:
            for key ∈ nonzero_keys(newarr)
                #     #generate new key that flips the parity of mode k:
                    kidx = falses(nmodes*2)
                    kidx[kl] = true #flip on the left
                    kidx[nmodes + kr] = true #flip on the right
                    newkey = flip_parity(key, Tuple(kidx), size(newarr)) #get the flipped index
                    arr[newkey.I..., kl, kr] = newarr[key] #assign
            end
        end
    end
end


function leftright_derivative_matrix(α::Vector{T};  ord::Vector{Int}=zeros(length(α))) where T <: SymbolicUtils.Sym
    shape = 2 .* (ord .+ 1)
    arr = SparseArray{Any}(undef, shape... , shape..., length(α), length(α))
    leftright_derivative_matrix!(arr, α)
    return arr
end



function reshape_basis(arr::AbstractArray{K} where K, nmodes)
    shape_arr = size(arr)
    out_dim = prod(shape_arr[1:nmodes]) #multiply all dimensions of the basis together
    remain_dim = shape_arr[2*nmodes+1:end] #remaining dimension of the array
    return reshape(arr, out_dim, out_dim, remain_dim...)
end







### TDVP type declarations and constructors ###



function TDVPMatrices(sys::TDVPSystem)
    type = ComplexF64
    #α = zeros(type, sys.nmodes)
    dims = 2 .* (sys.ord .+1)
    prdims = prod(dims)
    #ρ = zeros(type, prdims, prdims)

    S, H = [SparseArray{type}(undef, prdims, prdims) for i ∈ 1:2]
    L, Ldag, LdagL = [SparseArray{type}(undef, prdims, prdims, length(sys.lossoperators)) for i ∈ 1:3]
    κᵣ, a, aH = [SparseArray{type}(undef, prdims, prdims, sys.nmodes) for i ∈ 1:3]
    aL, aLdagL = [SparseArray{type}(undef, prdims, prdims, length(sys.lossoperators), sys.nmodes) for i ∈ 1:2]
    κₗᵣ = SparseArray{type}(undef, prdims, prdims, sys.nmodes, sys.nmodes)

    return TDVPMatrices(sys, S, H, κᵣ, κₗᵣ, L, Ldag, LdagL, a, aH, aL, aLdagL)
end



function TDVPFunction(sys::TDVPSystem)
    H = sys.hamiltonian
    J = sys.lossoperators
    ord = sys.ord
    nmodes = sys.nmodes

    resh = reshape_basis
    to_b = tobasisaugmented
    build_func_kwargs = Dict(:expression => false,
                            :lineumbers =>false,
                            :skipzeros => true,
                            :cse => true)
                           # force_SA => true

    a_ops = [QuantumCumulants.Destroy(FockSpace(:bosonic_mode), Symbol("a",aon), aon) for aon in 1:nmodes];

    #strangely, doing it with the macro throws an error...
    α = [SymbolicUtils.Sym{Number, Nothing}(Symbol("α",i), nothing) for i ∈ 1:nmodes]


    @info("Compiling runtime-generated functions. This might take a while...")
    S_f    = build_function(resh(to_b(1, α; ord = ord), nmodes), α; build_func_kwargs...)[2]
    H_f    = build_function(resh(to_b(H, α; ord = ord), nmodes), α; build_func_kwargs...)[2]

    κᵣ_f   = build_function(resh(right_derivative_matrix(α; ord = ord), nmodes), α; build_func_kwargs...)[2]
    κₗᵣ_f  = build_function(resh(leftright_derivative_matrix(α; ord = ord), nmodes), α; build_func_kwargs...)[2]

    L_f    = build_function(resh(to_b(J, α; ord = ord), nmodes), α; build_func_kwargs...)[2] #special ndims = 3
    Ldag_f = build_function(resh(to_b(adjoint.(J), α; ord = ord), nmodes), α; build_func_kwargs...)[2] #special ndims =3
    LdagL_f= build_function(resh(to_b(adjoint.(J) .* J, α; ord = ord), nmodes), α; build_func_kwargs...)[2] #special ndims =3

    a_f    = build_function(resh(to_b(a_ops, α; ord = ord), nmodes), α; build_func_kwargs...)[2] #special ndims = 3
    aH_f   = build_function(resh(to_b(map(x -> simplify(x*H), a_ops), α; ord = ord), nmodes), α; build_func_kwargs...)[2] #special ndims = 3

    aL_f   = build_function(resh(to_b([simplify(a*L) for L ∈ J, a ∈ a_ops], α; ord = ord), nmodes), α; build_func_kwargs...)[2] #special ndims = 4
    aLdagL_f =build_function(resh(to_b([simplify(a*L'*L) for L ∈ J, a ∈ a_ops], α; ord = ord), nmodes), α; build_func_kwargs...)[2] #special ndims = 4
    @info("Done")
    return TDVPFunction(sys,S_f,H_f, κᵣ_f, κₗᵣ_f, L_f, Ldag_f, LdagL_f, a_f, aH_f, aL_f, aLdagL_f)
end
